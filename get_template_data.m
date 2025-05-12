function [mri, leadfield, inverse_filters_list, source_positions, label_vector_inside, ...
    atlas_labels_reduced] = get_template_data(bids_dir_path, meth, lambda)

    if isempty(meth)
        meth = 'eloreta';
    end
    
    if isempty(lambda)
        lambda = {0.05};
    end

    % MRI (optional, for visualisation)
    mri_path = dir(fullfile(bids_dir_path, 'sub-template', 'ses-preop', 'anat', ...
        'sub-template_ses-preop_T1w.nii.gz'));
    mri = ft_read_mri(fullfile(mri_path.folder, mri_path.name));
    mri.coordsys = 'ras';

    % leadfield matrix (optional, only if you want to specify a SNR at the
    % scalp level or/and use it for testing inverse solutions)
    leadfield_path = fullfile(bids_dir_path, 'derivatives', 'source_modelling_v3.2.0', ...
            'sub-template', 'ses-preop', 'eeg', 'sub-template_ses-preop_leadfield.mat');
    if exist(leadfield_path, 'file')
        load(leadfield_path)
    else % create leadfield
        % headmodel
        headmodel_path = dir(fullfile(bids_dir_path, 'derivatives', 'mri_preprocessing_v3.2.0', ...
            'sub-template', 'ses-preop', 'anat', '*headmodel.mat' ));
        load(fullfile(headmodel_path.folder, headmodel_path.name), 'headmodel')
        % source model
        sourcemodel_path = dir(fullfile(bids_dir_path, 'derivatives', 'mri_preprocessing_v3.2.0', ...
            'sub-template', 'ses-preop', 'anat', '*sourcemodel.mat' ));
        load(fullfile(sourcemodel_path.folder, sourcemodel_path.name), 'sourcemodel')
        % sensor
        electrodes_tsv = ft_read_tsv(fullfile(bids_dir_path, 'sub-template', 'ses-preop', ...
            'eeg', 'sub-template_ses-preop_electrodes.tsv'));
        elec         = [];
        elec.label   = electrodes_tsv.name;
        elec.elecpos = [electrodes_tsv.x electrodes_tsv.y electrodes_tsv.z];
        % create leadfield
        cfg         = [];
        cfg.elec    = elec;  
        cfg.grid    = sourcemodel;   
        cfg.headmodel = headmodel;   
        leadfield   = ft_prepare_leadfield(cfg);
        save(leadfield_path, 'leadfield')
    end

    % inverse solution (optional, only if you want to test inverse solution 
    % and you're using inverse solution that does not depend on the data 
    % such as eLORETA)
    inverse_filters_list = cell(length(lambda), 1);
    if ismember('eloreta', meth)
        for i = 1:length(lambda)
            inverse_filters_path = fullfile(bids_dir_path, 'derivatives', 'source_modelling_v3.2.0', ...
                'sub-template', 'ses-preop', 'eeg', ['sub-template_ses-preop_desc-eloreta_lambda-' ...
                num2str(lambda{i}) '_inversefilters.mat']);
            if exist(inverse_filters_path, 'file')
                load(inverse_filters_path)
            else % creating inverse operator
                C = eye(length(leadfield.label));
                inverse_filters = ft_inverse_eloreta(leadfield, [], [], [], C, ...
                    'keepfilter', 'yes', 'lambda', lambda{i});
                save(inverse_filters_path, 'inverse_filters')
            end
            inverse_filters_list{i} = inverse_filters;
        end
    end

    % atlas
    atlas_path = fullfile(bids_dir_path, 'derivatives', 'cmp-v3.2.0', 'sub-template', ...
                'ses-preop', 'anat', 'sub-template_ses-preop_atlas-L2018_res-scale2_dseg.nii.gz');
    atlas_labels = readtable(fullfile(bids_dir_path, 'derivatives', 'cmp-v3.2.0', ...
        'sub-template', 'ses-preop', 'anat', 'sub-template_ses-preop_atlas-L2018_res-scale2_dseg.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    [atlas_path, temp_folder] = decompress_mri(atlas_path);
    atlas = ft_read_mri(atlas_path);
    atlas.coordsys = 'ras'; 
    if exist(temp_folder, 'dir')
        delete(atlas_path)
        rmdir(temp_folder)
        clear temp_folder
    end
    ROIidx2remove = [58:62 122:126 129];
    atlas.anatomy(ismember(atlas.anatomy, ROIidx2remove)) = 0;
    atlas_labels_reduced = atlas_labels; 
    atlas_labels_reduced(ROIidx2remove, :) = []; % labels and id
    % overlap the atlas and the source points
    [label_vector, ~]= assigned_label2source_points(leadfield, atlas, atlas_labels);
    label_vector_inside = label_vector(leadfield.inside);
    source_positions = leadfield.pos(leadfield.inside, :); % source positions
    
end