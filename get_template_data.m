function [mri, leadfield, inverse_filters, source_positions, sources_in_the_target_region] = get_template_data(bids_dir_path, spike_region_label)

    % MRI (optional, for visualisation)
    mri_path = dir(fullfile(bids_dir_path, 'sub-template', 'anat', 'sub-tmpT1_T1w.nii'));
    mri = ft_read_mri(fullfile(mri_path.folder, mri_path.name));
    mri.coordsys = 'ras';

    % leadfield matrix (optional, only if you want to specify a SNR at the
    % scalp level or/and use it for testing inverse solutions)
    leadfield_path = dir(fullfile(bids_dir_path, 'derivatives', 'source_modelling', ...
            'sub-template', 'eeg', 'sub-template_leadfield.mat'));
    load(fullfile(leadfield_path.folder, leadfield_path.name))

    % inverse solution (optional, only if you want to test inverse solution 
    % and you're using inverse solution that does not depend on the data 
    % such as eLORETA)
    inverse_filters_path = dir(fullfile(bids_dir_path, 'derivatives', 'source_modelling', ...
        'sub-template', 'eeg', 'sub-template_desc-eloreta_inversefilters.mat'));
    load(fullfile(inverse_filters_path.folder, inverse_filters_path.name))

    % atlas
    atlas_path = fullfile(bids_dir_path, 'derivatives', 'cmp-v3.1.0', 'sub-template', ...
                'anat', 'sub-template_atlas-L2018_res-scale2_dseg.nii.gz');
    atlas_labels = readtable(fullfile(bids_dir_path, 'derivatives', 'cmp-v3.1.0', ...
        'sub-template', 'anat', 'sub-template_atlas-L2018_res-scale2_dseg.tsv'), ...
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
    source_positions = leadfield.pos(leadfield.inside, :); % source positions
    
    sources_in_the_target_region = [];
    if ~isempty(spike_region_label)
        sources_in_the_target_region = find(strcmp(label_vector(leadfield.inside), spike_region_label)); % subset of sources inside the ROI
    end
        

end