function [source_activity, spiking_source_id, a, spike_polarity, weights] = assign_spike_to_baseline_activity(source_baseline_activity, ...
    spike_trace, source_positions, varargin)
    % INPUTS
    % - required - 
    %   source_baseline_activity = set of time-series containing the baseline
    %       activity (nb_of_sources*nb_of_orientations x T*Fs). Output of
    %       simulate_spikes().
    %   spike_trace = one time-series containing the spike pattern (1 x
    %       T*Fs). Output of simulate_spikes().
    %   source_positions = positions of the sources in the ROIs or the source space 
    %       (nb_of_sources x 3).
    % - optional -
    %   nb_of_sources = number of sources (default = number of time-series
    %       /nb_of_orientations).
    %   nb_of_orientations = number of orientations (default = 3).
    %   spike_source = how to assign the spike in the sources. 'none' = no
    %       spike, 'all' = in all source, 'random' = the spike centre is 
    %       randomly chosen among the sources. spiking_source can also be 
    %       an integer, corresponding to the spike centre. (Default =
    %       'random').
    %   spike_ori = orientation of the spike. It should be 'random' (meaning the
    %       orientation will be random) or between [1, nb_of_orientations].
    %       (Default = 'ramdom').
    %	weighting_func = spike weighting function, i.e., how to
    %       weight the spike amplitude with distance from the spike centre.
    %       (Default = semi-hann window)
    %   spread = spatial spread of the spike, i.e., radius around the spike 
    %       centre that constrain the spike presence (cm) (Default = 15cm).
    %   sources_in_the_target_region = subset of the source space
    %       corresponding for example to one ROI. If defined, the spike will 
    %       be constrained to this region.
    %   source_SNR = SNR desired at the source level. 
    %   scalp_SNR = SNR desired at the scalp level. The leadfield needs to
    %       be input as well.
    %   L = leadfield matrix (nb_of_sensors x
    %       nb_of_sources*nb_of_orientations).
    %   speak = whether to let the function speak or not.
    %
    % OUTPUTS
    %   source_activity = set of time-series containing the baseline
    %       activity + the spike pattern (nb_of_sources*nb_of_orientations 
    %       x T*Fs).
    %   spiking_source_id = ID of the spike centre.
    %   a = factor with wich the spike amplitude was multiplied, either 1
    %       or calculated based on the desired scalp SNR.
    %   spike_polarity = direction of the spike along the orientation (1 or 
    %       -1).
    %   weights = weights of the spike amplitude in each source, between
    %       [0, 1] (nb_of_sources x 1).
    
    % get the optional arguments, and using default values if not included
    % (based on Fieltrip function)
    nb_of_orientations = ft_getopt(varargin, 'nb_of_orientations', 3);
    nb_of_sources = ft_getopt(varargin, 'nb_of_sources', size(source_baseline_activity, 1)/nb_of_orientations);
    spike_source = ft_getopt(varargin, 'spike_source', 'random');
    spike_ori = ft_getopt(varargin, 'spike_ori', 'random');
    spread = ft_getopt(varargin, 'spread', 15);
    weighting_func = ft_getopt(varargin, 'spike_ori', @(dist) 0.5*(1-cos(2*pi.*dist./(2*spread)+pi)));
    sources_in_the_target_region = ft_getopt(varargin, 'sources_in_the_target_region', []);
    source_SNR = ft_getopt(varargin, 'source_SNR', []);
    scalp_SNR = ft_getopt(varargin, 'scalp_SNR', []);
    L = ft_getopt(varargin, 'L', []);
    speak = ft_getopt(varargin, 'speak', 'true');
    
    % check inputs
    if size(source_baseline_activity, 1) ~= nb_of_sources*nb_of_orientations
        error(['source_baseline_activity has ' num2str(size(source_baseline_activity, 1)) ...
            ' sources, which does not correspond to nb_of_sources x nb_of_orientations = ' ...
            num2str(nb_of_sources*nb_of_orientations)])
    elseif isinteger(spike_ori) && spike_ori>nb_of_orientations
        error('The spike orientation (ori) cannot be greater than the number of orientations.')
    elseif size(source_positions, 1) ~= nb_of_sources
        error('Mismatch with the number of sources in the positions.')
    elseif ~isempty(source_SNR) && ~isempty(scalp_SNR)
        error('Please specify either source_SNR or scalp_SNR but not both')
    elseif ~isempty(L) && ~isempty(scalp_SNR) && size(L, 2) ~= nb_of_sources*nb_of_orientations
        error('Mismatch with the number of sources in the leadfield matrix L.')   
    elseif ~isempty(scalp_SNR) && isempty(L)
        warning('The SNR at the scalp level is defined but the leadfield matrix is missing. The scalp SNR will be ignored.')
    elseif isempty(scalp_SNR) && ~isempty(L)
        warning('This is useless to define the leadfield matrix L without the desired scalp SNR. L will be ignored.')   
    end
    
    if speak
        fprintf('\nAssigning spike to baseline activity...')
    end
    
    spike_ids = find(spike_trace);
    spike_start = spike_ids(1);
    spike_end = spike_ids(end);

    % initialise the spiking source id and the weights in case the spike is
    % not assigned to a specific source (i.e., in case spiking_source = 'none' 
    % or 'all')
    spiking_source_id = NaN; 
    weights = NaN; 
    
    % add the spike in the chosen sources
    if strcmp(spike_source, 'none') % the spike is not assigned to any source
        source_activity = source_baseline_activity;
    elseif strcmp(spike_source, 'all') % the spike is assigned to all sources
        source_activity = source_baseline_activity + spike_trace;
    elseif isnumeric(spike_source) || strcmp(spike_source, 'random') 
    % the spike is assigned to one source, and spread in the neighbourhood
        if isnumeric(spike_source)
            spiking_source_id = spike_source;
        else % randomly selecting the main spiking source and the orientation
            spiking_source_id = randi(nb_of_sources); 
        end
        % chose the orientation of the spike
        ori_weights = zeros(nb_of_orientations, 1); 
        if strcmp(spike_ori, 'random') % random orientation
%             ori = rand(1,3)-0.5; ori_weights = ori./norm(ori); % attention, difficult to handle if you need a ground truth
            ori_weights(randi(nb_of_orientations)) = 1; % x, y, or z
        elseif ismember(spike_ori, 1:nb_of_orientations)
            ori_weights(spike_ori) = 1; % x, y, or z
        else
            error('Wrong definition of the orientation, It should be an integer between 1 and the number of orientations.')
        end
        spike_polarity = sign(rand-0.5);
        spike_trace = spike_polarity.*spike_trace; % positive or negative spike
        % using the distance from the spiking source to weight the spike's amplitude
        dist = sqrt(sum((source_positions-source_positions(spiking_source_id,:)).^2, 2)); % distances from spiking source
        weights = weighting_func(dist);
        % the spread corresponds to half of the window
        weights(dist>spread) = 0; % set what is outside the window to 0 
        % (otherwise it will oscillate since it's a cosine)
        
        % if a target region is specified (sources_in_the_target_region),
        % it implies that we don't want the spike outside this region, so
        % we force the weights outside this region to be null. If the input
        % positions are already correponding to only one ROI, this is 
        % obviously useless to specify the ROI. 
        if speak
            fprintf(['\nNumber of spiking sources = ' int2str(length(find(weights))) '.'])
        end
        if ~isempty(sources_in_the_target_region) % if there is a target ROI
            weights(setdiff(1:nb_of_sources, sources_in_the_target_region)) = 0;
            if speak
                fprintf([' Number of spiking sources after being constrained in the ROI = ' int2str(length(find(weights))) '.'])
            end
        elseif speak
            fprintf(' Not constraining the spatial extend of the spike.')
        end
        
        % weight the spike and create the source spiking activities
        weights_all = repelem(weights, nb_of_orientations);
        ori_weights_all = repmat(ori_weights, nb_of_sources, 1);
        spike_activity = ori_weights_all.*weights_all.*spike_trace;
        
        % chose the amplitude of the spike pattern according to the scalp SNR
        a = 1; % initialising
        if ~isempty(source_SNR)
            baseline_pow = mean(source_baseline_activity(:,spike_start:spike_end).^2, [1 2]);
            spike_pow = mean(spike_activity(:,spike_start:spike_end).^2, [1 2]);
            a = sqrt(source_SNR*(baseline_pow/spike_pow));
            if speak
                fprintf(['\nConsidering the SNR at the source level. a = ' num2str(a) '.'])
            end
        elseif ~isempty(scalp_SNR) && ~isempty(L)
            % project baseline and spike during the spike
            baseline_projection = L*source_baseline_activity(:,spike_start:spike_end); 
            spike_projection =  L*spike_activity(:,spike_start:spike_end); 
            % calculate SNR on this channel
            baseline_pow = mean(baseline_projection.^2, [1 2]);
            spike_pow = mean(spike_projection.^2, [1 2]);
            a = sqrt(scalp_SNR*(baseline_pow/spike_pow));
            if speak
                fprintf(['\nConsidering the SNR at the scalp level. a = ' num2str(a) '.'])
            end
        end
        
        % multiplying the spike trace by a
        source_activity = source_baseline_activity + a.*spike_activity;
        
    else 
        error('Incorrect definition of the variable spiking_source. Is should be none, all, random, or an integer.')
    end   
    
    if speak
        fprintf('\n')
    end
end