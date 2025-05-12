function [source_baseline_activity, spike_trace, source_baseline_activity_before_filtering, ...
    spike_start, spike_end] = simulate_spikes(varargin)
    % INPUTS
    % - optional -
    %   nb_of_sources = number of sources (default = 10).
    %   nb_of_orientations = number of orientations (default = 3).
    %   T = duration of the source time-series (s) (default = 4s).
    %   Fs = sampling frequency in (Hz) (default = 250Hz).
    %   spike_duration = duration of the spike (s) (default = 0.1s).
    %   single_time_series_SNR = SNR wanted in only one time-series (spike 
    %       vs 1 baseline). It does not account for the number of sources
    %       that are spiking, i.e., the spatial extend of the spike 
    %       (default = []).
    %   a = factor with which the spike will be multiplied (default = []).
    %   noise_colour = colour of the noise ('white' or 'pink') (default = 
    %   'pink').
    %   filter_coeffs = coefficients for the filter. filter_coeffs = [A,
    %       B] (default filter built with order = 4 and highpass frequency 
    %       = 1Hz).
    %   speak = whether to let the function talks or not (default = true).
    %
    % OUTPUTS
    %   source_baseline_activity = set of time-series containing the baseline
    %       activity (nb_of_sources*nb_of_orientations x T*Fs).
    %   spike_trace = one time-series containing the spike pattern (1 x
    %       T*Fs).
    %   source_baseline_activity_before_filtering = source baseline
    %       activity before the HP filter was applied.
    %   spike_start = data point when the spike starts.
    %   spike_end = data point when the spike ends.
    
    % get the optional arguments, and using default values if not included
    % (based on Fieltrip function)
    nb_of_sources = ft_getopt(varargin, 'nb_of_sources', 10);
    nb_of_orientations = ft_getopt(varargin, 'nb_of_orientations', 3);
    T = ft_getopt(varargin, 'T', 4);
    Fs = ft_getopt(varargin, 'Fs', 250);
    spike_duration = ft_getopt(varargin, 'spike_duration', 0.1);
    single_time_series_SNR = ft_getopt(varargin, 'single_time_series_SNR', []);
    a = ft_getopt(varargin, 'a', []);
    noise_colour = ft_getopt(varargin, 'noise_colour', 'pink');
    filter_coeffs = ft_getopt(varargin, 'filter_coeffs', []);
    speak = ft_getopt(varargin, 'speak', true);
    
    % check inputs
    if spike_duration>T
        error('The spike duration cannot be longer than the time-series duration T.')
    end
    if ~isempty(single_time_series_SNR) && ~isempty(a)
        error('Please define either the SNR at the source level or the spike amplitude factor a, but not both.')
    end
    if isempty(filter_coeffs)
        [B, A] = butter(4, 1/(Fs/2), 'high');
        filter_coeffs = [B; A];
    end
    
    if speak
        fprintf('\nCreating baseline activity and spike trace...')
    end

    % baseline activity
    try
        cn = dsp.ColoredNoise('Color', noise_colour, 'SamplesPerFrame', T*Fs, ...
            'NumChannels', nb_of_sources*nb_of_orientations);
        source_baseline_activity = cn()';
        % detrending
        source_baseline_activity = detrend(source_baseline_activity', 0)'; 
        source_baseline_activity = detrend(source_baseline_activity', 1)';
    catch % if error due to license
        if speak
            warning('No license available for signal_blocks. Using random white noise instead.')
        end
        source_baseline_activity = randn(nb_of_sources*nb_of_orientations, T*Fs);
    end
    
    % baseline filtering
    source_baseline_activity_before_filtering = source_baseline_activity;
    if ~isempty(filter_coeffs)
        B = filter_coeffs(1,:);
        A = filter_coeffs(2,:);
        for i=1:size(source_baseline_activity, 1)
            source_baseline_activity(i,:) = filtfilt(B, A, source_baseline_activity(i,:));
        end
    end
    
    % common pattern (spike)
    one_third = fix(spike_duration*Fs/3);
    spike_pattern = [hann(2*one_third)', zeros(1,one_third)];
    spike_pattern = spike_pattern - 0.5*flip(spike_pattern);
    spike_pattern = spike_pattern/max(abs(spike_pattern));
    spike_start = fix(T*Fs/2-spike_duration*Fs/2)+1; % spike start sample in the main vector
    spike_end = fix(T*Fs/2-spike_duration*Fs/2)+length(spike_pattern); % spike end sample in the main vector
    
    % spike amplitude
    % powers in the baseline and in the spike
    baseline_pow = mean(source_baseline_activity(1,:).^2); % first time-series is used as baseline
    spike_pow = mean(spike_pattern.^2);
    % determining the spike amplitude
    if isempty(a)
        a = 1; % initialisation
        if ~isempty(single_time_series_SNR)
            a = sqrt(single_time_series_SNR*(baseline_pow/spike_pow));
            if speak
                fprintf(['\nConsidering the SNR (' num2str(single_time_series_SNR) ') at the source level. a = ' ...
                    num2str(a) '.'])
            end
        end
    elseif speak
        single_time_series_SNR = a^2*(spike_pow/baseline_pow);
        fprintf(['\nThe amplitude a is fixed (' num2str(a) '). Estimated source SNR = ' ...
            num2str(single_time_series_SNR) '.'])
    end
    
    spike_trace = zeros(1, T*Fs);
    spike_trace(spike_start:spike_end) = spike_pattern;
    spike_trace = spike_trace.*a;
    
    if speak
        fprintf('\n')
    end
end