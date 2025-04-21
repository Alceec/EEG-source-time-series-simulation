function SNR = get_SNR(time_series, spike_borders)
    % INPUTS
    %   time_series = n x t matrix, where n = number of time_series
    %       (nb_of_sensors or nb_of_sources*nb_of_orientations) and t = 
    %       number of time points.
    %   spike_borders = [spike_start spike_end], time samples of the start
    %       and end of the spike.
    % OUTPUT
    %   SNR = signal-to-noise ratio.

    spike_start = spike_borders(1);
    spike_end = spike_borders(2);
    nb_of_time_points = size(time_series, 2);
    
    pow_baseline = mean(time_series(:,[1:spike_start-1 spike_end+1:nb_of_time_points]).^2, [1 2]); % mean along time
    pow_spike = mean(time_series(:,spike_start:spike_end).^2, [1 2]); % mean along time
    SNR = pow_spike/pow_baseline;
    
end