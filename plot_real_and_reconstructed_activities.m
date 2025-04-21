function plot_real_and_reconstructed_activities(spike_source, spiking_source_id, ...
	nb_of_sources, T, Fs, real_source_activity, eeg_activity, reconstructed_source_activity, ...
    spike_borders, save_data, saving_prefix, folder2save, format)

    time = (0:1/Fs:T-1/Fs)-T/2;
    nb2plot = 20;
    fontsize = 20;

    spike_start = spike_borders(1);
    spike_end = spike_borders(2);
    nb_of_channels = size(eeg_activity, 1);
    
    % chose the scalp channels to plot
    eeg_activity_pow = mean(eeg_activity(:,spike_start:spike_end).^2, 2);
    spk_ch = find(eeg_activity_pow==max(eeg_activity_pow)); % find the channel where the best spike is
    ch2plot = [spk_ch-3*nb2plot/2, spk_ch+3*nb2plot/2];
    ch2plot = ch2plot + (ch2plot(1)<0)*(abs(ch2plot(1))) + 1; % avoid a value below 1
    ch2plot = ch2plot - (ch2plot(2)>nb_of_channels)*(ch2plot(2)-nb_of_channels);
    
    % chose the sources to plot
    if strcmp(spike_source, 'none') || strcmp(spike_source, 'all')
        s2plot = [1, 20]; % interval for sources
    elseif isnumeric(spike_source) || strcmp(spike_source, 'random')
        s2plot = round([spiking_source_id*3-3*nb2plot/2, spiking_source_id*3+2+3*nb2plot/2]); % interval
        s2plot = s2plot + (s2plot(1)<0)*(abs(s2plot(1))) + 1; % avoid a value below 1
        s2plot = s2plot - (s2plot(2)>nb_of_sources*3)*(s2plot(2)-nb_of_sources*3); % avoid a value above the total number of signals
    else
        error('Incorrect definition of the variable spiking_source')
    end
    
    figure('Units', 'normalized', 'Position', [0 0.04 1 .88]);
    plot_title = 'Real and reconstructed activities';
    sgtitle(plot_title, 'fontsize', fontsize, 'fontweight', 'bold')
    subplot(1,3,1); hold on; title('Real source time-series')
    gap = max(real_source_activity, [], [1 2]);
    for i=s2plot(1):s2plot(2) 
        plot(time, real_source_activity(i,:)-gap*(i-1), 'k', 'Linewidth', .1)
    end
    xlabel('Time [s]'); ylabel('Sources'); yticks([]); set(gca, 'fontsize', fontsize)
    ylim([min(real_source_activity(s2plot(2),:)-gap*(s2plot(2)-1)), ...
        max(real_source_activity(s2plot(1),:)-gap*(s2plot(1)-1))]);
    subplot(1,3,2); hold on; title('Projected scalp activity')
    gap = 2*max(eeg_activity(ch2plot(1):ch2plot(2),:), [], [1 2]);
    for i=ch2plot(1):ch2plot(2) 
        plot(time, eeg_activity(i,:)-gap*(i-1), 'k', 'Linewidth', .1)
    end
    xlabel('Time [s]'); ylabel('Sensors'); yticks([]); set(gca, 'fontsize', fontsize)
    ylim([min(eeg_activity(ch2plot(2),:)-gap*(ch2plot(2)-1)), max(eeg_activity(ch2plot(1),:)-gap*(ch2plot(1)-1))]);
    subplot(1,3,3); hold on; title('Reconstructed source time-series')
    gap = max(reconstructed_source_activity, [], [1 2]);
    for i=s2plot(1):s2plot(2) 
        plot(time, reconstructed_source_activity(i,:)-gap*(i-1), 'k', 'Linewidth', .1)
    end
    xlabel('Time [s]'); ylabel('Sources'); yticks([]); set(gca, 'fontsize', fontsize)
    ylim([min(reconstructed_source_activity(s2plot(2),:)-gap*(s2plot(2)-1)), ...
        max(reconstructed_source_activity(s2plot(1),:)-gap*(s2plot(1)-1))]);
    
    if save_data
        if strcmp(format, 'all')
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.png']))
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.svg']))
        else
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.' format]))
        end
    end
    
end