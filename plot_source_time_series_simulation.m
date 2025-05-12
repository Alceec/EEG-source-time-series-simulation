function plot_source_time_series_simulation(spike_source, spiking_source_id, ...
    nb_of_sources, T, Fs, source_baseline_activity_before_filtering, source_baseline_activity, ...
    source_activity, spike_trace, source_positions, sources_in_the_target_region, ...
    weighting_func, spread, weights, leadfield, mri, save_data, saving_prefix, ...
    folder2save, format)

    time = (0:1/Fs:T-1/Fs)-T/2;
    nb2plot = 8;
    fontsize = 20;
    
    if isempty(sources_in_the_target_region)
        sources_in_the_target_region = (1:size(source_positions, 1))';
    end
    nb_of_sources_in_ROI = length(sources_in_the_target_region);
        
    if strcmp(spike_source, 'none') || strcmp(spike_source, 'all')
        s2plot = [1, 20]; % interval
        colors = zeros(nb_of_sources_in_ROI, 3);
    elseif isnumeric(spike_source) || strcmp(spike_source, 'random')
        if nb_of_sources<=nb2plot/2
            s2plot = [1, nb_of_sources*3]; % interval
        else
            s2plot = [spiking_source_id*3-3*nb2plot/2, spiking_source_id*3+2+3*nb2plot/2]; % interval
            s2plot = s2plot + (s2plot(1)<0)*(abs(s2plot(1))) + 1; % avoid a value below 1
            s2plot = s2plot - (s2plot(2)>nb_of_sources*3)*(s2plot(2)-nb_of_sources*3); % avoid a value above the total number of signals
        end
        colors = zeros(nb_of_sources_in_ROI, 3); colors(:,1) = weights(sources_in_the_target_region);
    else
        error('Incorrect definition of the variable spiking_source')
    end
    figure('Units', 'normalized', 'Position', [0 0.04 1 .88]);
    plot_title = 'Simulated source time-series';
    sgtitle(plot_title, 'fontsize', fontsize, 'fontweight', 'bold')
    subplot(2,3,1); hold on; title('Baseline activity')
    gap = max(source_baseline_activity, [], [1 2]);
    for i=s2plot(1):s2plot(2) 
        plot(time, source_baseline_activity(i,:)-gap*(i-1), 'k', 'Linewidth', .1)
    end
    xlabel('Time [s]'); ylabel('Sources'); yticks([])
    set(gca,'fontsize', fontsize)
    subplot(2,3,2); hold on; title('Power spectrum of the baseline activity')
    [p, freqs] = pspectrum(source_baseline_activity_before_filtering', Fs);
    plot(freqs, p(:,s2plot(1):s2plot(2)), 'c')
    [p, freqs] = pspectrum(source_baseline_activity', Fs);
    plot(freqs, p(:,s2plot(1):s2plot(2)), 'k')
    legend('before filtering')
    xlabel('Frequency [Hz]'); ylabel('Power');
    set(gca,'fontsize', fontsize)
    subplot(2,3,3); hold on; title('Common pattern')
    plot(time, spike_trace, 'k')
    xlabel('Time [s]'); ylabel('Amplitude')
    set(gca,'fontsize', fontsize)
    subplot(2,3,4); hold on; title('Addition')
    gap = max(source_activity, [], [1 2]);
    for i=s2plot(1):s2plot(2) 
        plot(time, source_activity(i,:)-gap*(i-1), 'k', 'Linewidth', .1)
    end
    xlabel('Time [s]'); ylabel('Sources'); yticks([])
    set(gca,'fontsize', fontsize) 
    subplot(2,3,5); hold on; title('Spike taper with distance')
    plot(0:spread/100:spread, weighting_func(0:spread/100:spread), 'k')
    xlabel('Distance from the spiking source [mm]'); ylabel('Weights');
    set(gca,'fontsize', fontsize)
    % plot in 3D the ROI with the spiking source
    if nb_of_sources==length(sources_in_the_target_region) % if the source input is only the ROI
        sources2plot = source_positions;
    else % if the source input is the whole source space
        sources2plot = source_positions(sources_in_the_target_region, :);
    end
    subplot(2,3,6); axis equal; hold on
    for i=1:length(sources2plot)
        scatter3(sources2plot(i,1), sources2plot(i,2), sources2plot(i,3), 500, 'MarkerEdgeColor', ...
            [1 1 1], 'MarkerFaceColor', colors(i,:))
    end
    title('Localisation of the spiking activity')
    set(gca,'fontsize', fontsize) 
    
    if save_data
        if strcmp(format, 'all')
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.png']))
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.svg']))
        else
            saveas(gcf, fullfile(folder2save, [saving_prefix ' ' plot_title '.' format]))
        end
    end
    
    if ~isempty(leadfield) && nb_of_sources == length(find(leadfield.inside)) && ~isempty(mri)
        % plot on the MRI when dealing with the whole source space
        weights_all = zeros(length(leadfield.inside), 1); % need the entire head space
        weights_all(leadfield.inside) = weights;  
        source2plot = rmfield(leadfield, {'leadfield', 'label', 'leadfielddimord', 'cfg'});
        source2plot.pow = weights_all;
        plot_distribution_on_individual_mri(source2plot, mri, 'pow', 'nearest', 'ortho')
        saveas(gcf, fullfile(folder2save, [saving_prefix 'Localisation of the simulated source spike.png']))
        if save_data
            if strcmp(format, 'all')
                saveas(gcf, fullfile(folder2save, [saving_prefix 'Localisation of the simulated source spike.png']))
                saveas(gcf, fullfile(folder2save, [saving_prefix 'Localisation of the simulated source spike.svg']))
            else
                saveas(gcf, fullfile(folder2save, [saving_prefix 'Localisation of the simulated source spike.' format]))
            end
        end
    end
    
end