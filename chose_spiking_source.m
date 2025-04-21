function spike_source = chose_spiking_source(source_positions, sources_in_the_target_region, ...
    spread, speak)
    % INPUTS
    %   source_positions = positions of the sources in the ROIs or the source space 
    %       (nb_of_sources x 3).
    %   sources_in_the_target_region = subset of the source space
    %       corresponding for example to one ROI. If defined, the spike will 
    %       be constrained to this region.
    %   spread = spatial spread of the spike, i.e., radius around the spike 
    %       centre that constrain the spike presence (cm).
    %   speak = whether to let the function speak or not.
    % OUTPUT
    %   spike_source = ID of the source to which the spike centre will be 
    %       assigned.
    
    if speak
        fprintf('\nChosing the spike centre...')
    end
    
    if isempty(sources_in_the_target_region) % the target region is the source space input
        sources_in_the_target_region = (1:size(source_positions, 1))';
    end
    
    nb_of_potential_sources = length(sources_in_the_target_region);
   
    % find the spiking source that have the most neighbours inside the
    % target region
    positions_in_the_target_region = source_positions(sources_in_the_target_region,:);
    potential_spike_sources = [];
    nb_of_neighbours_inside_max = 0;

    for i=1:nb_of_potential_sources
        % find all neighbours in the spread perimeter
        distance_to_source = sqrt(sum((source_positions-positions_in_the_target_region(i,:)).^2, 2));
        % all neighbours
        neighbours = find(distance_to_source<spread);
        % neighbours outside the ROI
        nb_of_neighbours_inside = length(find(ismember(neighbours, sources_in_the_target_region)));
        if nb_of_neighbours_inside > nb_of_neighbours_inside_max
            potential_spike_sources = sources_in_the_target_region(i);
            nb_of_neighbours_inside_max = nb_of_neighbours_inside;
        elseif nb_of_neighbours_inside==nb_of_neighbours_inside_max
            potential_spike_sources(end+1) = sources_in_the_target_region(i);
        end
    end
    
    spike_source = potential_spike_sources(1); 
    
    if speak
        fprintf(['\nFound a spike source with ' num2str(nb_of_neighbours_inside_max) ' sources inside the ROI.\n'])
    end
end