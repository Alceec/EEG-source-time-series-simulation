function plot_distribution_on_individual_mri(source2plot, mri, parameter, interpmethod, method)
    % INPUTS
    %   source2plot = Fieldtrip structure containing the positions (pos
    %       nb_tot_of_sources x 3), the source space (inside nb_tot_of_sources
    %       x 1), the value to plot (e.g., pow nb_tot_of_sources x 1).
    %   mri = MRI scan read with ft_read_mri().
    %   parameter = parameter to plot (e.g., 'pow')
    %   interpmethod = interpolation method (e.g., 'spline' or 'nearest')
    %   method = how to plot (e.g., 'ortho' or 'slice')

    cfg = [];    
    cfg.parameter = parameter;
    cfg.interpmethod = interpmethod;
    source_int = ft_sourceinterpolate(cfg, source2plot, mri);
    cfg = [];
    cfg.method  = method; 
    cfg.funparameter = parameter;
    ft_sourceplot(cfg, source_int);
    
end