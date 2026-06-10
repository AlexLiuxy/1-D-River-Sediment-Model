function Grid = Build_Grid_1D(Config, Params)
% BUILD_GRID_1D
% Build the 1-D sediment grid and depth-dependent physical fields.

Grid.z = linspace(0, Config.Lbottom, Config.n).';
Grid.n = numel(Grid.z);
Grid.dz = Grid.z(2) - Grid.z(1);

z = Grid.z;

Grid.poros = Config.porosbottom + ...
    (Config.porostop - Config.porosbottom) .* exp(-z ./ Config.porosscale);

Grid.D_solid_mix = Config.Bioturbbottom + ...
    (Config.Bioturbtop - Config.Bioturbbottom) .* exp(-z ./ Config.bioturbscale);

Grid.Alpha_exchange = Config.Bioirrig_bottom + ...
    (Config.Bioirrig_top - Config.Bioirrig_bottom) .* exp(-z ./ Config.Bioirrig_scale);

Grid.v_solid = Config.vbottom .* ...
    (1 - Config.porosbottom) ./ max(1 - Grid.poros, 1e-6);

Grid.v_fluid = Config.vbottom_fluid .* ...
    (1 + Config.porosbottom) ./ max(1 + Grid.poros, 1e-6);

age = Config.ageinit + cumsum(Grid.dz ./ max(Grid.v_solid, 1e-12));
Grid.k_sed = 10.^(-0.95 .* log10(age) - 0.81);
Grid.k_sed = Config.k_sed_scale .* Grid.k_sed;

Grid.rho = Params.rho;
end