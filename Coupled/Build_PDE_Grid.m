function Grid = Build_PDE_Grid(Config)

Grid.z = linspace(0, Config.Lbottom, Config.n_pde).';
Grid.n = numel(Grid.z);
Grid.dz = Grid.z(2) - Grid.z(1);

Grid.poros = Config.porosbottom + ...
    (Config.porostop - Config.porosbottom) .* exp(-Grid.z ./ Config.porosscale);

Grid.D_solid_mix = Config.Bioturbbottom + ...
    (Config.Bioturbtop - Config.Bioturbbottom) .* exp(-Grid.z ./ Config.bioturbscale);

Grid.Alpha_exchange = Config.Bioirrig_bottom + ...
    (Config.Bioirrig_top - Config.Bioirrig_bottom) .* exp(-Grid.z ./ Config.Bioirrig_scale);

Grid.v_burial = Config.vbottom .* ...
    (1 - Config.porosbottom) ./ max(1 - Grid.poros, 1e-6);

Grid.v_burial_Fluid = Config.vbottom_fluid .* ...
    (1 + Config.porosbottom) ./ (1 + Grid.poros);

age = Config.ageinit + cumsum(Grid.dz ./ max(Grid.v_burial, 1e-12));
Grid.k_sed = 10.^(-0.95 .* log10(age) - 0.81);

end