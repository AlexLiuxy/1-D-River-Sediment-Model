function State = Initialize_PDE_State_From_Steady(Steady, Grid)

z0 = Steady.z_sed(:);
z  = Grid.z(:);

ip = @(v, floorv) max(interp1(z0, v(:), z, 'linear', 'extrap'), floorv);

State.OM_lab = ip(Steady.OM_lab, 1e-12);
State.OM_ref = ip(Steady.OM_ref, 0);

State.O2    = ip(Steady.O2, 0);
State.Fe2   = ip(Steady.Fe2, 0);
State.FeOOH = ip(Steady.FeOOH, 0);
State.SO4   = ip(Steady.SO4, 0);
State.HS    = ip(Steady.HS, 0);
State.CH4   = ip(Steady.CH4, 0);
State.DIC   = ip(Steady.DIC, 1e-12);
State.ALK   = ip(Steady.ALK, 1e-12);
State.CaCO3 = ip(Steady.CaCO3, 0);

end