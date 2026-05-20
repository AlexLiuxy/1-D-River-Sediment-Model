function Forcing = Build_Forcing_PDE(Config)
% BUILD_FORCING_PDE
% Constant forcing for the first FV-MOL PDE benchmark.

Forcing.T = Config.T_future;
Forcing.Salinity = Config.Salinity;

Forcing.O2_top  = Config.O2init;
Forcing.SO4_top = Config.SO4init;
Forcing.DIC_top = Config.DICinit;
Forcing.ALK_top = Config.HCO3init;
Forcing.Ca_top  = Config.Calcium;
Forcing.CH4_top = Config.CH4init;
Forcing.Fe_top  = Config.Feinit;
Forcing.HS_top  = Config.HSinit;

Forcing.F_FeOx  = Config.F_FeOx;     % mmol/m2/d
Forcing.F_CaCO3 = Config.F_CaCO3;    % g/m2/yr

% Main OM forcing.
% Current stable ODE convention: BE * NPP * 1e-4 = g OM / cm2 / yr.
if isfield(Config, 'F_OM_total') && ~isempty(Config.F_OM_total)
    F_OM_total = Config.F_OM_total;
else
    F_OM_total = Config.BE * Config.NPP * 1e-4;
end

f_lab = max(0, min(1, Config.f_lab));

Forcing.F_OM_total = F_OM_total;
Forcing.F_lab_OM   = f_lab .* F_OM_total;
Forcing.F_ref_OM   = (1 - f_lab) .* F_OM_total;
end