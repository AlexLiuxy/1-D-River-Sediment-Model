function Forcing = Build_Forcing_Constant(Config, Grid)

nt_dummy = 1;

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

if isfield(Config, 'F_OM_total') && ~isempty(Config.F_OM_total)
    F_OM_total = Config.F_OM_total;
else
    F_OM_total = Config.BE * Config.NPP * 1e-4;
end

Forcing.F_OM_total = F_OM_total;
Forcing.F_lab_OM = Config.f_lab * F_OM_total;
Forcing.F_ref_OM = (1 - Config.f_lab) * F_OM_total;

Forcing.time = 0;
Forcing.nt_dummy = nt_dummy; %#ok<NASGU>
end