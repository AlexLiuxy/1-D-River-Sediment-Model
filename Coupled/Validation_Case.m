function [Config, Params, Case] = Validation_Case(name)

Config = Config_Baseline();
Params = Params_Static();

Case = struct();
Case.name = name;

switch lower(name)

    case 'geneva_1'
        Case.note = 'Lake Geneva site #1, Rhone delta. Methane-rich delta sediment.';

        Config.Lbottom = 21;

        Config.T_future = 4;
        Config.Salinity = 0.2;

        Config.O2init = 125;
        Config.SO4init = 500;
        Config.CH4init = 0;
        Config.DICinit = 2500;
        Config.HCO3init = 4000;
        Config.Calcium = 2000;
        Config.Feinit = 0;
        Config.HSinit = 0;

        Config.vbottom = 2.3;
        Config.porostop = 0.57;
        Config.porosbottom = 0.55;
        Config.porosscale = 5;

        Config.Bioturbtop = 1;
        Config.Bioturbbottom = 0.05;
        Config.bioturbscale = 3;

        Config.Bioirrig_top = 10;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 1;

        Config.F_OM_total = 131.7 * 1e-4;
        Config.f_lab = 0.6;

        Config.F_FeOx = 0.1;
        Config.F_CaCO3 = 0;

        Config.water_depth_m = 130;
        Config.use_CH4_bubbling = true;

        Params.DO2 = 357;
        Params.DSO4 = 164;
        Params.DHCO3 = 161;
        Params.DCa = Params.DHCO3;
        Params.DCH4 = 274;


    case 'geneva_4'
        Case.note = 'Lake Geneva site #4, profundal sediment. Lower OC accumulation than delta sites.';

        Config.Lbottom = 27;

        Config.T_future = 4;
        Config.Salinity = 0.2;

        Config.O2init = 125;
        Config.SO4init = 500;
        Config.CH4init = 0;
        Config.DICinit = 1500;
        Config.HCO3init = 2500;
        Config.Calcium = 2000;
        Config.Feinit = 0;
        Config.HSinit = 0;

        Config.vbottom = 0.5;
        Config.porostop = 0.75;
        Config.porosbottom = 0.72;
        Config.porosscale = 5;

        Config.Bioturbtop = 1;
        Config.Bioturbbottom = 0.05;
        Config.bioturbscale = 3;

        Config.Bioirrig_top = 10;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 1;

        Config.F_OM_total = 31.2 * 1e-4;
        Config.f_lab = 0.6;

        Config.F_FeOx = 0.1;
        Config.F_CaCO3 = 0;

        Config.water_depth_m = 210;
        Config.use_CH4_bubbling = true;

        Params.DO2 = 357;
        Params.DSO4 = 164;
        Params.DHCO3 = 161;
        Params.DCa = Params.DHCO3;
        Params.DCH4 = 274;


    case 'michigan'
        Case.note = 'Lake Michigan. Redox partition and low methane validation case.';

        Config.Lbottom = 18;

        Config.T_future = 4.1;
        Config.Salinity = 0;

        Config.O2init = 315;
        Config.SO4init = 180;
        Config.CH4init = 0;
        Config.DICinit = 2500;
        Config.HCO3init = 2500;
        Config.Calcium = 1200;
        Config.Feinit = 7;
        Config.HSinit = 0;

        Config.vbottom = 0.08;
        Config.porostop = 0.92;
        Config.porosbottom = 0.87;
        Config.porosscale = 10;

        Config.Bioturbtop = 6.2;
        Config.Bioturbbottom = 0.37;
        Config.bioturbscale = 5;

        Config.Bioirrig_top = 20;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 2;

        Config.F_OM_total = 8.6 * 365 * 12e-3 * 1e-4;
        Config.f_lab = 0.6;

        Config.F_FeOx = 1.5;
        Config.F_CaCO3 = 0.5;

        Config.water_depth_m = 101;
        Config.use_CH4_bubbling = true;

        Params.rho = 2.45;
        Params.k_SO4 = 30;

    otherwise
        error('Unknown validation case: %s', name);
end
end