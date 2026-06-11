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
        Config.HCO3init = 3000;     % ALK, not HCO3
        Config.F_OM_total = 131.7 * 1e-4;
        Config.f_lab = 0.8;
        Config.Calcium = 2000;
        Config.Feinit = 0;
        Config.HSinit = 0;

        Config.vbottom = 2.3;
        Config.porostop = 0.57;
        Config.porosbottom = 0.55;
        Config.porosscale = 5;

        Config.Bioturbtop = 0.22;
        Config.Bioturbbottom = 0.05;
        Config.bioturbscale = 3;

        Config.Bioirrig_top = 2;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 1;

        Config.F_OM_total = 131.7 * 1e-4;
        Config.k_sed_scale = 2;
        Config.f_lab = 0.8;

        Config.F_FeOx = 0.1;
        Config.F_CaCO3 = 0;

        Config.water_depth_m = 130;
        Config.use_CH4_bubbling = true;

        Params.DO2 = 357;
        Params.DSO4 = 164;
        Params.DHCO3 = 161;
        Params.DCa = Params.DHCO3;
        Params.DCH4 = 274;
    
    case 'geneva_1_fe'
        [Config, Params, Case] = Validation_Case('geneva_1');
        Case.name = name;
    
        Config.F_FeOx = 0.5;
        Params.KFEMonod = 50;

    case 'geneva_1_trend'
        [Config, Params, Case] = Validation_Case('geneva_1');
        Case.name = name;
        Case.note = 'Lake Geneva site #1, stronger mineralization trend test.';

        Config.F_OM_total = 2.0 * 131.7 * 1e-4;
        Config.f_lab = 0.8;

        Config.Bioirrig_top = 0.5;
        Config.Bioturbtop = 0.05;

        Config.n = 101;
        Config.n_pde = Config.n;
        Config.t_spinup = 150;
        Config.dt_out_spinup = 5;

    case 'geneva_1_test'
        [Config, Params, Case] = Validation_Case('geneva_1_final');
        Case.name = name;
        Case.note = 'Lake Geneva site #1, stronger Fe availability test.';

        Config.k_sed_scale = Config.k_sed_scale * 1.2;   % 稍微增加总 mineralization
        Params.k_SO4 = 50;                                % 抑制低 SO4 下的 SRR，减少 HS
        Config.F_FeOx = 0.6;                              % 提供更多 Fe source
        Params.Fe_inventory_factor = 0.8;                 % 增加 reactive Fe pool
        Params.KFEMonod = 150;                            % 避免 Fe reduction 过强
        Params.kFeS = 1;                                  % 减弱 Fe2 被 HS 清除

     case 'geneva_1_final'
        [Config, Params, Case] = Validation_Case('geneva_1');
        Case.name = name;
        Case.note = 'Lake Geneva site 1 trend case: delta, high OC accumulation, anaerobic-dominated.';
    
        Config.k_sed_scale = 2.5;
        Config.f_lab = 0.65;
    
        Config.Bioirrig_top = 2;
        Config.Bioirrig_scale = 1.5;
    
        Params.DSO4 = 150;
        Params.k_SO4 = 15;
    
        Params.Fe_inventory_factor = 0.075;
        Params.KFEMonod = 100;
        Params.kFeS = 8;

    case 'geneva_2_proxy'
        [Config, Params, Case] = Validation_Case('geneva_1_final');
        Case.name = name;
        Case.note = 'Lake Geneva site 2 proxy based on site 1 delta parameter family.';
    
        Config.water_depth_m = 90;
        Config.Lbottom = 44;
        Config.vbottom = 2.3;
        Config.porostop = 0.53;
        Config.porosbottom = 0.51;
        Config.F_OM_total = 129.8 * 1e-4;

    case 'geneva_5_proxy'
        [Config, Params, Case] = Validation_Case('geneva_4_final');
        Case.name = name;
        Case.note = 'Lake Geneva site 5 proxy based on site 4 profundal parameter family.';
    
        Config.water_depth_m = 240;
        Config.Lbottom = 30;
        Config.vbottom = 0.4;
        Config.porostop = 0.77;
        Config.porosbottom = 0.74;
        Config.F_OM_total = 19.8 * 1e-4;

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
        Config.k_sed_scale = 10;
        Config.f_lab = 0.8;

        Config.F_FeOx = 0.1;
        Config.F_CaCO3 = 0;

        Config.water_depth_m = 210;
        Config.use_CH4_bubbling = true;

        Params.DO2 = 357;
        Params.DSO4 = 164;
        Params.DHCO3 = 161;
        Params.DCa = Params.DHCO3;
        Params.DCH4 = 274;

%         Params.KFEMonod = 200;


    case 'geneva_4_final'
        [Config, Params, Case] = Validation_Case('geneva_4');
        Case.name = name;
        Case.note = 'Lake Geneva site 4 trend case: profundal, lower OC accumulation, more aerobic than delta sites.';
    
        Config.k_sed_scale = 6;
        Config.f_lab = 0.80;
    
        Config.Bioirrig_top = 8;
        Config.Bioirrig_scale = 2;
    
        Params.DSO4 = 180;
        Params.k_SO4 = 15;
    
        Params.Fe_inventory_factor = 0.8;
        Config.F_FeOx = 0.6;
        Params.KFEMonod = 150;
        Params.kFeS = 1;

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
        Config.k_sed_scale = 10;
        Config.f_lab = 0.8;


        Config.F_FeOx = 1.5;
        Config.F_CaCO3 = 0.5;

        Config.water_depth_m = 101;
        Config.use_CH4_bubbling = true;

        Params.rho = 2.45;
        Params.k_SO4 = 30;
        Params.Fe_inventory_factor = 0.12;

%     case 'michigan_v2'
%         [Config, Params, Case] = Validation_Case('michigan');
%         Case.name = name;
%         Case.note = 'Lake Michigan, lower OM inventory and stronger sulfate supply.';
% 
%         Config.F_OM_total = 0.0007;
%         Config.k_sed_scale = 8;
%         Config.f_lab = 0.7;
% 
%         Params.DSO4 = 300;
% %         Params.KFEMonod = 3000;
%     case 'michigan_v3'
%         [Config, Params, Case] = Validation_Case('michigan');
%         Case.name = name;
%         Case.note = 'Lake Michigan, lower OM inventory with stronger reactivity.';
% 
%         Config.F_OM_total = 0.0007;
%         Config.k_sed_scale = 50;
%         Config.f_lab = 0.75;
% 
%         Params.DSO4 = 150;
%         Params.Fe_inventory_factor = 0.08;
%         Params.kFeS = 3;
%     case 'michigan_v4'
%         [Config, Params, Case] = Validation_Case('michigan');
%         Case.name = name;
%         Case.note = 'Lake Michigan, low OM stock with stronger turnover and shallower O2.';
% 
%         Config.F_OM_total = 0.0007;
%         Config.k_sed_scale = 180;
%         Config.f_lab = 0.75;
% 
%         Config.Bioirrig_top = 5;
%         Config.Bioirrig_scale = 1.5;
% 
%         Params.DO2 = 120;
%         Params.DSO4 = 150;
% 
%         Params.Fe_inventory_factor = 0.07;
%         Params.kFeS = 4;
    
    case 'michigan_v5'
        [Config, Params, Case] = Validation_Case('michigan');
        Case.name = name;
        Case.note = 'Lake Michigan, stronger anoxic demand with low OM stock.';
    
        Config.F_OM_total = 0.0012;
        Config.k_sed_scale = 80;
        Config.f_lab = 0.75;
    
        Config.Bioirrig_top = 2;
        Config.Bioirrig_scale = 1.5;
    
        Params.DO2 = 80;
        Params.DSO4 = 150;
    
        Params.Fe_inventory_factor = 0.07;
        Params.KFEMonod = 100;
        Params.kFeS = 5; 

    %  case 'michigan_v6'
    %     [Config, Params, Case] = Validation_Case('michigan_v5');
    %     Case.name = name;
    %     Case.note = 'Lake Michigan, higher FeOOH pool with similar Fe gate.';
    % 
    %     Params.DSO4 = 120;
    %     Params.KFEMonod = 700;
    % 
    %     Params.Fe_inventory_factor = 0.45;
    %     Params.kFeS = 8;      
    % 
    % case 'michigan_v7'
    %     [Config, Params, Case] = Validation_Case('michigan_v5');
    %     Case.name = name;
    %     Case.note = 'Lake Michigan, moderate FeOOH pool.';
    % 
    %     Params.DSO4 = 180;
    %     Params.KFEMonod = 1000;
    %     Params.Fe_inventory_factor = 0.15;
    %     Params.kFeS = 12;
    % 
    % case 'michigan_fe_shape'
    %     [Config, Params, Case] = Validation_Case('michigan_v7');
    %     Case.name = name;
    % 
    %     % Params.kFeOx = 50;
    %     Params.kFeS = 20;
    % 
    %     case 'michigan_v5b'
    % [Config, Params, Case] = Validation_Case('michigan_v5');
    % Case.name = name;
    % Case.note = 'Lake Michigan v5 with slightly higher FeOOH and stronger FeS sink.';
    % 
    % Params.Fe_inventory_factor = 0.09;
    % Params.KFEMonod = 200;
    % Params.kFeS = 8;
    % 
    % case 'michigan_v5c'
    % [Config, Params, Case] = Validation_Case('michigan_v5');
    % Case.name = name;
    % Case.note = 'Lake Michigan v5 with stronger FeS sink.';
    % 
    % Params.Fe_inventory_factor = 0.06;
    % Params.KFEMonod = 200;
    % Params.kFeS = 20;
    % Config.Bioirrig_top = 2;
    % Config.Bioirrig_scale = 3;
    % 
    % case 'michigan_v5_shape2'
    % [Config, Params, Case] = Validation_Case('michigan_v5');
    % Case.name = name;
    % 
    % Config.Bioirrig_top = 2;
    % Config.Bioirrig_scale = 3;
    % 
    % Params.DSO4 = 120;
    % case 'michigan_v5_shape1'
    % [Config, Params, Case] = Validation_Case('michigan_v5');
    % Case.name = name;
    % 
    % Params.k_SO4 = 15;
    % Params.kFeS = 8;

    case 'michigan_final'
        [Config, Params, Case] = Validation_Case('michigan_v5');
        Case.name = name;
        Case.note = 'Lake Michigan v5 with stronger mid-depth FeS sink and slightly higher FeOOH.';
    
        Params.k_SO4 = 10;
        Params.DSO4 = 150;
        
        % Params.Fe_inventory_factor = 0.08;
        % Config.F_FeOx = 2;
        % Params.KFEMonod = 80;
        % Params.kFeS = 10;
        Params.kFeS = 8;

        Params.Fe_inventory_factor = 0.075;
        Params.KFEMonod = 100;
    
    case 'georgia_stl2_trial'
        Case.name = name;
        Case.note = 'Georgia STL2 trial: estuarine creek-bank, sulfate-dominated, high-metabolism site.';
    
        Config.Lbottom = 50;
        Config.T_future = 25;
        Config.Salinity = 25;
    
        Config.O2init = 200;
        Config.SO4init = 22000;
        Config.CH4init = 0;
        Config.DICinit = 2500;
        Config.HCO3init = 2500;
        Config.Calcium = 9000;
        Config.Feinit = 0;
        Config.HSinit = 0;
    
        Config.vbottom = 0.10;
        Config.porostop = 0.88;
        Config.porosbottom = 0.75;
        Config.porosscale = 10;
    
        Config.Bioturbtop = 1;
        Config.Bioturbbottom = 0.05;
        Config.bioturbscale = 3;
    
        Config.Bioirrig_top = 1;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 3;
    
        Config.F_OM_total = 0.003;
        Config.k_sed_scale = 30;
        Config.f_lab = 0.75;
    
        Config.F_FeOx = 1.0;
        Config.F_CaCO3 = 0;
        Config.water_depth_m = 1;
        Config.use_CH4_bubbling = true;
    
        Params.DO2 = 200;
        Params.DSO4 = 120;
        Params.DH2S = 300;
        Params.DCH4 = 300;
        Params.DHCO3 = 250;
        Params.DCa = Params.DHCO3;
    
        Params.k_SO4 = 100;
        Params.K_CH4_SO4 = 500;
        Params.k_AOM = 3;
    
        Params.Fe_inventory_factor = 0.10;
        Params.KFEMonod = 100;
        Params.kFeS = 5;

        Config.k_sed_scale = 50;


    case 'georgia_stl2_v2'
        Case.name = name;
        Case.note = 'Georgia STL2 trial: high-turnover estuarine creek-bank sediment, sulfate-dominated.';
    
        Config.Lbottom = 50;
        Config.T_future = 25;
        Config.Salinity = 20;
    
        Config.O2init = 180;
        Config.SO4init = 22000;
        Config.CH4init = 0;
    
        Config.DICinit = 4000;
        Config.HCO3init = 3000;
        Config.Calcium = 8000;
    
        Config.Feinit = 0;
        Config.HSinit = 0;
    
        Config.vbottom = 0.10;
    
        Config.porostop = 0.88;
        Config.porosbottom = 0.75;
        Config.porosscale = 10;
    
        Config.Bioturbtop = 0.5;
        Config.Bioturbbottom = 0.02;
        Config.bioturbscale = 3;
    
        Config.Bioirrig_top = 0.2;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 5;
    
        Config.F_OM_total = 0.006;
        Config.k_sed_scale = 200;
        Config.f_lab = 0.85;
    
        Config.F_FeOx = 5.0;
        Config.F_CaCO3 = 0;
    
        Config.water_depth_m = 0.5;
        Config.use_CH4_bubbling = true;
    
        Params.DO2 = 120;
        Params.DSO4 = 40;
        Params.DH2S = 150;
        Params.DCH4 = 200;
        Params.DHCO3 = 180;
        Params.DCa = Params.DHCO3;
    
        Params.k_SO4 = 50;
        Params.K_CH4_SO4 = 300;
        Params.k_AOM = 5;
    
        Params.Fe_inventory_factor = 0.5;
        Params.KFEMonod = 80;
        Params.kFeS = 20;
        Params.kFeOx = 5;

    case 'georgia_stl2_v3'
        Case.name = name;
        Case.note = 'Georgia STL2 trial: distributed sulfate reduction and shallow Fe2+ peak.';
    
        Config.Lbottom = 50;
        Config.T_future = 25;
        Config.Salinity = 20;
    
        Config.O2init = 180;
        Config.SO4init = 22000;
        Config.CH4init = 0;
    
        Config.DICinit = 3500;
        Config.HCO3init = 4200;
        Config.Calcium = 8000;
    
        Config.Feinit = 0;
        Config.HSinit = 0;
    
        Config.vbottom = 0.35;
    
        Config.porostop = 0.88;
        Config.porosbottom = 0.75;
        Config.porosscale = 10;
    
        Config.Bioturbtop = 0.5;
        Config.Bioturbbottom = 0.02;
        Config.bioturbscale = 4;
    
        Config.Bioirrig_top = 0.8;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 4;
    
        Config.F_OM_total = 0.020;
        Config.k_sed_scale = 40;
        Config.f_lab = 0.45;
    
        Config.F_FeOx = 15.0;
        Config.F_CaCO3 = 0;
    
        Config.water_depth_m = 0.5;
        Config.use_CH4_bubbling = true;
    
        Params.DO2 = 120;
        Params.DSO4 = 80;
        Params.DH2S = 180;
        Params.DCH4 = 200;
        Params.DHCO3 = 180;
        Params.DCa = Params.DHCO3;
    
        Params.k_SO4 = 60;
        Params.K_CH4_SO4 = 300;
        Params.k_AOM = 5;
    
        Params.Fe_inventory_factor = 1.0;
        Params.KFEMonod = 150;
        Params.kFeS = 10;
        Params.kFeOx = 5;
    
    case 'georgia_stl2_trial2'
        [Config, Params, Case] = Validation_Case('georgia_stl2_trial');
        Case.name = name;
        Case.note = 'Georgia STL2 adjusted from trial: lower sulfate boundary, broader sulfate reduction, delayed Fe2+ peak.';
    
        % Boundary conditions
        Config.SO4init = 22000;
    
        % Keep trial-like pH structure
        Config.DICinit = 2500;
        Config.HCO3init = 2500;
        Config.Calcium = 9000;
    
        % Broaden OM degradation without creating huge OM stock
        Config.vbottom = 0.25;
        Config.F_OM_total = 0.006;
        Config.k_sed_scale = 45;
        Config.f_lab = 0.55;
    
        % Slightly stronger shallow exchange to delay Fe2+ buildup
        Config.Bioirrig_top = 2.0;
        Config.Bioirrig_scale = 2.5;
    
        % Sulfate transport: lower than trial, not as low as v2
        Params.DSO4 = 80;
        Params.k_SO4 = 100;
    
        % Fe source: stronger than trial, weaker than v3
        Config.F_FeOx = 6.0;
        Params.Fe_inventory_factor = 0.6;
        Params.KFEMonod = 150;
    
        % Fe sink: moderate, not too strong
        Params.kFeS = 8;
        Params.kFeOx = 5;
    
    case 'georgia_stl2_trial3'
        [Config, Params, Case] = Validation_Case('georgia_stl2_trial');
        Case.name = name;
        Case.note = 'Georgia STL2: broader sulfate reduction with delayed Fe2+ peak.';
    
        Config.SO4init = 22000;
    
        % keep trial-like carbonate condition
        Config.DICinit = 2500;
        Config.HCO3init = 2500;
        Config.Calcium = 9000;
    
        % spread OM degradation deeper
        Config.F_OM_total = 0.008;
        Config.k_sed_scale = 20;
        Config.f_lab = 0.30;
        Config.vbottom = 0.25;
    
        % mild shallow exchange to delay Fe2+ peak, not too strong
        Config.Bioirrig_top = 1.5;
        Config.Bioirrig_scale = 3.5;
    
        % sulfate transport
        Params.DSO4 = 60;
        Params.k_SO4 = 80;
    
        % Fe source: stronger than trial, weaker/less top-heavy than v3
        Config.F_FeOx = 6.0;
        Params.Fe_inventory_factor = 0.6;
        Params.KFEMonod = 250;
        Params.kFeS = 8;
        Params.kFeOx = 5;
    
    case 'georgia_stl2_trial4'
        [Config, Params, Case] = Validation_Case('georgia_stl2_trial3');
        Case.name = name;
        Case.note = 'Georgia STL2: preserve delayed Fe2+ peak, strengthen mid-depth Fe sink and pH recovery.';
    
        Params.KFEMonod = 700;
        Params.kFeS = 20;
    
        Params.k_SO4 = 40;
    
        Config.DICinit = 2700;
        Config.HCO3init = 2600;

    case 'georgia_sap1_trial'

        Case.name = name;
        Case.note = 'Georgia SAP1 trial: sulfidic sulfate-reducing marsh creek-bank sediment.';
    
        % Domain / environment
        Config.Lbottom = 50;
        Config.T_future = 25;
        Config.Salinity = 20;
    
        % Top boundary / initial porewater
        Config.O2init = 180;
        Config.SO4init = 20000;
        Config.CH4init = 0;
    
        Config.DICinit = 2000;
        Config.HCO3init = 3200;   % ALK top, legacy field name
        Config.Calcium = 8000;
    
        Config.Feinit = 0;
        Config.HSinit = 0;
    
        % Physical structure
        Config.vbottom = 0.45;
    
        Config.porostop = 0.90;
        Config.porosbottom = 0.76;
        Config.porosscale = 12;
    
        Config.Bioturbtop = 0.3;
        Config.Bioturbbottom = 0.02;
        Config.bioturbscale = 4;
    
        Config.Bioirrig_top = 0.15;
        Config.Bioirrig_bottom = 0;
        Config.Bioirrig_scale = 6;
    
        % OM supply and reactivity
        Config.F_OM_total = 0.035;
        Config.k_sed_scale = 12;
        Config.f_lab = 0.35;
    
        % Keep Fe weak: SAP1 observation Fe2+ ~0
        Config.F_FeOx = 0.05;
        Params.Fe_inventory_factor = 0.02;
        Params.KFEMonod = 1000;
    
        % Carbonate / bubbling
        Config.F_CaCO3 = 0;
        Config.water_depth_m = 0.5;
        Config.use_CH4_bubbling = true;
    
        % Transport
        Params.DO2 = 120;
        Params.DSO4 = 100;
        Params.DH2S = 250;
        Params.DCH4 = 250;
        Params.DHCO3 = 220;
        Params.DCa = Params.DHCO3;
    
        % Reactions
        Params.k_SO4 = 100;
        Params.K_CH4_SO4 = 1000;
        Params.k_AOM = 0.5;
    
        Params.kFeS = 30;
        Params.kFeOx = 5;
    
    case 'georgia_sap1_final'
        [Config, Params, Case] = Validation_Case('georgia_sap1_trial');
        Case.name = name;
        Case.note = 'Georgia SAP1 trial2: stronger distributed sulfate reduction with lower alkalinity.';
    
        % Boundary carbonate: lower pH / avoid supersaturation
        Config.DICinit = 2200;
        Config.HCO3init = 2300;
        Config.Calcium = 6500;

        % Config.n=501;
    
        % Stronger but more distributed OM throughput
        Config.F_OM_total = 0.060;
        Config.k_sed_scale = 8;
        Config.f_lab = 0.25;
        Config.vbottom = 0.70;
    
        % Keep weak exchange
        Config.Bioirrig_top = 0.05;
        Config.Bioirrig_scale = 8;
    
        % Make sulfate depletion easier and retain products
        Params.DSO4 = 50;
        Params.DH2S = 120;
        Params.DHCO3 = 150;
        Params.DCH4 = 180;
    
        % Allow SRR and methane transition
        Params.k_SO4 = 60;
        Params.K_CH4_SO4 = 500;
        Params.k_AOM = 0.2;
    
        % Keep Fe suppressed
        Config.F_FeOx = 0.02;
        Params.Fe_inventory_factor = 0.01;
        Params.KFEMonod = 1500;
        Params.kFeS = 50;
    
    otherwise
        error('Unknown validation case: %s', name);
end
end