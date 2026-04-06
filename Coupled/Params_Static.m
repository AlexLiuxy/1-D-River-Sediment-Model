function Params = Params_Static()
% PARAMS_STATIC
% Constants that should not change from site to site unless explicitly audited.

    Params.rho = 2.73;           % g / cm3
    Params.k_O2 = 2;             % uM
    Params.k_SO4 = 20;           % uM
    Params.KFEMonod = 200;       % umol / g

    Params.DSO4 = 310;%300;           % cm2 / yr
    Params.DCH4 = 300;           % cm2 / yr
    Params.DH2S = 300;           % cm2 / yr
    Params.DO2  = 300;           % cm2 / yr
    Params.DHCO3 = 400;          % cm2 / yr
    Params.DPO4  = 400;          % cm2 / yr

    Params.Kreox = 500;          % 1 / umol / L / yr
    Params.kFeOx = 10;           % 1 / umol / L / yr
    Params.kFeS  = 10;           % 1 / umol / L / yr

    Params.K_CH4_SO4   = 100;    % uM
    Params.K_CH4_O2    = 1;      % uM
    Params.k_AOM       = 1.0;    % 1 / yr
    Params.k_aerobic_CH4 = 6;    % 1 / yr

    Params.Ksp_ca = 3000;        % uM^2
    Params.k_calcite = 1;
    Params.k_calcite_dis1 = 0.005;
    Params.k_calcite_dis2 = 10;
    Params.n_power_CaCO31 = 1.76;
    Params.n_power_CaCO32 = 0.11;
    Params.n_power_CaCO33 = 4;

    Params.Calcium_activity = 0.6;
    Params.CO3_activity     = 0.6;

    Params.P_C_ratio = 0.0094;
    Params.Q10   = 2;
    Params.T_ref = 25;

    % extras already used by old core
    Params.KFeS = 2500;
    Params.K_HS = 7;
    Params.kapatite = 0.05;
    Params.Kviv = 3e6;
    Params.alpha_viv = 1.5;
    Params.kviv = 1.7e-22;
end