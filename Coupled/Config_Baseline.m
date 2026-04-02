function Config = Config_Baseline()
% CONFIG_BASELINE
% Site / scenario specific settings for the OLD sequential core.

    % ---------------- Domain ----------------
    Config.Lbottom = 30;          % cm
    Config.n = 101;
    Config.nmesh = 1000;

    % ---------------- Physical structure ----------------
    Config.vbottom = 0.5;         % cm / yr   (river-audited baseline)
    Config.vbottom_fluid = 0;     % cm / yr

    Config.porostop = 0.9;
    Config.porosbottom = 0.7;
    Config.porosscale = 3;

    Config.Bioturbtop = 10;       % cm2 / yr
    Config.Bioturbbottom = 1;     % cm2 / yr
    Config.bioturbscale = 3;

    Config.Bioirrig_top = 100;    % 1 / yr
    Config.Bioirrig_bottom = 0;
    Config.Bioirrig_scale = 0.75;

    % ---------------- Boundary concentrations ----------------
    Config.O2init   = 150;        % uM
    Config.SO4init  = 200;        % uM
    Config.DICinit  = 1000;       % uM
    Config.HCO3init = 950;        % uM
    Config.Calcium  = 1000;       % uM
    Config.CH4init  = 0;          % uM
    Config.Feinit   = 0;          % uM
    Config.HSinit   = 0;          % uM
    Config.Pinitial = 0;          % uM

    % ---------------- Fluxes ----------------
    Config.NPP = 200;             % g / m2 / yr
    Config.BE  = 0.1;
    Config.F_FeOx  = 2;          % mmol / m2 / d
    Config.F_CaCO3 = 10;          % g / m2 / yr

    % ---------------- Temperature / OM age ----------------
    Config.T_future = 25;
    Config.ageinit = 0.1;
    Config.age_root = 1;

    % ---------------- Root-zone extras ----------------
    Config.DOC_root_1 = 0;
    Config.O2_root_1  = 0;
    Config.POC_root_1 = 0;

    % ---------------- Step 2 hydro inputs ----------------
    % keep these even if old core only weakly uses them for now
    Config.X_f     = 0.3;
    Config.phi_c   = 0.4;
    Config.phi_f   = 0.6;
    Config.H       = 2.0;         % m
    Config.S       = 1e-4;        % m / m
    Config.alpha_L = 1.0;         % cm
    Config.w_s     = 0.34;        % mm / s
    Config.tau_c   = 0.3;         % N / m2

    % ---------------- switches ----------------
    Config.use_constant_porosity = false;
    Config.constant_porosity = 0.75;

    Config.use_hydro_phi = false;
    Config.use_hydro_npp_multiplier = false;
    Config.use_hydro_diffusion_multiplier = false;

    % safety cap for weak-coupling injection
    Config.max_diffusion_multiplier = 2.0;
end