function Params = Params()
    % PARAMS_STATIC: Universal thermodynamic, kinetic, and physical constants.
    % These values remain locked during Monte Carlo/USGS upscaling.
    
    % -- Physical Constants --
    Params.rho = 2.73;               % Solid density [g/cm3]
    Params.Mineral_Mass = 215;       % Background mineral mass factor
    
    % -- Diffusion Coefficients (cm2/yr) --
    Params.DO2   = 300;
    Params.DSO4  = 300;
    Params.DH2S  = 300;
    Params.DCH4  = 300;
    Params.DHCO3 = 400;
    Params.DPO4  = 400;
    
    % -- Half-Saturation Constants (uM) --
    Params.k_O2        = 2;          % Oxygen Monod
    Params.k_SO4       = 20;         % Sulfate Monod
    Params.KFEMonod    = 200;        % FeOOH Monod (umol/g)
    Params.K_CH4_SO4   = 100;        % AOM SO4 half-saturation
    Params.K_CH4_O2    = 1;          % Aerobic methanotrophy O2 half-saturation
    
    % -- Kinetic Rate Constants --
    Params.Kreox         = 500;      % Sulfide oxidation rate [1/umol/L/yr]
    Params.kFeOx         = 10;       % Fe(II) oxidation [1/umol/L/yr]
    Params.kFeS          = 10;       % FeS precipitation [1/umol/L/yr]
    Params.k_AOM         = 0.1;      % Anaerobic Ox. of Methane [1/yr]
    Params.k_aerobic_CH4 = 1;        % Aerobic Ox. of Methane [1/yr]
    
    % -- Carbonate System Thermodynamics --
    Params.Ksp_ca         = 3000;    % Calcite solubility product
    Params.k_calcite      = 1;       % Precipitation rate constant
    Params.k_calcite_dis1 = 0.005;   % Dissolution rate 1
    Params.k_calcite_dis2 = 10;      % Dissolution rate 2
    Params.n_power_CaCO31 = 1.76;    % Precipitation order
    Params.n_power_CaCO32 = 0.11;    % Dissolution order 1
    Params.n_power_CaCO33 = 4;       % Dissolution order 2
    Params.Calcium_activity = 0.6;
    Params.CO3_activity     = 0.6;
    
    % -- Stoichiometry & Organic Matter --
    Params.P_C_ratio = 0.0094;
    Params.Q10       = 2;            % Temperature scaling
    Params.T_ref     = 25;           % Reference temperature for Q10
end