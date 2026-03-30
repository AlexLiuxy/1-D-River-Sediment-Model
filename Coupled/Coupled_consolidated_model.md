# MATLAB Model Source Code

## File: BC_River_System.m
```matlab
function res = BC_River_System(Ya, Yb, Config, Hydro)
    % Top boundary fluxes for solid phases derived from NPP and Burial Efficiency
%     phi_top = Config.porostop;
    phi_top = Hydro.phi;
    rho = 2.73; % Params.rho
    v_s_top = Config.vbottom * (1 - Config.porosbottom) / (1 - phi_top);
%     CorgInit = (Config.BE * Config.NPP * 1E-4) / (v_s_top * rho * (1 - phi_top));
%     FeOxInit = 36.5 * Config.F_FeOx * (phi_top / (1 - phi_top)) / v_s_top / rho;
%     CaCO3Init = 1E-4 * Config.F_CaCO3 * (phi_top / (1 - phi_top)) / v_s_top / rho;
    % Apply scouring switch here as well for absolute consistency
    Effective_NPP = Config.NPP * Hydro.POM_Flux_Multiplier;
    CorgInit = (Config.BE * Effective_NPP * 1E-4) / (v_s_top * rho * (1 - phi_top));
    FeOxInit = 36.5 * Config.F_FeOx * (phi_top / (1 - phi_top)) / v_s_top / rho;
    CaCO3Init = 1E-4 * Config.F_CaCO3 * (phi_top / (1 - phi_top)) / v_s_top / rho;
    res = [
        % 1. Top Boundaries (Dirichlet)
        Ya(1)  - Config.O2init;
        Ya(3)  - Config.SO4init;
        Ya(5)  - Config.CH4init;
        Ya(7)  - Config.DICinit;
        Ya(9)  - Config.HCO3init;
        Ya(11) - Config.Calcium;
        Ya(13) - Config.Feinit;
        Ya(15) - Config.HSinit;
        % Solid top boundaries
        Ya(17) - CorgInit;
        Ya(19) - FeOxInit;
        Ya(21) - CaCO3Init;
        % 2. Bottom Boundaries (Neumann, dC/dx = 0)
        Yb(2);
        Yb(4);
        Yb(6);
        Yb(8);
        Yb(10);
        Yb(12);
        Yb(14);
        Yb(16);
        Yb(18);
        Yb(20);
        Yb(22);
    ];
end
```

## File: Config.m
```matlab
function Config = Config()
    % CONFIG_USGS_BASELINE: Site-specific parameters mapped from USGS data.
    % These will be randomized via Latin Hypercube Sampling in Step 4.
    % -- Domain & Hydrology --
    Config.Lbottom       = 30;       % Max depth [cm]
    Config.n_grid        = 101;      % Number of spatial nodes
    Config.vbottom       = 0.5;      % Sedimentation rate [cm/yr] (AUDITED: 100 -> 0.5)
    Config.vbottom_fluid = 0;        % Groundwater upwelling [cm/yr]
    % -- Sediment Physical Properties --
    Config.porostop      = 0.9;      % Porosity at SWI
    Config.porosbottom   = 0.7;      % Porosity at depth
    Config.porosscale    = 3;        % e-folding depth for porosity [cm]
    Config.Bioturbtop    = 10;       % SWI Bioturbation [cm2/yr]
    Config.Bioturbbottom = 1;        % Deep Bioturbation [cm2/yr]
    Config.bioturbscale  = 3;        % e-folding depth for bioturbation [cm]
    Config.Bioirrig_top  = 100;      % SWI Bioirrigation [1/yr]
    Config.Bioirrig_scale= 0.75;     % e-folding depth for bioirrigation [cm]
    % -- Top Boundary Conditions (SWI Concentrations) --
    Config.O2init   = 150;           % Oxygen [uM]
    Config.SO4init  = 200;           % Sulfate [uM] (AUDITED: 3000 -> 200)
    Config.DICinit  = 1000;          % DIC [uM]
    Config.HCO3init = 950;           % Alkalinity [uM]
    Config.Calcium  = 1000;          % Calcium [uM]
    Config.CH4init  = 0;             % Methane [uM]
    Config.Feinit   = 0;             % Dissolved Fe [uM]
    Config.HSinit   = 0;             % Sulfide [uM]
    % -- Depositional Fluxes --
    Config.NPP      = 600;           % Net Primary Production [g/m2/yr]
    Config.BE       = 0.2;           % Burial Efficiency
    Config.F_FeOx   = 10;            % FeOx flux [mmol/m2/d]
    Config.F_CaCO3  = 10;            % CaCO3 flux [g/m2/yr]
    % -- Climate / Environmental --
    Config.T_future = 25;            % Local water temperature [C]
    Config.ageinit  = 0.1;           % Initial age of OM at SWI [yr]
    % -- Hydrodynamics & PPT Parameters (Step 2 Integration) --
    Config.X_f     = 0.3;      % Volumetric fraction of fine grains
    Config.phi_c   = 0.4;      % Porosity of coarse framework
    Config.phi_f   = 0.6;      % Porosity of fine matrix
    Config.H       = 2.0;      % River depth [m]
    Config.S       = 0.0001;   % River slope [m/m]
    Config.alpha_L = 1.0;      % Longitudinal dispersivity [cm]
    Config.w_s     = 0.34;     % Characteristic settling velocity [mm/s]
    Config.tau_c   = 0.3;      % Critical shear stress [N/m2]
end
```

## File: Hydro_Preprocessor.m
```matlab
function Hydro = Hydro_Preprocessor(Config)
% HYDRO_PREPROCESSOR
% Re-engineered based on Krone (1962) Deposition Law and Boano (2014) Hyporheic transport
    %% 1. Unified Packing Model for Porosity (El-Husseiny, 2021)
    X_f   = Config.X_f;
    phi_c = Config.phi_c;
    phi_f = Config.phi_f;
    if X_f < phi_c
        phi_mix = phi_c - X_f * (1 - phi_f);
    else
        phi_mix = X_f * phi_f;
    end
    Hydro.phi = phi_mix;
    %% 2. Hydrodynamics & Bed Shear Stress
    g = 9.81;
    H = Config.H;
    S = Config.S;
    rho_w = 1000; % Water density [kg/m3]
    u_star_m = sqrt(g * H * S); % Friction velocity [m/s]
    tau_b = rho_w * (u_star_m^2); % Bed shear stress [N/m2]
    %% 3. Literature-Corrected POM Deposition (Krone 1962)
    % Determines net effective deposition based on critical shear stress.
    tau_cd = Config.tau_c; % Critical shear stress for deposition [N/m2]
    if tau_b >= tau_cd
        Hydro.POM_Flux_Multiplier = 0; % Fully erosive, no net deposition
    else
        Hydro.POM_Flux_Multiplier = 1 - (tau_b / tau_cd); % Partial deposition
    end
    %% 4. Literature-Corrected Hyporheic Diffusion (Boano 2014 / Chen 2021)
    % Mechanical dispersion uses pore-water velocity (v_pore), attenuated from u_star.
    attenuation_factor = 1e-3;
    v_pore_m = u_star_m * attenuation_factor; % [m/s]
    v_pore_cm = v_pore_m * 100; % [cm/s]
    alpha_L = Config.alpha_L; % Longitudinal dispersivity [cm]
    Seconds_per_year = 3600 * 24 * 365;
    % Surface hyporheic dispersion coefficient [cm2/yr]
    Hydro.D_hyp_surface = (alpha_L * v_pore_cm) * Seconds_per_year;
    % Penetration depth for turbulence decay (e.g., decays over top 5 cm)
    Hydro.D_hyp_decay_depth = 5.0;
    %% 5. Export hydrodynamics for logging
    Hydro.u_star = u_star_m;
    Hydro.tau_b  = tau_b;
    Hydro.v_pore = v_pore_m;
end
```

## File: Main_River_Model.m
```matlab
clear; clc; close all;
tic
% --- 1. Load Configurations and Parameters ---
% Replacing legacy 'global' variables with structured data arrays
Config = Config();
Params = Params();
% --- 2. Define Spatial Grid ---
x_mesh = linspace(0, Config.Lbottom, Config.n_grid);
% --- 3. Establish the Index Map and Execute Hydro Preprocessor ---
% bvp4c requires converting 2nd-order ODEs into a system of 1st-order ODEs.
% Y(odd)  = Concentration of species (uM)
% Y(even) = Concentration gradient dC/dx
%
% Species Order:
% 1: O2, 2: SO4, 3: CH4, 4: DIC, 5: HCO3(Alk), 6: Ca, 7: Fe, 8: HS
num_species = 11;
total_vars  = num_species * 2;
% Directly pass the Config struct to extract hydro parameters
Hydro = Hydro_Preprocessor(Config);
% % --- 4. Initialize bvp4c Solver ---
% % Flat initial guess for all 16 variables (8 species + 8 gradients)
% initial_guess = zeros(1, total_vars);
%
% % Set surface boundary conditions as initial guess to help convergence
% initial_guess(1)  = Config.O2init;
% initial_guess(3)  = Config.SO4init;
% initial_guess(5)  = Config.CH4init;
% initial_guess(7)  = Config.DICinit;
% initial_guess(9)  = Config.HCO3init;
% initial_guess(11) = Config.Calcium;
% initial_guess(13) = Config.Feinit;
% initial_guess(15) = Config.HSinit;
%
% % Solid phase dynamic calculation
% phi_top = Config.porostop;
% v_s_top = Config.vbottom * (1 - Config.porosbottom) / (1 - phi_top);
%
% initial_guess(17) = (Config.BE * Config.NPP * 1E-4) / (v_s_top * Params.rho * (1 - phi_top)); % C_org
% initial_guess(19) = 36.5 * Config.F_FeOx * (phi_top / (1 - phi_top)) / v_s_top / Params.rho;  % FeOx
% initial_guess(21) = 1E-4 * Config.F_CaCO3 * (phi_top / (1 - phi_top)) / v_s_top / Params.rho; % CaCO3
%
% solinit = bvpinit(x_mesh, initial_guess);
% --- 4. Initialize bvp4c Solver ---
% % Calculate top boundary fluxes for solid phase initial guesses
% phi_top = Config.porostop;
% v_s_top = Config.vbottom * (1 - Config.porosbottom) / (1 - phi_top);
%
% Corg_top  = (Config.BE * Config.NPP * 1E-4) / (v_s_top * Params.rho * (1 - phi_top));
% FeOx_top  = 36.5 * Config.F_FeOx * (phi_top / (1 - phi_top)) / v_s_top / Params.rho;
% CaCO3_top = 1E-4 * Config.F_CaCO3 * (phi_top / (1 - phi_top)) / v_s_top / Params.rho;
phi_top = Hydro.phi; % Use unified porosity from PPT
v_s_top = Config.vbottom * (1 - Config.porosbottom) / (1 - phi_top);
% Apply the Scouring Switch (POM_Flux_Multiplier) from Lamb (2020) logic
Effective_NPP = Config.NPP * Hydro.POM_Flux_Multiplier;
Corg_top  = (Config.BE * Effective_NPP * 1E-4) / (v_s_top * Params.rho * (1 - phi_top));
FeOx_top  = 36.5 * Config.F_FeOx * (phi_top / (1 - phi_top)) / v_s_top / Params.rho;
CaCO3_top = 1E-4 * Config.F_CaCO3 * (phi_top / (1 - phi_top)) / v_s_top / Params.rho;
% Create a depth-dependent initial guess matrix to prevent stiff boundary layer blowups
% Rows = variables (22), Columns = spatial nodes (n_grid)
initial_guess = zeros(total_vars, length(x_mesh));
for i = 1:length(x_mesh)
    depth = x_mesh(i);
    % Dissolved phases (Flat initial guess is mathematically stable due to high diffusion)
    initial_guess(1, i)  = Config.O2init;
    initial_guess(3, i)  = Config.SO4init;
    initial_guess(5, i)  = Config.CH4init;
    initial_guess(7, i)  = Config.DICinit;
    initial_guess(9, i)  = Config.HCO3init;
    initial_guess(11, i) = Config.Calcium;
    initial_guess(13, i) = Config.Feinit;
    initial_guess(15, i) = Config.HSinit;
    % Solid phases (Exponential decay guess required to match bioturbation physics)
    % This strictly suppresses the 'Unable to meet tolerance' warning and massive residual
    decay_factor = exp(-depth / Config.bioturbscale);
    initial_guess(17, i) = Corg_top * decay_factor;
    initial_guess(19, i) = FeOx_top * decay_factor;
    initial_guess(21, i) = CaCO3_top * decay_factor; % Fix: Apply decay to CaCO3 as well; CaCO3_top; % CaCO3 dissolution is slower, linear/flat guess is acceptable
end
solinit.x = x_mesh;
solinit.y = initial_guess;
% Tight tolerances required for stiff biogeochemical reaction fronts
options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'NMax', 5000);
% --- 5. Execute Fully Coupled Solver ---
fprintf('Running coupled BVP solver...\n');
sol = bvp4c(@(x, Y) ODE_River_System(x, Y, Config, Params, Hydro), ...
            @(Ya, Yb) BC_River_System(Ya, Yb, Config, Hydro), ...
            solinit, options);
fprintf('Simulation converged successfully.\n');
% --- 6. Data Extraction (Unpacking all 8 species) ---
depth_cm = sol.x;
O2_prof  = sol.y(1, :);
SO4_prof = sol.y(3, :);
CH4_prof = sol.y(5, :);
DIC_prof = sol.y(7, :);
ALK_prof = sol.y(9, :);
Ca_prof  = sol.y(11, :);
Fe_prof  = sol.y(13, :);
HS_prof  = sol.y(15, :);
% --- 7. Post-Processing: Recalculate pH and Omega for plotting ---
% Initialize arrays for derived variables
pH_prof = zeros(size(depth_cm));
Omega_prof = zeros(size(depth_cm));
% Constants for carbonate system (must match ODE_River_System)
K1 = 1.18e-6; K2 = 4.36e-11; Kw = 1e-14; Kb = 2.3e-9; B_T = 4e-4;
for i = 1:length(depth_cm)
    current_DIC = DIC_prof(i);
    current_ALK = ALK_prof(i);
    current_Ca  = Ca_prof(i);
    % Polynomial roots for H+
    p5 = -1;
    p4 = -current_ALK - Kb - K1;
    p3 = current_DIC*K1 - current_ALK*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = current_DIC*(Kb*K1 + 2*K1*K2) - current_ALK*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*current_DIC*Kb*K1*K2 - current_ALK*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    r = roots([p5 p4 p3 p2 p1 p0]);
    H = max(real(r));
    pH_prof(i) = -log10(H);
    % Carbonate speciation and Omega
    CO3 = current_DIC / (1 + H/K2 + H^2/(K1*K2));
    Omega_prof(i) = (current_Ca * CO3) / Params.Ksp_ca;
end
toc
% =========================================================================
% --- Post-Processing & EXACT Original Visualization ---
% =========================================================================
% 1. Map new arrays to your exact original variable names
z_sed     = sol.x;
Oxygen    = sol.y(1, :);
Sulfate   = sol.y(3, :);
CH4       = sol.y(5, :);
C_DIC     = sol.y(7, :);
C_alka    = sol.y(9, :);
C_Fe      = sol.y(13, :);
C_HS      = sol.y(15, :);
C_organic = sol.y(17, :);
FeooH     = sol.y(19, :);
CaCO3     = sol.y(21, :);
% 2. Recalculate intermediate rates for plotting
pH = zeros(size(z_sed));
sigma_carb = zeros(size(z_sed));
C_H2CO3 = zeros(size(z_sed));
RC = zeros(size(z_sed));
R_SRR = zeros(size(z_sed));
K1 = 1.18e-6; K2 = 4.36e-11; Kw = 1e-14; Kb = 2.3e-9; B_T = 4e-4;
for i = 1:length(z_sed)
    x_curr = z_sed(i);
    phi = Config.porosbottom + (Config.porostop - Config.porosbottom) * exp(-x_curr / Config.porosscale);
    v_burial_s = Config.vbottom * (1 - Config.porosbottom) / (1 - phi);
    age = Config.ageinit + x_curr / v_burial_s;
    k_sed = 10^(-0.95 * log10(age) - 0.81);
    Temp_factor = Params.Q10^((Config.T_future - Params.T_ref)/10);
    RC_val = Temp_factor * k_sed * C_organic(i) * Params.rho * ((1 - phi) / 12);
    RC(i) = RC_val; % keeping in base units for now
    Inhib_O2 = (Params.k_O2 / (Oxygen(i) + Params.k_O2));
    R_SRR(i) = (RC_val * 1e9) * Inhib_O2 * (Sulfate(i) / (Sulfate(i) + Params.k_SO4));
    p5 = -1; p4 = -C_alka(i) - Kb - K1;
    p3 = C_DIC(i)*K1 - C_alka(i)*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = C_DIC(i)*(Kb*K1 + 2*K1*K2) - C_alka(i)*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*C_DIC(i)*Kb*K1*K2 - C_alka(i)*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    r = roots([p5 p4 p3 p2 p1 p0]);
    H = max(real(r));
    pH(i) = -log10(H);
    C_H2CO3(i) = C_DIC(i) / (1 + K1/H + K1*K2/H^2);
    CO3 = C_DIC(i) / (1 + H/K2 + H^2/(K1*K2));
    sigma_carb(i) = (Params.Calcium_activity * Config.Calcium * Params.CO3_activity * CO3) / Params.Ksp_ca - 1;
end
% BEsed_org = C_organic ./ C_organic(1);
BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
POC_root = zeros(size(z_sed)); % Was 0 in your code
% --- EXACT ORIGINAL PLOTTING LOGIC ---
clf;
n_plot = 6; % number of plots in each row
m_plot = 3; % number of total rows
% Organic
subplot(m_plot,n_plot,1);
plot((C_organic + POC_root).*100,z_sed,'lineWidth',2); axis ij
title('Organic (%gDw)')
ylabel('Depth (cm)');
box on
% Oxygen
subplot(m_plot,n_plot,2);
plot(Oxygen,z_sed,'lineWidth',2); axis ij
title('[O_2] (\muM)')
box on
grid on
 ax.LineWidth = 2;
% Iron
subplot(m_plot,n_plot,3);
plot(C_Fe,z_sed,'lineWidth',2); axis ij
title('[Fe^{2+}] (\muM)')
box on
grid on
 ax.LineWidth = 2;
% Sulfate
subplot(m_plot,n_plot,4);
plot(Sulfate,z_sed,'lineWidth',2); axis ij
title('[SO_4] (\muM)')
box on
grid on
 ax.LineWidth = 2;
% Sulfide
subplot(m_plot,n_plot,5);
plot(C_HS,z_sed,'lineWidth',2); axis ij
title('[H_2S] (\muM)')
box on
grid on
 ax.LineWidth = 2;
% Methane
subplot(m_plot,n_plot,6);
plot(CH4,z_sed,'lineWidth',2); axis ij
title('[CH_4] (\muM)')
box on
grid on
 ax.LineWidth = 2;
% CaCO3
subplot(m_plot,n_plot,7);
plot(CaCO3.*100,z_sed,'lineWidth',2); axis ij
title('CaCO3')
ylabel('Depth (cm)');
box on
grid on
 ax.LineWidth = 2;
% DIC
subplot(m_plot,n_plot,8);
plot(C_DIC,z_sed,'lineWidth',2); axis ij
title('DIC (\muM)')
box on
grid on
 ax.LineWidth = 2;
% ALK
subplot(m_plot,n_plot,9);
plot(C_alka,z_sed,'lineWidth',2); axis ij
title('ALK (\muM)')
box on
grid on
 ax.LineWidth = 2;
% Carbonic Acid
subplot(m_plot,n_plot,10);
plot(C_H2CO3,z_sed,'lineWidth',2); axis ij
title('Carb Acid (\muM)')
box on
grid on
ax.LineWidth = 2;
% pH
subplot(m_plot,n_plot,11);
plot(pH,z_sed,'lineWidth',2); axis ij
title('pH')
box on
grid on
 ax.LineWidth = 2;
% Organic degradation rate
subplot(m_plot,n_plot,12);
plot(RC.* 1E9,z_sed,'lineWidth',2); axis ij  %umol/l/year
title('Mineralization Rate (\mumol/l/year)')
box on
grid on
 ax.LineWidth = 2;
% Burial efficiency of Organic
subplot(m_plot,n_plot,13);
plot(BEsed_org.*100,z_sed,'lineWidth',2); axis ij  %umol/l/year
title('Burial Efficiency of organic matter')
ylabel('Depth (cm)');
box on
grid on
 ax.LineWidth = 2;
% Sulfate Reduction Rate
subplot(m_plot,n_plot,14);
plot((0.5.*R_SRR)./365,z_sed,'lineWidth',2); axis ij %umol/l/year
title('Sulfate Reduction Rate (nmol/cm3/d)')
box on
grid on
 ax.LineWidth = 2;
% Caclite saturation
subplot(m_plot,n_plot,16);
plot(sigma_carb,z_sed,'lineWidth',2); axis ij
title('Calciite saturation (\Omega - 1)')
box on
% Fe(III)
subplot(m_plot,n_plot,17);
plot(FeooH,z_sed,'lineWidth',2); axis ij
title('Fe(III) (\mumol/gr)')
box on
```

## File: ODE_River_System.m
```matlab
function dYdx = ODE_River_System(x, Y, Config, Params, Hydro)
    % 1. Physical & Hydrological Parameters
%     phi = Config.porosbottom + (Config.porostop - Config.porosbottom) * exp(-x / Config.porosscale);
    phi = Hydro.phi;
    alpha_bio = Config.Bioirrig_top * exp(-x / Config.Bioirrig_scale);
    v_burial_s = Config.vbottom * (1 - Config.porosbottom) / (1 - phi); % Solid burial velocity
    v_burial_f = Config.vbottom_fluid * (1 + Config.porosbottom) / (1 + phi); % Fluid burial velocity
    Db = Config.Bioturbbottom + (Config.Bioturbtop - Config.Bioturbbottom) * exp(-x / Config.bioturbscale);
% Enhance liquid diffusion with Hyporheic mechanical dispersion (decaying with depth)
    % Assume hyporheic influence penetrates ~5 cm deep
    local_D_hyp = Hydro.D_hyp_surface * exp(-x ./ Hydro.D_hyp_decay_depth);
    D_O2_eff  = Params.DO2 + local_D_hyp;
    D_SO4_eff = Params.DSO4 + local_D_hyp;
    D_CH4_eff = Params.DCH4 + local_D_hyp;
    D_DIC_eff = Params.DHCO3 + local_D_hyp;
    D_HS_eff  = Params.DH2S + local_D_hyp;
    % Solid phase volume fraction ratio
    solid_ratio = (1 - phi) / phi;
    solid_rho_factor = solid_ratio * Params.rho * 1e3; % Conversion factor for solid to dissolved
    % Age of sediment for reactivity (approximation without iterative loop)
    age = Config.ageinit + x / v_burial_s;
    k_sed = 10^(-0.95 * log10(age) - 0.81);
    Temp_factor = Params.Q10^((Config.T_future - Params.T_ref)/10);
    % 2. Unpack Variables (11 Species * 2 = 22 variables)
    O2    = max(Y(1), 1e-9);  dO2dx    = Y(2);
    SO4   = max(Y(3), 1e-9);  dSO4dx   = Y(4);
    CH4   = max(Y(5), 1e-9);  dCH4dx   = Y(6);
    DIC   = max(Y(7), 1e-9);  dDICdx   = Y(8);
    ALK   = max(Y(9), 1e-9);  dALKdx   = Y(10);
    Ca    = max(Y(11), 1e-9); dCadx    = Y(12);
    Fe    = max(Y(13), 1e-9); dFedx    = Y(14);
    HS    = max(Y(15), 1e-9); dHSdx    = Y(16);
    C_org = max(Y(17), 1e-12);dCorgdx  = Y(18); % Solid Organic Carbon (g/gDw)
    FeOx  = max(Y(19), 1e-9); dFeOxdx  = Y(20); % Solid Fe(OH)3 (umol/g)
    CaCO3 = max(Y(21), 0);    dCaCO3dx = Y(22); % Solid CaCO3
    % 3. Carbonate System (Recalculate real-time pH)
    K1 = 1.18e-6; K2 = 4.36e-11; Kw = 1e-14; Kb = 2.3e-9; B_T = 4e-4;
    p5 = -1; p4 = -ALK - Kb - K1;
    p3 = DIC*K1 - ALK*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = DIC*(Kb*K1 + 2*K1*K2) - ALK*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*DIC*Kb*K1*K2 - ALK*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    r = roots([p5 p4 p3 p2 p1 p0]);
    H = max(real(r));
    CO3 = DIC / (1 + H/K2 + H^2/(K1*K2));
    Omega_ca = (Params.Calcium_activity * Ca * Params.CO3_activity * CO3) / Params.Ksp_ca;
    % CaCO3 precipitation/dissolution
    unit_conv = 1 / (1E3 * 1E6 * 1E-2 * Params.rho * (1 - phi));
    if Omega_ca > 1
        R_carb_form = (Omega_ca - 1)^Params.n_power_CaCO31 * Params.k_calcite * unit_conv;
        R_carb_disso = 0;
    else
        R_carb_form = 0;
        if Omega_ca > 0.8
            R_carb_disso = (1 - Omega_ca)^Params.n_power_CaCO32 * Params.k_calcite_dis1 * CaCO3;
        else
            R_carb_disso = (1 - Omega_ca)^Params.n_power_CaCO33 * Params.k_calcite_dis2 * CaCO3;
        end
    end
    R1_carb_total = R_carb_form - R_carb_disso * (1E3 * 1E6 * 1E-2 * Params.rho * (1 - phi));
    % 4. Primary Organic Matter Mineralization (Redox Cascade)
    % Total Mineralization Rate (molCorg/cm3/yr -> scaled to umol/L/yr below)
    RC = Temp_factor * k_sed * C_org * Params.rho * ((1 - phi) / 12);
    % Inhibition terms based on standard Monod kinetics
    Inhib_O2 = (Params.k_O2 / (O2 + Params.k_O2));
    Inhib_Fe = (Params.KFEMonod / (FeOx + Params.KFEMonod));
    % Partitioning the total rate (RC) into specific electron acceptors [umol/L/yr]
    R_respi = (RC * 1e9) * (O2 / (O2 + Params.k_O2));
    R_iron  = 4 * (RC * 1e9) * Inhib_O2 * (FeOx / (FeOx + Params.KFEMonod));
    R_SRR   = (RC * 1e9) * Inhib_O2 * (SO4 / (SO4 + Params.k_SO4));
    Rate_Meth = 0.5 * (RC * 1e9) * Inhib_O2 * Inhib_Fe * (Params.k_SO4 / (SO4 + Params.k_SO4));
    % 5. Secondary Redox Reactions [umol/L/yr]
    R_AOM         = Params.k_AOM * CH4 * (SO4 / (SO4 + Params.K_CH4_SO4));
    R_aerobic_CH4 = Params.k_aerobic_CH4 * CH4 * (O2 / (O2 + Params.K_CH4_O2));
    R_HS_ox       = Params.Kreox * HS * O2;
    R_Fe_ox       = Params.kFeOx * Fe * O2;
    R_FeS         = Params.kFeS * Fe * HS;
    % 6. Assemble ODEs (dYdx)
    dYdx = zeros(22, 1);
    % --- Aqueous Species (Transport: Diffusion + Advection + Bioirrigation) ---
    % O2
    dYdx(1) = dO2dx;
    dYdx(2) = (v_burial_f * dO2dx - alpha_bio*(Config.O2init - O2) + R_respi + 2*R_aerobic_CH4 + 2*R_HS_ox + 0.25*R_Fe_ox) / (phi * D_O2_eff);
    % SO4
    dYdx(3) = dSO4dx;
    dYdx(4) = (v_burial_f * dSO4dx - alpha_bio*(Config.SO4init - SO4) + 0.5*R_SRR + R_AOM - R_HS_ox) / (phi * D_SO4_eff);
    % CH4
    dYdx(5) = dCH4dx;
    dYdx(6) = (v_burial_f * dCH4dx - alpha_bio*(Config.CH4init - CH4) - Rate_Meth + R_aerobic_CH4 + R_AOM) / (phi * D_CH4_eff);
    % DIC
    dYdx(7) = dDICdx;
    dYdx(8) = (v_burial_f * dDICdx - alpha_bio*(Config.DICinit - DIC) - (RC*1e9 - Rate_Meth) + R1_carb_total) / (phi * D_DIC_eff);
    % ALK
    dYdx(9) = dALKdx;
    dYdx(10)= (v_burial_f * dALKdx - alpha_bio*(Config.HCO3init - ALK) - 2*R_FeS + 2*R1_carb_total) / (phi * D_DIC_eff);
    % Ca
    dYdx(11)= dCadx;
    dYdx(12)= (v_burial_f * dCadx - alpha_bio*(Config.Calcium - Ca) + R1_carb_total) / (phi * D_DIC_eff);
    % Fe2+
    dYdx(13)= dFedx;
    dYdx(14)= (v_burial_f * dFedx - alpha_bio*(Config.Feinit - Fe) - R_iron + R_Fe_ox + R_FeS) / (phi * D_HS_eff);
    % HS-
    dYdx(15)= dHSdx;
    dYdx(16)= (v_burial_f * dHSdx - alpha_bio*(Config.HSinit - HS) - 0.5*R_SRR + R_HS_ox + R_FeS) / (phi * D_HS_eff);
    % --- Solid Phase Species (Transport: Bioturbation + Burial, NO molecular diffusion) ---
    % Solid mass balance uses Db and v_burial_s. Conversions apply.
    % C_organic (g/gDw)
    dYdx(17) = dCorgdx;
    % Governing: d/dx(Db * dC/dx) - v_s * dC/dx - Rate = 0
    dYdx(18) = (v_burial_s * dCorgdx + Temp_factor * k_sed * C_org) / Db;
    % Fe(OH)3 (umol/g)
    dYdx(19) = dFeOxdx;
    R_iron_solid_unit = - R_iron * (phi/(1-phi)) * 1e-3 * (1/Params.rho); % convert aqueous rate to solid unit
    R_FeOx_solid_unit = R_Fe_ox * (phi/(1-phi)) * 1e-3 * (1/Params.rho);
    dYdx(20) = (v_burial_s * dFeOxdx - R_iron_solid_unit - R_FeOx_solid_unit) / Db;
%     % CaCO3 (solid)
%     dYdx(21) = dCaCO3dx;
%     NR_CaCO3 = R_carb_form - R_carb_disso;
%     dYdx(22) = (NR_CaCO3) / v_burial_s; % Assuming no bioturbation for CaCO3 to match original code logic
% --- Solid Phase: CaCO3 ---
    dYdx(21) = dCaCO3dx;
    NR_CaCO3 = R_carb_form - R_carb_disso;
    % Apply Bioturbation (Db) to maintain 2nd-order ODE structural consistency
    dYdx(22) = (v_burial_s * dCaCO3dx - NR_CaCO3) / Db;
end
```

## File: Params.m
```matlab
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
```

