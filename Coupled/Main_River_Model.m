% =========================================================================
% Script: Main_River_Model.m 
% Purpose: Fully coupled 1-D Reactive Transport Model for River Sediments
% =========================================================================
clear; clc; close all;

% --- 1. Load Configurations and Parameters ---
% Replacing legacy 'global' variables with structured data arrays
Config = Config();
Params = Params();

% --- 2. Define Spatial Grid ---
x_mesh = linspace(0, Config.Lbottom, Config.n_grid);

% --- 3. Establish the Index Map (CRITICAL) ---
% bvp4c requires converting 2nd-order ODEs into a system of 1st-order ODEs.
% Y(odd)  = Concentration of species (uM)
% Y(even) = Concentration gradient dC/dx
%
% Species Order: 
% 1: O2, 2: SO4, 3: CH4, 4: DIC, 5: HCO3(Alk), 6: Ca, 7: Fe, 8: HS
num_species = 8;
total_vars  = num_species * 2; 

% --- 4. Initialize bvp4c Solver ---
% Flat initial guess for all 16 variables (8 species + 8 gradients)
initial_guess = zeros(1, total_vars); 

% Set surface boundary conditions as initial guess to help convergence
initial_guess(1)  = Config.O2init;
initial_guess(3)  = Config.SO4init;
initial_guess(5)  = Config.CH4init;
initial_guess(7)  = Config.DICinit;
initial_guess(9)  = Config.HCO3init;
initial_guess(11) = Config.Calcium;
initial_guess(13) = Config.Feinit;
initial_guess(15) = Config.HSinit;

solinit = bvpinit(x_mesh, initial_guess);

% Tight tolerances required for stiff biogeochemical reaction fronts
options = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NMax', 5000);

% --- 5. Execute Fully Coupled Solver ---
fprintf('Running coupled BVP solver...\n');
sol = bvp4c(@(x, Y) ODE_River_System(x, Y, Config, Params), ...
            @(Ya, Yb) BC_River_System(Ya, Yb, Config), ...
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

% --- 8. Original Visualization Logic ---
% [!!! PASTE YOUR ORIGINAL PLOTTING CODE HERE !!!]
% You now have all variables available: 
% depth_cm, O2_prof, SO4_prof, CH4_prof, DIC_prof, ALK_prof, Ca_prof, Fe_prof, HS_prof, pH_prof, Omega_prof