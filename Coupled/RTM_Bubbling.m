function [R_Bubble, Bubble] = RTM_Bubbling(CH4, z, T, Config)
% RTM_BUBBLING
% Simple dissolved-CH4 ebullition loss.
%
% This module does not track gas-phase storage or bubble transport.
% It removes dissolved CH4 only when porewater CH4 exceeds an effective
% pressure/solubility-based threshold.
%
% Inputs:
%   CH4    uM, dissolved methane profile
%   z      cm, sediment depth below sediment-water interface
%   T      deg C
%   Config bubbling parameters
%
% Outputs:
%   R_Bubble   uM/yr, bulk-volume CH4 loss rate
%   Bubble     diagnostics

CH4 = max(CH4(:), 0);
z = z(:);

% --- physical constants ---
P_atm_Pa = 101325;
rho_w = 1000;          % kg/m3
g = 9.80665;           % m/s2
sigma_w = 0.072;       % N/m, water-air surface tension approximation

% --- methane solubility approximation ---
% C_ref is approximate dissolved CH4 concentration in equilibrium with
% pure CH4 gas at 1 atm and 25 C.
% This is intentionally simple for the sediment RTM core.
C_ref_25C = 1400;      % uM per atm CH4
beta_T = 1700;         % K, van't Hoff-like temperature sensitivity

T_K = T + 273.15;
C_CH4_1atm = C_ref_25C .* exp(beta_T .* (1 ./ T_K - 1 ./ 298.15));
C_CH4_1atm = C_CH4_1atm .* ones(size(CH4));

% --- pressure correction ---
z_m = max(z, 0) ./ 100;
water_depth_m = Config.water_depth_m;
bubble_radius_m = max(Config.bubble_radius_m, 1e-6);

P_hydro_atm = rho_w .* g .* (water_depth_m + z_m) ./ P_atm_Pa;
P_laplace_atm = 2 .* sigma_w ./ bubble_radius_m ./ P_atm_Pa;

P_total_atm = 1 + P_hydro_atm + P_laplace_atm;

% Effective CH4 partial pressure in bubbles.
% This approximates mixed-gas bubbles without tracking N2/O2/CO2 explicitly.
CH4_threshold = Config.bubble_CH4_fraction .* C_CH4_1atm .* P_total_atm;

if Config.use_CH4_bubbling
    R_Bubble = Config.k_bubble .* max(CH4 - CH4_threshold, 0);
else
    R_Bubble = zeros(size(CH4));
end

Bubble = struct();
Bubble.CH4_threshold = CH4_threshold;
Bubble.P_total_atm = P_total_atm;
Bubble.P_hydro_atm = P_hydro_atm;
Bubble.P_laplace_atm = P_laplace_atm;
Bubble.C_CH4_1atm = C_CH4_1atm;
Bubble.use_CH4_bubbling = Config.use_CH4_bubbling;
Bubble.k_bubble = Config.k_bubble;
Bubble.bubble_CH4_fraction = Config.bubble_CH4_fraction;
end