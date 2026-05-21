# MATLAB Model Source Code

## File: Build_Forcing_PDE.m
```matlab
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
```

## File: Build_Grid_1D.m
```matlab
function Grid = Build_Grid_1D(Config, Params)
% BUILD_GRID_1D
% Build the 1-D sediment grid and depth-dependent physical fields.
Grid.z = linspace(0, Config.Lbottom, Config.n).';
Grid.n = numel(Grid.z);
Grid.dz = Grid.z(2) - Grid.z(1);
z = Grid.z;
Grid.poros = Config.porosbottom + ...
    (Config.porostop - Config.porosbottom) .* exp(-z ./ Config.porosscale);
Grid.D_solid_mix = Config.Bioturbbottom + ...
    (Config.Bioturbtop - Config.Bioturbbottom) .* exp(-z ./ Config.bioturbscale);
Grid.Alpha_exchange = Config.Bioirrig_bottom + ...
    (Config.Bioirrig_top - Config.Bioirrig_bottom) .* exp(-z ./ Config.Bioirrig_scale);
Grid.v_solid = Config.vbottom .* ...
    (1 - Config.porosbottom) ./ max(1 - Grid.poros, 1e-6);
Grid.v_fluid = Config.vbottom_fluid .* ...
    (1 + Config.porosbottom) ./ max(1 + Grid.poros, 1e-6);
age = Config.ageinit + cumsum(Grid.dz ./ max(Grid.v_solid, 1e-12));
Grid.k_sed = 10.^(-0.95 .* log10(age) - 0.81);
Grid.rho = Params.rho;
end
```

## File: CH4_bc.m
```matlab

function res = CH4_bc(CH4a,CH4b)
global CH4init
  res = [ CH4a(1)-CH4init
          CH4b(2) ];
end
```

## File: CH4_ODE.m
```matlab

function dydx = CH4_ODE(x,CH4)
global k_SO4 RC Oxygen k_O2 DCH4 v_burial_Fluid Alpha_Bioirrig z_sed poros Sulfate CH4init
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2
global Rate_Meth
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
% SO4 = interp1(z_sed,Sulfate,x);
SO4_raw = interp1(z_sed, Sulfate, x, 'linear', 'extrap');
SO4_pos = max(real(SO4_raw), 0);
if SO4_pos <= 1e-6
    f_AOM_SO4 = 0;
else
    f_AOM_SO4 = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
end
% Inh = (k_O2./(O2+k_O2));
% Inh1 = (k_SO4./(SO4+k_SO4));
% RC1 = interp1(z_sed,RC,x);
R_Meth_current = double(interp1(z_sed, Rate_Meth, x));
% NR = + v_burial_f.* CH4(2) - 0.5.*RC1.*Inh.*Inh1.* 1E9 - (Alpha_Bioirrig_1.*(CH4init-CH4(1)))...
%       + k_AOM.* CH4(1).* (SO4./(SO4+K_CH4_SO4)) + k_aerobic_CH4.* CH4(1).* (O2./(O2+K_CH4_O2)); % umol/l/year
NR = + v_burial_f.* (CH4(2)/(fi*DCH4)) - R_Meth_current - (Alpha_Bioirrig_1.*(CH4init-CH4(1)))...
    + k_AOM .* CH4(1) .* f_AOM_SO4 ...%       + k_AOM.* CH4(1).* (SO4./(SO4+K_CH4_SO4))...
    + k_aerobic_CH4.* CH4(1).* (O2./(O2+K_CH4_O2)); % umol/l/year
dydx = [ CH4(2) /fi/DCH4
           NR];
end
```

## File: Config_Baseline.m
```matlab
function Config = Config_Baseline()
% CONFIG_BASELINE
% Site / scenario specific settings for the OLD sequential core.
Config.Corg_top = 0.02;    % g/gDw, i.e. 1.2 % dry weight at sediment surface
    % ---------------- Domain ----------------
    Config.Lbottom = 30;          % cm
    Config.n = 101;
    Config.nmesh = 1000;
    % ---------------- Physical structure ----------------
    Config.vbottom = 0.2;         % cm / yr   (river-audited baseline)
    Config.vbottom_fluid = 0;     % cm / yr
    Config.porostop = 0.9;
    Config.porosbottom = 0.7;
    Config.porosscale = 3;
    Config.Bioturbtop = 1;%1;%10;       % cm2 / yr
    Config.Bioturbbottom = 0.05;%0.05;%1;     % cm2 / yr
    Config.bioturbscale = 3;
    Config.Bioirrig_top = 20;%20;%100;    % 1 / yr
    Config.Bioirrig_bottom = 0;
    Config.Bioirrig_scale =1;%1.0;%0.75;
    Config.f_lab = 0.6;   % fraction of total OM input entering labile/reactive pool
    Config.k_ref_factor = 0;   % first test: k_ref = 0.03 * k_sed
    % ---------------- Boundary concentrations ----------------
    Config.O2init   = 250;%150;        % uM
    Config.SO4init  = 200;%200;        % uM
    Config.DICinit  = 1200;%1200;%3400;       % uM
    Config.HCO3init = 1100;%1100;%3500;        % uM
    Config.Calcium  = 1000;%1000;%5700;       % uM
    Config.CH4init  = 0;          % uM
    Config.Feinit   = 0;          % uM
    Config.HSinit   = 0;          % uM
    Config.Pinitial = 0;          % uM
    % ---------------- Fluxes ----------------
    Config.NPP = 400;             % g / m2 / yr
    Config.BE  = 0.1;
    Config.F_FeOx  = 1.5;%2;          % mmol / m2 / d , external reactive Fe input forcing
    Config.F_CaCO3 = 2;%;          % g / m2 / yr
    % ---------------- Temperature / OM age ----------------
    Config.T_future = 25;%25;
    Config.ageinit = 0.1;
    Config.age_root = 1;
    Config.Salinity = 0.1;%0.1;
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
    % ---------------- Transient PDE settings ----------------
    Config.PDE_use_steady_IC = true;   % development mode only
    Config.PDE_do_spinup     = true;   % formal transient runs should spin up
    Config.t_spinup = 100;             % yr, first candidate
    Config.t_final  = 3;               % yr, transient experiment length after spin-up
    Config.dt_out_spinup = 1;          % yr
    Config.dt_out_event  = 1/365;      % yr, daily output
    Config.n_pde = Config.n;
end
```

## File: Coupled_Carbonate_bc.m
```matlab
function res = Coupled_Carbonate_bc(Ya, Yb)
global DICinit HCO3init CaCO3_init
  res = [ Ya(1) - DICinit;       % 1: DIC top
          Yb(2);                 % 2: DIC bottom flux = 0
          Ya(3) - HCO3init;      % 3: ALK top
          Yb(4);                 % 4: ALK bottom flux = 0
          Ya(5) - CaCO3_init;    % 5: CaCO3 top
          Yb(6) ];               % 6: CaCO3 bottom flux = 0
end
% function res = Coupled_Carbonate_bc(Ya, Yb)
% global DICinit HCO3init Calcium CaCO3_init
%
% % State order:
% % 1 DIC, 2 DIC flux
% % 3 ALK, 4 ALK flux
% % 5 Ca,  6 Ca flux
% % 7 CaCO3, 8 CaCO3 dummy flux
%
% res = [ Ya(1) - DICinit;     % DIC top
%         Yb(2);               % DIC bottom dissolved flux = 0
%
%         Ya(3) - HCO3init;    % ALK top (keep your current boundary for now)
%         Yb(4);               % ALK bottom dissolved flux = 0
%
%         Ya(5) - Calcium;     % dissolved Ca top
%         Yb(6);               % dissolved Ca bottom flux = 0
%
%         Ya(7) - CaCO3_init;  % solid CaCO3 top
%         Yb(8) ];             % solid CaCO3 bottom dummy flux = 0
% end
```

## File: Coupled_Carbonate_ODE.m
```matlab
function dYdx = Coupled_Carbonate_ODE(x, Y)
% State order:
% Y(1) = DIC,   Y(2) = DIC flux variable
% Y(3) = ALK,   Y(4) = ALK flux variable
% Y(5) = CaCO3, Y(6) = CaCO3 dummy flux variable
global DHCO3 DICinit HCO3init
global Alpha_Bioirrig v_burial_Fluid v_burial z_sed poros rho
global k_calcite k_calcite_dis1 k_calcite_dis2
global n_power_CaCO31 n_power_CaCO32 n_power_CaCO33
global Calcium Calcium_activity CO3_activity Ksp_ca
global T_future Salinity
global R_DIC_prod R_ALK_prod
v_burial_f = double(interp1(z_sed, v_burial_Fluid, x, 'linear', 'extrap'));
v_burial_s = double(interp1(z_sed, v_burial,       x, 'linear', 'extrap'));
fi         = double(interp1(z_sed, poros,          x, 'linear', 'extrap'));
Alpha      = double(interp1(z_sed, Alpha_Bioirrig, x, 'linear', 'extrap'));
R_DIC_prod_1 = double(interp1(z_sed, R_DIC_prod, x, 'linear', 'extrap'));
R_ALK_prod_1 = double(interp1(z_sed, R_ALK_prod, x, 'linear', 'extrap'));
DIC   = max(real(Y(1)), 1e-12);
ALK   = max(real(Y(3)), 1e-12);
CaCO3 = max(real(Y(5)), 0);
[~, CO3_current, ~] = River_Carbonate(ALK, DIC, T_future, Salinity, 1);
CO3_current = max(real(CO3_current), 1e-12);
sigma_carb = (Calcium_activity * Calcium * CO3_activity * CO3_current) / Ksp_ca - 1;
unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho * max(1 - fi, 1e-6));
if sigma_carb > 0
    R_carb_form = abs(sigma_carb)^n_power_CaCO31 * k_calcite * unit_conversion;
else
    R_carb_form = 0;
end
if (sigma_carb < 0) && (sigma_carb > -0.2)
    R_carb_disso = abs(sigma_carb)^n_power_CaCO32 * k_calcite_dis1 * CaCO3;
elseif sigma_carb <= -0.2
    R_carb_disso = abs(sigma_carb)^n_power_CaCO33 * k_calcite_dis2 * CaCO3;
else
    R_carb_disso = 0;
end
R1_carb_total = R_carb_form - R_carb_disso * (1E3 * 1E6 * 1E-2 * rho * max(1 - fi, 1e-6));
Advection_DIC = v_burial_f .* (Y(2) / (fi * DHCO3));
NR_DIC = Advection_DIC ...
       - R_DIC_prod_1 ...
       + R1_carb_total ...
       - Alpha .* (DICinit - DIC);
Advection_ALK = v_burial_f .* (Y(4) / (fi * DHCO3));
NR_ALK = Advection_ALK ...
       - R_ALK_prod_1 ...
       + 2 * R1_carb_total ...
       - Alpha .* (HCO3init - ALK);
NR_CaCO3 = R_carb_form - R_carb_disso;
dYdx = zeros(6,1);
dYdx(1) = Y(2) / (fi * DHCO3);
dYdx(2) = NR_DIC;
dYdx(3) = Y(4) / (fi * DHCO3);
dYdx(4) = NR_ALK;
dYdx(5) = NR_CaCO3 / max(v_burial_s, 1e-8);
dYdx(6) = 0;
end
```

## File: Fe3_bc.m
```matlab
function res = Fe3_bc(Fe3a,Fe3b)
global Fe_3_init
  res = [ Fe3a(1)-Fe_3_init
          Fe3b(2) ];
end

```

## File: Fe3_ODE.m
```matlab

function dydx = Fe3_ODE(x,Fe3)
global z_sed Bioturb R_FeRed
global rho Oxygen v_burial C_Fe kFeOx poros
R_FeRed_1  = interp1(z_sed, R_FeRed, x, 'linear', 'extrap');   % umol/L/yr
poros_1    = interp1(z_sed, poros, x, 'linear', 'extrap');
O2         = interp1(z_sed, Oxygen, x, 'linear', 'extrap');
C_Fe_1     = interp1(z_sed, C_Fe, x, 'linear', 'extrap');
Db         = max(interp1(z_sed, Bioturb, x, 'linear', 'extrap'), 1e-6);
v_burial_1 = max(interp1(z_sed, v_burial, x, 'linear', 'extrap'), 1e-6);
sigh       = max(1 - poros_1, 1e-6);
R_FeRed_solid = R_FeRed_1 .* (poros_1 ./ max(1 - poros_1, 1e-6)) .* 1e-3 ./ rho;
R_FeOx_solid  = (kFeOx .* C_Fe_1 .* O2) .* (poros_1 ./ max(1 - poros_1, 1e-6)) .* 1e-3 ./ rho;
NR = - R_FeRed_solid + R_FeOx_solid;   % umol/g/yr
dFe3dx = Fe3(2) ./ (sigh .* Db);
dydx = [ dFe3dx
         v_burial_1 .* dFe3dx + NR ];
end

```

## File: Fe_bc.m
```matlab
function res = Fe_bc(Fea,Feb)
global Feinit
  res = [ Fea(1)-Feinit
          Feb(2) ];
end
```

## File: Fe_ODE.m
```matlab
function dydx = Fe_ODE(x,Fe2)
global Oxygen DH2S v_burial_Fluid Alpha_Bioirrig z_sed poros Feinit
global kFeOx R_FeRed kFeS C_HS pH K_HS
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
C_HS_tot = interp1(z_sed, C_HS, x);
pH_x = interp1(z_sed, pH, x);
C_HS_free = C_HS_tot ./ (1 + ((10.^(6 - pH_x)) ./ K_HS));
R_FeRed_1 = interp1(z_sed,R_FeRed,x);
NR = + v_burial_f .* (Fe2(2)/(fi*DH2S)) ...
     - R_FeRed_1 ...
     - (Alpha_Bioirrig_1 .* (Feinit-Fe2(1))) ...
+ (kFeS .* Fe2(1) .* C_HS_free) ...
     + (kFeOx .* Fe2(1) .* O2);   % umol/l/year
dydx = [ Fe2(2) /fi/DH2S
         NR ];
end
```

## File: HS_bc.m
```matlab

function res = HS_bc(HSa,HSb)
global HSinit
  res = [ HSa(1)-HSinit
          HSb(2) ];
end
```

## File: HS_ODE.m
```matlab

function dydx = HS_ODE(x, H2S)
global Oxygen v_burial_Fluid HSinit
global Kreox Alpha_Bioirrig z_sed poros DH2S C_Fe kFeS pH K_HS
global R_SRR R_AOM_actual
v_burial_f = interp1(z_sed, v_burial_Fluid, x, 'linear', 'extrap');
Alpha_Bioirrig_1 = interp1(z_sed, Alpha_Bioirrig, x, 'linear', 'extrap');
fi = interp1(z_sed, poros, x, 'linear', 'extrap');
O2 = interp1(z_sed, Oxygen, x, 'linear', 'extrap');
C_Fe_1 = interp1(z_sed, C_Fe, x, 'linear', 'extrap');
pH_x = interp1(z_sed, pH, x, 'linear', 'extrap');
R_SRR_1 = interp1(z_sed, R_SRR, x, 'linear', 'extrap');
R_AOM_1 = interp1(z_sed, R_AOM_actual, x, 'linear', 'extrap');
HS_free = H2S(1) ./ (1 + ((10.^(6 - pH_x)) ./ K_HS));
NR = + v_burial_f .* (H2S(2)/(fi*DH2S)) ...
     - R_SRR_1 ...
     - R_AOM_1 ...
     + (kFeS .* C_Fe_1 .* HS_free) ...
     + (Kreox .* HS_free .* O2) ...
     - (Alpha_Bioirrig_1 .* (HSinit - H2S(1)));
dydx = [ H2S(2) / fi / DH2S
         NR ];
end
```

## File: Hydro_Preprocessor.m
```matlab
function Hydro = Hydro_Preprocessor(Config, Params)
% HYDRO_PREPROCESSOR
% Output only weak-coupling modifiers for now.
    % ---------- 1. Packing-based porosity predictor ----------
    X_f   = Config.X_f;
    phi_c = Config.phi_c;
    phi_f = Config.phi_f;
    if X_f < phi_c
        phi_mix = phi_c - X_f * (1 - phi_f);
    else
        phi_mix = X_f * phi_f;
    end
    Hydro.phi_mix = phi_mix;
    % ---------- 2. Hydrodynamics ----------
    g = 9.81;
    rho_w = 1000;
    u_star = sqrt(g * Config.H * Config.S);   % m / s
    tau_b  = rho_w * u_star^2;                % N / m2
    Hydro.u_star = u_star;
    Hydro.tau_b = tau_b;
    % ---------- 3. Weak POM multiplier ----------
    % For old-core sanity runs, do NOT hard-zero NPP.
    if tau_b >= Config.tau_c
        raw_mult = 0.25;   % keep nonzero to avoid killing OM input
    else
        raw_mult = 1 - tau_b / Config.tau_c;
    end
    Hydro.NPP_multiplier = max(0.25, min(1.0, raw_mult));
    % ---------- 4. Weak diffusion multiplier ----------
    attenuation_factor = 1e-3;
    v_pore_m  = u_star * attenuation_factor;
    v_pore_cm = v_pore_m * 100;
    seconds_per_year = 3600 * 24 * 365;
    D_hyp_surface = Config.alpha_L * v_pore_cm * seconds_per_year;
    % Convert to a bounded multiplier relative to old molecular D ~ 300-400
    D_ref = mean([Params.DO2, Params.DSO4, Params.DH2S, Params.DCH4, Params.DHCO3]);
    raw_mult = 1 + D_hyp_surface / D_ref;
    Hydro.diffusion_multiplier = min(raw_mult, Config.max_diffusion_multiplier);
    Hydro.D_hyp_surface = D_hyp_surface;
end
```

## File: Initialize_PDE_State.m
```matlab
function State = Initialize_PDE_State(Grid, Forcing, Config, Params)
% INITIALIZE_PDE_STATE
% Simple positive initial condition for PDE spin-up.
%
% This is not a formal steady state. It only gives ode15s a reasonable
% nonnegative starting point.
z = Grid.z;
n = Grid.n;
rho = Params.rho;
A_solid = rho .* max(1 - Grid.poros, 1e-6);
solid_flux = A_solid .* max(Grid.v_solid, 1e-12);
% Labile OM: depositional top value with first-order burial decay.
OM_lab_top = Forcing.F_lab_OM ./ max(solid_flux(1), 1e-12);
decay_int = cumsum((Grid.k_sed ./ max(Grid.v_solid, 1e-12)) .* Grid.dz);
State.OM_lab = max(OM_lab_top .* exp(-decay_int), 1e-12);
% Refractory OM: preserved inventory approximation.
State.OM_ref = max(Forcing.F_ref_OM ./ max(solid_flux, 1e-12), 0);
% FeOOH top inventory from external Fe flux.
% F_FeOx: mmol/m2/d -> umol/cm2/yr by *36.5.
F_FeOOH = Forcing.F_FeOx .* 36.5;
FeOOH_top = F_FeOOH ./ max(solid_flux(1), 1e-12);
if isfield(Params, 'Fe_inventory_factor')
    FeOOH_top = Params.Fe_inventory_factor .* FeOOH_top;
end
State.FeOOH = max(FeOOH_top .* exp(-z ./ 3), 1e-12);
% CaCO3 solid top inventory from depositional flux.
F_CaCO3_cm2 = Forcing.F_CaCO3 .* 1e-4;  % g/m2/yr -> g/cm2/yr
CaCO3_top = F_CaCO3_cm2 ./ max(solid_flux(1), 1e-12);
State.CaCO3 = max(CaCO3_top .* exp(-z ./ 10), 0);
% Solutes.
State.O2  = max(Forcing.O2_top .* exp(-z ./ 1.5), 0);
State.Fe2 = zeros(n,1);
State.SO4 = max(Forcing.SO4_top .* ones(n,1), 0);
State.HS  = zeros(n,1);
State.CH4 = max(Forcing.CH4_top .* ones(n,1), 0);
State.DIC = max(Forcing.DIC_top .* ones(n,1), 1e-12);
State.ALK = max(Forcing.ALK_top .* ones(n,1), 1e-12);
end
```

## File: O2_bc.m
```matlab

function res = O2_bc(O2a,O2b)
global O2init
  res = [ O2a(1)-O2init
          O2b(2) ];
end
```

## File: O2_ODE.m
```matlab

function dydx = O2_ODE(x,O2)
global DO2 k_O2 z_sed RC Alpha_Bioirrig O2init v_burial_Fluid poros O2_root
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2_root_ODE = interp1(z_sed,O2_root,x);
RC1 = interp1(z_sed,RC,x);
% NR = + v_burial_f.* O2(2) + RC1 * (O2(1)/(O2(1)+k_O2)) * 1E9 - (Alpha_Bioirrig_1.*(O2init-O2(1))) - O2_root_ODE; % umol/l/year
NR = + v_burial_f.* (O2(2)/(fi*DO2)) + RC1 * (O2(1)/(O2(1)+k_O2)) * 1E9 - (Alpha_Bioirrig_1.*(O2init-O2(1))) - O2_root_ODE; % umol/l/year
dydx = [ O2(2) /fi/DO2
           NR];
end

```

## File: organicbc.m
```matlab
function res = organicbc(C_orga,C_orgb)
global F_lab_OM v_burial poros rho Bioturb
NPP1 = F_lab_OM; %gram/cm2/year
v_burial1 = v_burial(1,1);  %cm/year
poros1 = poros(1,1);
A1 = rho * (1-poros1);
Bioturb1 = Bioturb(1,1);
BC_1 = - Bioturb1 * A1 * C_orga(2) + A1 * v_burial1 * C_orga(1) - NPP1;
res = [ BC_1
        C_orgb(2)];
end
% function res = organicbc(C_orga,C_orgb)
% global Corg_top
% res = [ C_orga(1) - Corg_top     % top: fixed solid-phase OM concentration
%         C_orgb(2) ];             % bottom: zero gradient / zero diffusive flux
% end
```

## File: organicbc_1.m
```matlab
function res = organicbc_1(C_orga,C_orgb)
global NPP v_burial poros rho Bioturb BE
NPP1 = BE * NPP * 1E-4; %gram/cm2/year
v_burial1 = v_burial(1,1);  %cm/year
poros1 = poros(1,1);
A1 = rho * (1-poros1);
Bioturb1 = Bioturb(1,1);
BC_1 = - Bioturb1 * A1 * C_orga(2) + A1 * v_burial1 * C_orga(1) - NPP1;
res = [ BC_1
        C_orgb(2)];
end
% function res = organicbc_1(C_orga,C_orgb)
% global Corg_top
% res = [ C_orga(1) - Corg_top     % top: fixed solid-phase OM concentration
%         C_orgb(2) ];             % bottom: zero gradient / zero diffusive flux
% end
```

## File: organicODE.m
```matlab
function dydx = organicODE(x,C_org)
global k_sed v_burial z_sed Bioturb poros Temp_factor
v_burial_1 = interp1(z_sed,v_burial,x);
k_sed1 = interp1(z_sed,k_sed,x);
Db = interp1(z_sed,Bioturb,x);
fi = interp1(z_sed,poros,x);
sigh = 1 - fi;
% RC = k_sed1.*C_org(1); %k_sed1.*u.*rho.*((1-phi)/phi)*12; % molCorg/cm3sed/yr mineralization rate
NR = + v_burial_1.* (C_org(2)/sigh) + Temp_factor.*k_sed1.*C_org(1);  %g/gDw/year
% NR = + v_burial_1.* (C_org(2)) + Temp_factor.*k_sed1.*C_org(1);  %g/gDw/year
dydx = [ C_org(2) /sigh/Db
         NR];
end
```

## File: organicODE_1.m
```matlab
function dydx = organicODE_1(x,C_org)
global k_sed v_burial z_sed Bioturb Temp_factor
v_burial_1 = interp1(z_sed,v_burial,x);
k_sed1 = interp1(z_sed,k_sed,x);
Db = interp1(z_sed,Bioturb,x);
RC = Temp_factor.*k_sed1.*C_org(1); %k_sed1.*u.*rho.*((1-phi)/phi)*12; % molCorg/cm3sed/yr mineralization rate
NR = - RC;  %g/gDw/year
dydx = [ NR / v_burial_1
         0];
end
```

## File: Pack_State.m
```matlab
function Y = Pack_State(State)
% PACK_STATE
% Convert state struct to one column vector for ode15s.
Y = [
    State.OM_lab(:)
    State.OM_ref(:)
    State.FeOOH(:)
    State.CaCO3(:)
    State.O2(:)
    State.Fe2(:)
    State.SO4(:)
    State.HS(:)
    State.CH4(:)
    State.DIC(:)
    State.ALK(:)
];
end
```

## File: Params_Static.m
```matlab
function Params = Params_Static()
% PARAMS_STATIC
% Constants that should not change from site to site unless explicitly audited.
    Params.rho = 2.73;           % g / cm3
    Params.k_O2 = 2;             % uM
    Params.k_SO4 = 20;           % uM
    Params.KFEMonod = 1000;%200;       % umol / g
    Params.DSO4 = 150;%300;           % cm2 / yr
    Params.DCH4 = 300;           % cm2 / yr
    Params.DH2S = 300;           % cm2 / yr
    Params.DO2  = 300;           % cm2 / yr
    Params.DHCO3 = 400;          % cm2 / yr
    Params.DPO4  = 400;          % cm2 / yr
    Params.Kreox = 500;          % 1 / umol / L / yr
    Params.kFeOx = 10;%10;           % 1 / umol / L / yr
    Params.kFeS  = 1;%1;%10;           % 1 / umol / L / yr
    Params.K_CH4_SO4   = 100;    % uM
    Params.K_CH4_O2    = 1;      % uM
    Params.k_AOM       = 0.2;%0.2;%1.0;    % 1 / yr
    Params.k_aerobic_CH4 = 6;    % 1 / yr
    Params.Ksp_ca = 4.5e5;%3000;        % uM^2
    Params.k_calcite = 1;
    Params.k_calcite_dis1 = 0.005;
    Params.k_calcite_dis2 = 10;
    Params.n_power_CaCO31 = 1.76;
    Params.n_power_CaCO32 = 0.11;
    Params.n_power_CaCO33 = 4;
    Params.Calcium_activity = 0.6;%1.0;%0.6;
    Params.CO3_activity     = 0.6;%1.0;%0.6;
    Params.P_C_ratio = 0.0094;
    Params.Q10   = 2;
    Params.T_ref = 25;
    Params.FeC_frac_max = 0.35;
    Params.Fe_inventory_factor = 0.4;  % global fixed factor, first candidate
    % extras already used by old core
%     Params.KFeS = 2500;
    Params.K_HS = 7;
    Params.kapatite = 0.05;
    Params.Kviv = 3e6;
    Params.alpha_viv = 1.5;
    Params.kviv = 1.7e-22;
end
```

## File: Phos_bc.m
```matlab
function res = Phos_bc(phos1a,phos1b)
global Pinitial
  res = [ phos1a(1)-Pinitial
          phos1b(2) ];
end
```

## File: Phos_ODE.m
```matlab
function dydx = Phos_ODE(x,p1)
global poros z_sed P_C_ratio
global RC   %%molCorg/cm3sed/yr
global DPO4
global Rviv1 Rapat
% y(1) is Sulfide concentration in uM
  fi = interp1(z_sed,poros,x);
  Rcarbon = interp1(z_sed,RC,x);  %molCorg/cm3sed/yr
  Rviv = interp1(z_sed,Rviv1,x);
  Rapa = interp1(z_sed,Rapat,x);
  prate = Rapa + Rviv - (P_C_ratio.*Rcarbon.*1E9); %umolS/Lsed/yr
  dydx = [ p1(2) /fi/DPO4
           prate];
end
```

## File: Plot_PDE_Result.m
```matlab
function Plot_PDE_Result(Result)
% PLOT_PDE_RESULT
% Basic diagnostic plot for the first PDE spin-up.
Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;
z = Grid.z;
figure('Name', 'FV-MOL PDE Result', 'Color', 'w');
clf;
n_plot = 6;
m_plot = 3;
subplot(m_plot,n_plot,1);
plot((S.OM_lab + S.OM_ref) .* 100, z, 'LineWidth', 2); axis ij
title('Organic (%gDW)');
ylabel('Depth (cm)');
grid on; box on
subplot(m_plot,n_plot,2);
plot(S.O2, z, 'LineWidth', 2); axis ij
title('O_2 (\muM)');
grid on; box on
subplot(m_plot,n_plot,3);
plot(S.Fe2, z, 'LineWidth', 2); axis ij
title('Fe^{2+} (\muM)');
grid on; box on
subplot(m_plot,n_plot,4);
plot(S.SO4, z, 'LineWidth', 2); axis ij
title('SO_4 (\muM)');
grid on; box on
subplot(m_plot,n_plot,5);
plot(S.HS, z, 'LineWidth', 2); axis ij
title('H_2S total (\muM)');
grid on; box on
subplot(m_plot,n_plot,6);
plot(S.CH4, z, 'LineWidth', 2); axis ij
title('CH_4 (\muM)');
grid on; box on
subplot(m_plot,n_plot,7);
plot(S.CaCO3 .* 100, z, 'LineWidth', 2); axis ij
title('CaCO_3 (%gDW)');
ylabel('Depth (cm)');
grid on; box on
subplot(m_plot,n_plot,8);
plot(S.DIC, z, 'LineWidth', 2); axis ij
title('DIC (\muM)');
grid on; box on
subplot(m_plot,n_plot,9);
plot(S.ALK, z, 'LineWidth', 2); axis ij
title('ALK (\muM)');
grid on; box on
subplot(m_plot,n_plot,10);
plot(D.H2CO3, z, 'LineWidth', 2); axis ij
title('H_2CO_3 (\muM)');
grid on; box on
subplot(m_plot,n_plot,11);
plot(D.pH, z, 'LineWidth', 2); axis ij
title('pH');
grid on; box on
subplot(m_plot,n_plot,12);
plot(D.RC_uM, z, 'LineWidth', 2); axis ij
title('Mineralization (\muM/yr)');
grid on; box on
subplot(m_plot,n_plot,13);
plot(S.FeOOH, z, 'LineWidth', 2); axis ij
title('FeOOH (\mumol/g)');
ylabel('Depth (cm)');
grid on; box on
subplot(m_plot,n_plot,14);
plot(D.R_SRR ./ 365, z, 'LineWidth', 2); axis ij
title('SRR (\muM/d)');
grid on; box on
subplot(m_plot,n_plot,15);
plot(D.R_Meth ./ 365, z, 'LineWidth', 2); axis ij
title('Methanogenesis (\muM/d)');
grid on; box on
subplot(m_plot,n_plot,16);
plot(D.R_FeS ./ 365, z, 'LineWidth', 2); axis ij
title('FeS formation (\muM/d)');
grid on; box on
subplot(m_plot,n_plot,17);
plot(D.sigma_carb, z, 'LineWidth', 2); axis ij
title('\Omega - 1');
grid on; box on
subplot(m_plot,n_plot,18);
plot(D.R_carb_net_solid, z, 'LineWidth', 2); axis ij
title('Net CaCO_3 rxn (g/g/yr)');
grid on; box on
drawnow;
end
```

## File: River_Carbonate.m
```matlab
function [pH, CO3, H2CO3] = River_Carbonate(ALK, DIC, T, S, P)
    if nargin < 3 || isempty(T), T = 20; end
    if nargin < 4 || isempty(S), S = 0.1; end
    if nargin < 5 || isempty(P), P = 1; end %#ok<NASGU>
    ALK = max(real(ALK), 1e-12);
    DIC = max(real(DIC), 1e-12);
    T_K = T + 273.15;
    B_T = 400 * (S / 35);   % umol/kg
    lnK1 = 2.83655 - 2307.1266 / T_K - 1.5529413 * log(T_K) ...
         - (0.20760841 + 4.0484 / T_K) * sqrt(S) ...
         + 0.08468345 * S - 0.00654208 * S^(1.5) ...
         + log(1 - 0.001005 * S);
    K1 = exp(lnK1) * 1e6;
    lnK2 = -9.226508 - 3351.6106 / T_K - 0.2005743 * log(T_K) ...
         + (-0.106901773 - 23.9722 / T_K) * sqrt(S) ...
         + 0.1130822 * S - 0.00846934 * S^(1.5) ...
         + log(1 - 0.001005 * S);
    K2 = exp(lnK2) * 1e6;
    lnKb = (-8966.90 - 2890.53 * sqrt(S) - 77.942 * S + 1.728 * S^(1.5) - 0.0996 * S^2) / T_K ...
         + 148.0248 + 137.1942 * sqrt(S) + 1.62142 * S ...
         + (-24.4344 - 25.085 * sqrt(S) - 0.2474 * S) * log(T_K) ...
         + 0.053105 * sqrt(S) * T_K;
    Kb = exp(lnKb) * 1e6;
    Kw = exp(148.96502 - 13847.26 / T_K - 23.6521 * log(T_K) ...
       + (118.67 / T_K - 5.977 + 1.0495 * log(T_K)) * sqrt(S) - 0.01615 * S) * 1e12;
    f = @(logH) alk_balance_residual(10.^logH, ALK, DIC, B_T, K1, K2, Kb, Kw);
    % lighter bracket search
    logH_grid = linspace(-7, 2, 81);
    f_grid = zeros(size(logH_grid));
    for i = 1:numel(logH_grid)
        f_grid(i) = f(logH_grid(i));
    end
    idx = find(f_grid(1:end-1) .* f_grid(2:end) <= 0, 1, 'first');
    if ~isempty(idx)
        logH = fzero(f, [logH_grid(idx), logH_grid(idx+1)]);
    else
        [~, imin] = min(abs(f_grid));
        logH = logH_grid(imin);
    end
    H = max(10.^logH, 1e-12);
    pH = 6 - log10(H);
    denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
    denom = max(denom, 1e-30);
    H2CO3 = DIC ./ denom;
    CO3   = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
    H2CO3 = max(real(H2CO3), 1e-12);
    CO3   = max(real(CO3),   1e-12);
end
function res = alk_balance_residual(H, ALK, DIC, B_T, K1, K2, Kb, Kw)
    H = max(H, 1e-12);
    denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
    denom = max(denom, 1e-30);
    HCO3 = DIC .* (K1 ./ H) ./ denom;
    CO3  = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
    BOH4 = B_T .* (Kb ./ H) ./ (1 + Kb ./ H);
    OH   = Kw ./ H;
    TA_calc = HCO3 + 2 .* CO3 + BOH4 + OH - H;
    res = TA_calc - ALK;
end
```

## File: RTM_Budget.m
```matlab
function Budget = RTM_Budget(Result)
% RTM_BUDGET
% PDE-consistent Fe, CH4, and FeOOH budget diagnostics.
%
% Solute budgets use phi-weighted reaction integrals because PDE storage is:
%   d(phi*C)/dt = transport + phi*reaction + phi*exchange
%
% Solid budgets use A-weighted reaction integrals:
%   A = rho*(1-phi)
Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;
P = Result.Params;
F = Result.Forcing;
z = Grid.z(:);
dz = Grid.dz;
phi = Grid.poros(:);
A = Grid.rho .* max(1 - phi, 1e-6);
Budget = struct();
% =========================
% CH4 budget
% =========================
CH4 = S.CH4(:);
J_CH4_top_down = top_solute_flux(CH4, F.CH4_top, phi, P.DCH4, Grid.v_fluid, dz);
J_CH4_top_up = max(0, -J_CH4_top_down);
I_CH4_irrig_sink = trapz(z, phi .* Grid.Alpha_exchange(:) .* max(CH4 - F.CH4_top, 0)) .* 1e-3;
Budget.CH4.I_Meth = trapz(z, D.R_Meth(:)) .* 1e-3;
Budget.CH4.I_AOM  = trapz(z, D.R_AOM(:)) .* 1e-3;
Budget.CH4.I_Ox   = trapz(z, D.R_CH4Ox(:)) .* 1e-3;
Budget.CH4.I_IrrigSink = I_CH4_irrig_sink;
Budget.CH4.J_TopUp = J_CH4_top_up;
Budget.CH4.Source = Budget.CH4.I_Meth;
Budget.CH4.Sink = Budget.CH4.I_AOM + Budget.CH4.I_Ox + ...
                  Budget.CH4.I_IrrigSink + Budget.CH4.J_TopUp;
Budget.CH4.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dCH4dt = (S.CH4(:) - S_prev.CH4(:)) ./ max(dt, 1e-12);
    Budget.CH4.Storage = trapz(z, phi .* dCH4dt) .* 1e-3;
end
Budget.CH4.Residual = Budget.CH4.Source - Budget.CH4.Sink - nan0(Budget.CH4.Storage);
Budget.CH4.SinkSourceRatio = Budget.CH4.Sink ./ max(Budget.CH4.Source, 1e-12);
% =========================
% Fe2 budget
% =========================
Fe2 = S.Fe2(:);
J_Fe2_top_down = top_solute_flux(Fe2, F.Fe_top, phi, P.DH2S, Grid.v_fluid, dz);
J_Fe2_top_up = max(0, -J_Fe2_top_down);
I_Fe2_irrig_sink = trapz(z, phi .* Grid.Alpha_exchange(:) .* max(Fe2 - F.Fe_top, 0)) .* 1e-3;
Budget.Fe2.I_FeRed = trapz(z, D.R_FeRed(:)) .* 1e-3;
Budget.Fe2.I_FeOx  = trapz(z, D.R_FeOx(:)) .* 1e-3;
Budget.Fe2.I_FeS   = trapz(z, D.R_FeS(:)) .* 1e-3;
Budget.Fe2.I_IrrigSink = I_Fe2_irrig_sink;
Budget.Fe2.J_TopUp = J_Fe2_top_up;
Budget.Fe2.Source = Budget.Fe2.I_FeRed;
Budget.Fe2.Sink = Budget.Fe2.I_FeOx + Budget.Fe2.I_FeS + ...
                  Budget.Fe2.I_IrrigSink + Budget.Fe2.J_TopUp;
Budget.Fe2.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dFe2dt = (S.Fe2(:) - S_prev.Fe2(:)) ./ max(dt, 1e-12);
    Budget.Fe2.Storage = trapz(z, phi .* dFe2dt) .* 1e-3;
end
Budget.Fe2.Residual = Budget.Fe2.Source - Budget.Fe2.Sink - nan0(Budget.Fe2.Storage);
Budget.Fe2.SinkSourceRatio = Budget.Fe2.Sink ./ max(Budget.Fe2.Source, 1e-12);
% =========================
% FeOOH solid budget
% =========================
FeOOH = S.FeOOH(:);
F_FeOOH_top = F.F_FeOx .* 36.5;
if isfield(P, 'Fe_inventory_factor')
    F_FeOOH_top = P.Fe_inventory_factor .* F_FeOOH_top;
end
F_FeOOH_bottom = A(end) .* Grid.v_solid(end) .* FeOOH(end);
% FeOOH reaction terms are equivalent to phi-weighted Fe solute rates.
Budget.FeOOH.TopInput = F_FeOOH_top;
Budget.FeOOH.BottomBurial = F_FeOOH_bottom;
Budget.FeOOH.I_FeRed = trapz(z, D.R_FeRed(:)) .* 1e-3;
Budget.FeOOH.I_FeOx  = trapz(z, D.R_FeOx(:)) .* 1e-3;
Budget.FeOOH.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dFeOOHdt = (S.FeOOH(:) - S_prev.FeOOH(:)) ./ max(dt, 1e-12);
    Budget.FeOOH.Storage = trapz(z, A .* dFeOOHdt);
end
Budget.FeOOH.Net = Budget.FeOOH.TopInput + Budget.FeOOH.I_FeOx ...
                   - Budget.FeOOH.I_FeRed - Budget.FeOOH.BottomBurial ...
                   - nan0(Budget.FeOOH.Storage);
% =========================
% Print
% =========================
fprintf('\n================ RTM Budget ================\n');
fprintf('\n--- CH4 budget, umol/cm2/yr ---\n');
fprintf('Meth source:          %.3f\n', Budget.CH4.I_Meth);
fprintf('AOM sink:             %.3f\n', Budget.CH4.I_AOM);
fprintf('O2 oxidation sink:    %.3f\n', Budget.CH4.I_Ox);
fprintf('Irrigation sink:      %.3f\n', Budget.CH4.I_IrrigSink);
fprintf('Top diffusive efflux: %.3f\n', Budget.CH4.J_TopUp);
fprintf('Storage:              %.3f\n', Budget.CH4.Storage);
fprintf('Sink/source ratio:    %.3f\n', Budget.CH4.SinkSourceRatio);
fprintf('Residual:             %.3f\n', Budget.CH4.Residual);
fprintf('\n--- Fe2 budget, umol/cm2/yr ---\n');
fprintf('Fe reduction source:  %.3f\n', Budget.Fe2.I_FeRed);
fprintf('Fe oxidation sink:    %.3f\n', Budget.Fe2.I_FeOx);
fprintf('FeS sink:             %.3f\n', Budget.Fe2.I_FeS);
fprintf('Irrigation sink:      %.3f\n', Budget.Fe2.I_IrrigSink);
fprintf('Top diffusive efflux: %.3f\n', Budget.Fe2.J_TopUp);
fprintf('Storage:              %.3f\n', Budget.Fe2.Storage);
fprintf('Sink/source ratio:    %.3f\n', Budget.Fe2.SinkSourceRatio);
fprintf('Residual:             %.3f\n', Budget.Fe2.Residual);
fprintf('\n--- FeOOH budget, umol/cm2/yr ---\n');
fprintf('Top FeOOH input:      %.3f\n', Budget.FeOOH.TopInput);
fprintf('FeOOH regeneration:   %.3f\n', Budget.FeOOH.I_FeOx);
fprintf('FeOOH reduction sink: %.3f\n', Budget.FeOOH.I_FeRed);
fprintf('Bottom burial loss:   %.3f\n', Budget.FeOOH.BottomBurial);
fprintf('Storage:              %.3f\n', Budget.FeOOH.Storage);
fprintf('Residual:             %.3f\n', Budget.FeOOH.Net);
fprintf('============================================\n\n');
end
function J_top_down = top_solute_flux(C, C_top, phi, D, v, dz)
% Positive downward. Convert to umol/cm2/yr.
C = C(:);
phi = phi(:);
v = v(:);
if v(1) >= 0
    C_adv = C_top;
else
    C_adv = C(1);
end
F_top = -phi(1) .* D .* (C(1) - C_top) ./ (0.5 .* dz) ...
        + phi(1) .* v(1) .* C_adv;
J_top_down = F_top .* 1e-3;
end
function x = nan0(x)
if isnan(x)
    x = 0;
end
end
```

## File: RTM_Carbonate.m
```matlab
function Carb = RTM_Carbonate(ALK, DIC, T, S)
% RTM_CARBONATE
% Vectorized carbonate speciation solver for the transient PDE RHS.
%
% Uses fixed-iteration vectorized bisection in log10([H+] in uM).
% This avoids fzero and roots inside ode15s RHS calls.
%
% Inputs:
%   ALK, DIC : uM, vectors or scalars
%   T        : deg C
%   S        : salinity
%
% Outputs in Carb:
%   pH, CO3, H2CO3, HCO3, BOH4, OH
ALK = max(real(ALK(:)), 1e-12);
DIC = max(real(DIC(:)), 1e-12);
if nargin < 3 || isempty(T)
    T = 25;
end
if nargin < 4 || isempty(S)
    S = 0.1;
end
T_K = T + 273.15;
% River/low-salinity boron approximation inherited from current code.
B_T = 400 .* (S ./ 35);
lnK1 = 2.83655 - 2307.1266 ./ T_K - 1.5529413 .* log(T_K) ...
     - (0.20760841 + 4.0484 ./ T_K) .* sqrt(S) ...
     + 0.08468345 .* S - 0.00654208 .* S.^1.5 ...
     + log(1 - 0.001005 .* S);
K1 = exp(lnK1) .* 1e6;
lnK2 = -9.226508 - 3351.6106 ./ T_K - 0.2005743 .* log(T_K) ...
     + (-0.106901773 - 23.9722 ./ T_K) .* sqrt(S) ...
     + 0.1130822 .* S - 0.00846934 .* S.^1.5 ...
     + log(1 - 0.001005 .* S);
K2 = exp(lnK2) .* 1e6;
lnKb = (-8966.90 - 2890.53 .* sqrt(S) - 77.942 .* S ...
      + 1.728 .* S.^1.5 - 0.0996 .* S.^2) ./ T_K ...
     + 148.0248 + 137.1942 .* sqrt(S) + 1.62142 .* S ...
     + (-24.4344 - 25.085 .* sqrt(S) - 0.2474 .* S) .* log(T_K) ...
     + 0.053105 .* sqrt(S) .* T_K;
Kb = exp(lnKb) .* 1e6;
Kw = exp(148.96502 - 13847.26 ./ T_K - 23.6521 .* log(T_K) ...
   + (118.67 ./ T_K - 5.977 + 1.0495 .* log(T_K)) .* sqrt(S) ...
   - 0.01615 .* S) .* 1e12;
% Broadcast constants to vector size.
K1 = K1 .* ones(size(ALK));
K2 = K2 .* ones(size(ALK));
Kb = Kb .* ones(size(ALK));
Kw = Kw .* ones(size(ALK));
B_T = B_T .* ones(size(ALK));
% Bracket log10(H_uM).
% Current River_Carbonate used -7 to 2; keep same safe range.
lo = -7 .* ones(size(ALK));
hi =  2 .* ones(size(ALK));
flo = alk_residual(10.^lo, ALK, DIC, B_T, K1, K2, Kb, Kw);
fhi = alk_residual(10.^hi, ALK, DIC, B_T, K1, K2, Kb, Kw);
% If a rare point is not bracketed, keep the nearest endpoint.
bad = flo .* fhi > 0;
% Fixed iteration count is predictable and ode15s-friendly.
for iter = 1:50
    mid = 0.5 .* (lo + hi);
    fmid = alk_residual(10.^mid, ALK, DIC, B_T, K1, K2, Kb, Kw);
    % Residual generally decreases with increasing H.
    go_right = fmid > 0;
    lo(go_right) = mid(go_right);
    hi(~go_right) = mid(~go_right);
end
logH = 0.5 .* (lo + hi);
% Fallback for non-bracketed cells.
if any(bad)
    use_lo = abs(flo) < abs(fhi);
    logH(bad & use_lo) = lo(bad & use_lo);
    logH(bad & ~use_lo) = hi(bad & ~use_lo);
end
H = max(10.^logH, 1e-12);
denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
denom = max(denom, 1e-30);
H2CO3 = DIC ./ denom;
HCO3  = DIC .* (K1 ./ H) ./ denom;
CO3   = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
BOH4  = B_T .* (Kb ./ H) ./ (1 + Kb ./ H);
OH    = Kw ./ H;
Carb = struct();
Carb.pH    = 6 - log10(H);
Carb.H2CO3 = max(real(H2CO3), 1e-12);
Carb.HCO3  = max(real(HCO3),  1e-12);
Carb.CO3   = max(real(CO3),   1e-12);
Carb.BOH4  = max(real(BOH4),  0);
Carb.OH    = max(real(OH),    0);
end
function res = alk_residual(H, ALK, DIC, B_T, K1, K2, Kb, Kw)
    H = max(H, 1e-12);
    denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
    denom = max(denom, 1e-30);
    HCO3 = DIC .* (K1 ./ H) ./ denom;
    CO3  = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
    BOH4 = B_T .* (Kb ./ H) ./ (1 + Kb ./ H);
    OH   = Kw ./ H;
    TA_calc = HCO3 + 2 .* CO3 + BOH4 + OH - H;
    res = TA_calc - ALK;
end
```

## File: RTM_Compare.m
```matlab
function Compare = RTM_Compare(Ode, Pde)
% RTM_COMPARE
% Quantitative ODE-vs-PDE benchmark table.
%
% Usage:
%   Run ODE script first and export Ode struct.
%   Run PDE and get Result.
%   Compare = RTM_Compare(Ode, Result);
P = Pde.State_final;
D = Pde.Diag_final;
G = Pde.Grid;
z_pde = G.z(:);
Compare = table();
Compare = add_row(Compare, 'OPD_cm', ...
    first_depth_below(Ode.z, Ode.O2, 1, Ode.z(end)), ...
    first_depth_below(z_pde, P.O2, 1, z_pde(end)));
Compare = add_row(Compare, 'SO4_depth_10_cm', ...
    first_depth_below(Ode.z, Ode.SO4, 10, Ode.z(end)), ...
    first_depth_below(z_pde, P.SO4, 10, z_pde(end)));
Compare = add_row(Compare, 'O2_bottom_uM', ...
    Ode.O2(end), P.O2(end));
Compare = add_row(Compare, 'SO4_bottom_uM', ...
    Ode.SO4(end), P.SO4(end));
Compare = add_row(Compare, 'Fe2_max_uM', ...
    max(Ode.Fe2), max(P.Fe2));
Compare = add_row(Compare, 'Fe2_bottom_uM', ...
    Ode.Fe2(end), P.Fe2(end));
Compare = add_row(Compare, 'HS_max_uM', ...
    max(Ode.HS), max(P.HS));
Compare = add_row(Compare, 'CH4_max_uM', ...
    max(Ode.CH4), max(P.CH4));
Compare = add_row(Compare, 'CH4_bottom_uM', ...
    Ode.CH4(end), P.CH4(end));
Compare = add_row(Compare, 'FeOOH_top_umol_g', ...
    Ode.FeOOH(1), P.FeOOH(1));
Compare = add_row(Compare, 'FeOOH_max_umol_g', ...
    max(Ode.FeOOH), max(P.FeOOH));
Compare = add_row(Compare, 'FeOOH_bottom_umol_g', ...
    Ode.FeOOH(end), P.FeOOH(end));
Compare = add_row(Compare, 'DIC_bottom_uM', ...
    Ode.DIC(end), P.DIC(end));
Compare = add_row(Compare, 'ALK_bottom_uM', ...
    Ode.ALK(end), P.ALK(end));
Compare = add_row(Compare, 'pH_bottom', ...
    Ode.pH(end), D.pH(end));
Compare = add_row(Compare, 'sigma_top5_mean', ...
    mean(Ode.sigma(Ode.z <= 5)), ...
    mean(D.sigma_carb(z_pde <= 5)));
Compare = add_row(Compare, 'sigma_bottom', ...
    Ode.sigma(end), D.sigma_carb(end));
Compare = add_row(Compare, 'CaCO3_top_pct', ...
    100 .* Ode.CaCO3(1), ...
    100 .* P.CaCO3(1));
Compare = add_row(Compare, 'CaCO3_bottom_pct', ...
    100 .* Ode.CaCO3(end), ...
    100 .* P.CaCO3(end));
% Optional integrated diagnostics if available.
if isfield(Ode, 'RC') && isfield(D, 'RC_uM')
    Compare = add_row(Compare, 'I_RC_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.RC .* 1e9) .* 1e-3, ...
        trapz(z_pde,  D.RC_uM(:)) .* 1e-3);
end
if isfield(Ode, 'R_SRR') && isfield(D, 'R_SRR')
    Compare = add_row(Compare, 'I_SRR_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_SRR) .* 1e-3, ...
        trapz(z_pde,  D.R_SRR(:)) .* 1e-3);
end
if isfield(Ode, 'Rate_Meth') && isfield(D, 'R_Meth')
    Compare = add_row(Compare, 'I_Meth_umolCH4_cm2_yr', ...
        trapz(Ode.z, Ode.Rate_Meth) .* 1e-3, ...
        trapz(z_pde,  D.R_Meth(:)) .* 1e-3);
end
% Redox partition diagnostics.
if isfield(Ode, 'R_respi') && isfield(D, 'R_respi')
    Compare = add_row(Compare, 'I_O2_resp_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.R_respi) .* 1e-3, ...
        trapz(z_pde,  D.R_respi(:)) .* 1e-3);
end
if isfield(Ode, 'R_FeRed') && isfield(D, 'R_FeRed')
    Compare = add_row(Compare, 'I_FeRed_umolFe_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeRed) .* 1e-3, ...
        trapz(z_pde,  D.R_FeRed(:)) .* 1e-3);
    Compare = add_row(Compare, 'I_FeRed_C_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeRed ./ 4) .* 1e-3, ...
        trapz(z_pde,  D.R_FeRed(:) ./ 4) .* 1e-3);
end
if isfield(Ode, 'R_AOM_actual') && isfield(D, 'R_AOM')
    Compare = add_row(Compare, 'I_AOM_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_AOM_actual) .* 1e-3, ...
        trapz(z_pde,  D.R_AOM(:)) .* 1e-3);
elseif isfield(Ode, 'R_AOM') && isfield(D, 'R_AOM')
    Compare = add_row(Compare, 'I_AOM_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_AOM) .* 1e-3, ...
        trapz(z_pde,  D.R_AOM(:)) .* 1e-3);
end
if isfield(Ode, 'R_CH4Ox') && isfield(D, 'R_CH4Ox')
    Compare = add_row(Compare, 'I_CH4Ox_umolCH4_cm2_yr', ...
        trapz(Ode.z, Ode.R_CH4Ox) .* 1e-3, ...
        trapz(z_pde,  D.R_CH4Ox(:)) .* 1e-3);
end
if isfield(Ode, 'R_FeS') && isfield(D, 'R_FeS')
    Compare = add_row(Compare, 'I_FeS_umolFe_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeS) .* 1e-3, ...
        trapz(z_pde,  D.R_FeS(:)) .* 1e-3);
end
if isfield(Ode, 'R_HS_Ox') && isfield(D, 'R_HSOx')
    Compare = add_row(Compare, 'I_HSOx_umolS_cm2_yr', ...
        trapz(Ode.z, Ode.R_HS_Ox) .* 1e-3, ...
        trapz(z_pde,  D.R_HSOx(:)) .* 1e-3);
end
% Derived comparison columns.
Compare.AbsDiff = Compare.PDE - Compare.ODE;
Compare.RelDiff_pct = 100 .* Compare.AbsDiff ./ max(abs(Compare.ODE), 1e-12);
disp(Compare);
end
function T = add_row(T, name, ode_val, pde_val)
newrow = table(string(name), ode_val, pde_val, ...
    'VariableNames', {'Metric','ODE','PDE'});
T = [T; newrow];
end
function depth = first_depth_below(z, x, threshold, default_depth)
z = z(:);
x = x(:);
idx = find(x < threshold, 1, 'first');
if isempty(idx)
    depth = default_depth;
else
    depth = z(idx);
end
end
```

## File: RTM_ExportODE.m
```matlab
function Ode = RTM_ExportODE()
% RTM_EXPORTODE
% Export current ODE workspace variables into a benchmark struct.
%
% Run this from the command window after Run_RTM_1D.m finishes:
%   Ode = RTM_ExportODE();
required_vars = {'z_sed','Oxygen','Sulfate','C_Fe','C_HS','CH4', ...
                 'FeooH','C_DIC','ALK','pH','sigma_carb','CaCO3'};
caller_vars = evalin('base', 'who');
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, caller_vars)
        error('Missing variable "%s" in base workspace. Run Run_RTM_1D.m first.', required_vars{i});
    end
end
Ode = struct();
Ode.z     = evalin('base', 'z_sed(:)');
Ode.O2    = evalin('base', 'Oxygen(:)');
Ode.SO4   = evalin('base', 'Sulfate(:)');
Ode.Fe2   = evalin('base', 'C_Fe(:)');
Ode.HS    = evalin('base', 'C_HS(:)');
Ode.CH4   = evalin('base', 'CH4(:)');
Ode.FeOOH = evalin('base', 'FeooH(:)');
Ode.DIC   = evalin('base', 'C_DIC(:)');
Ode.ALK   = evalin('base', 'ALK(:)');
Ode.pH    = evalin('base', 'pH(:)');
Ode.sigma = evalin('base', 'sigma_carb(:)');
Ode.CaCO3 = evalin('base', 'CaCO3(:)');
optional_vars = {'RC','R_respi','R_SRR','Rate_Meth', ...
                 'R_FeRed','R_FeS','R_FeOx','R_HS_Ox', ...
                 'R_AOM','R_AOM_actual','R_CH4Ox'};
for i = 1:numel(optional_vars)
    name = optional_vars{i};
    if ismember(name, caller_vars)
        Ode.(name) = evalin('base', [name, '(:)']);
    end
end
fprintf('ODE benchmark exported: %d depth nodes.\n', numel(Ode.z));
end
```

## File: RTM_PDE_RHS.m
```matlab
function dYdt = RTM_PDE_RHS(t, Y, Grid, Forcing, Params, Config)
% RTM_PDE_RHS
% Coupled finite-volume method-of-lines RHS for transient 1-D RTM.
%
% All dynamic species are advanced together by ode15s.
% Transport is finite-volume conservative.
% Reactions are computed by RTM_Reaction_Rates.
State = Unpack_State(Y, Grid);
% Safety only for reaction calculation.
% Formal nonnegativity should be handled by ode15s NonNegative and
% depletion-aware reaction limiters.
State = floor_state(State);
[Rates, ~] = RTM_Reaction_Rates(State, Grid, Forcing, Params, Config);
dState = struct();
% ========================================================================
% Solids: conservative mixing + burial + reaction.
% Top boundary is depositional flux.
% ========================================================================
F_OM_lab_top = Forcing.F_lab_OM;               % g/cm2/yr
F_OM_ref_top = Forcing.F_ref_OM;               % g/cm2/yr
F_FeOOH_top  = Forcing.F_FeOx .* 36.5;         % umol/cm2/yr
F_CaCO3_top  = Forcing.F_CaCO3 .* 1e-4;        % g/cm2/yr
if isfield(Params, 'Fe_inventory_factor')
    F_FeOOH_top = Params.Fe_inventory_factor .* F_FeOOH_top;
end
dState.OM_lab = solid_rhs(State.OM_lab, Grid, Rates.OM_lab, F_OM_lab_top);
dState.OM_ref = solid_rhs(State.OM_ref, Grid, Rates.OM_ref, F_OM_ref_top);
dState.FeOOH  = solid_rhs(State.FeOOH,  Grid, Rates.FeOOH,  F_FeOOH_top);
dState.CaCO3  = solid_rhs(State.CaCO3,  Grid, Rates.CaCO3,  F_CaCO3_top);
% ========================================================================
% Solutes: conservative diffusion/advection + exchange + reaction.
% Top boundary is fixed concentration through face flux.
% Bottom boundary is zero diffusive gradient with advective outflow.
% ========================================================================
dState.O2 = solute_rhs(State.O2, Grid, Params.DO2, Rates.O2, Forcing.O2_top);
dState.Fe2 = solute_rhs(State.Fe2, Grid, Params.DH2S, Rates.Fe2, Forcing.Fe_top);
dState.SO4 = solute_rhs(State.SO4, Grid, Params.DSO4, Rates.SO4, Forcing.SO4_top);
dState.HS = solute_rhs(State.HS, Grid, Params.DH2S, Rates.HS, Forcing.HS_top);
dState.CH4 = solute_rhs(State.CH4, Grid, Params.DCH4, Rates.CH4, Forcing.CH4_top);
dState.DIC = solute_rhs(State.DIC, Grid, Params.DHCO3, Rates.DIC, Forcing.DIC_top);
dState.ALK = solute_rhs(State.ALK, Grid, Params.DHCO3, Rates.ALK, Forcing.ALK_top);
dYdt = Pack_State(dState);
end
% ========================================================================
% Local finite-volume helpers.
% ========================================================================
function dCdt = solute_rhs(C, Grid, D, R, C_top)
% d(phi*C)/dt = -div(F) + phi*R + phi*Alpha*(C_top-C)
C = C(:);
R = R(:);
n = Grid.n;
dz = Grid.dz;
phi = Grid.poros(:);
Alpha = Grid.Alpha_exchange(:);
v = Grid.v_fluid(:);
F = zeros(n+1,1);  % face flux, positive downward, unit uM*cm/yr
% Top face: Dirichlet concentration C_top.
F(1) = -phi(1) .* D .* (C(1) - C_top) ./ (0.5 .* dz) ...
       + phi(1) .* v(1) .* upwind_top(C(1), C_top, v(1));
% Internal faces.
for i = 1:n-1
    phi_f = 0.5 .* (phi(i) + phi(i+1));
    v_f   = 0.5 .* (v(i) + v(i+1));
    if v_f >= 0
        C_up = C(i);
    else
        C_up = C(i+1);
    end
    F(i+1) = -phi_f .* D .* (C(i+1) - C(i)) ./ dz ...
             + phi_f .* v_f .* C_up;
end
% Bottom face: zero diffusive gradient, advective outflow if any.
F(n+1) = phi(end) .* v(end) .* C(end);
storage_tendency = -(F(2:end) - F(1:end-1)) ./ dz ...
                   + R ...
                   + phi .* Alpha .* (C_top - C);
% d(phi*C)/dt = -div(F) + R_bulk + phi*Alpha*(C_top-C)
% R is a bulk-volume reaction term, in uM_bulk/yr.
dCdt = storage_tendency ./ max(phi, 1e-12);
end
function dCdt = solid_rhs(C, Grid, R, F_top)
% d(A*C)/dt = -div(F) + A*R
% A = rho*(1-phi), C is solid concentration per gDW.
%
% F_top is imposed depositional flux at sediment-water interface.
% Units must match A*v*C:
%   OM/CaCO3: g/cm2/yr
%   FeOOH:    umol/cm2/yr
C = C(:);
R = R(:);
n = Grid.n;
dz = Grid.dz;
rho = Grid.rho;
phi = Grid.poros(:);
A = rho .* max(1 - phi, 1e-6);
Db = Grid.D_solid_mix(:);
v = Grid.v_solid(:);
F = zeros(n+1,1);  % positive downward
% Top depositional flux.
F(1) = F_top;
% Internal faces.
for i = 1:n-1
    A_f  = 0.5 .* (A(i) + A(i+1));
    Db_f = 0.5 .* (Db(i) + Db(i+1));
    v_f  = 0.5 .* (v(i) + v(i+1));
    if v_f >= 0
        C_up = C(i);
    else
        C_up = C(i+1);
    end
    F(i+1) = -A_f .* Db_f .* (C(i+1) - C(i)) ./ dz ...
             + A_f .* v_f .* C_up;
end
% Bottom face: burial outflow, zero diffusive gradient.
F(n+1) = A(end) .* v(end) .* C(end);
storage_tendency = -(F(2:end) - F(1:end-1)) ./ dz + A .* R;
dCdt = storage_tendency ./ max(A, 1e-12);
end
function C_face = upwind_top(C_cell, C_top, v_face)
if v_face >= 0
    C_face = C_top;
else
    C_face = C_cell;
end
end
function State = floor_state(State)
names = fieldnames(State);
for i = 1:numel(names)
    name = names{i};
    if strcmp(name, 'DIC') || strcmp(name, 'ALK')
        State.(name) = max(real(State.(name)), 1e-12);
    else
        State.(name) = max(real(State.(name)), 0);
    end
end
end
```

## File: RTM_Reaction_Rates.m
```matlab
function [Rates, Diag] = RTM_Reaction_Rates(State, Grid, Forcing, Params, Config)
% RTM_REACTION_RATES
% Coupled reaction-rate calculator for the transient FV-MOL RTM.
%
% This file intentionally extracts mechanism logic from the stable ODE model,
% but does NOT copy bvp4c-specific flux variables or Picard iteration logic.
%
% Required State fields, all column vectors with length Grid.n:
%   Solids:
%     OM_lab      g/gDW
%     OM_ref      g/gDW
%     FeOOH       umol/gDW
%     CaCO3       g/gDW
%
%   Solutes:
%     O2          uM
%     Fe2         uM
%     SO4         uM
%     HS          uM total sulfide
%     CH4         uM
%     DIC         uM
%     ALK         uM
%
% Returned Rates fields have units:
%   Solids: same concentration unit per yr
%   Solutes: uM/yr
%
% Transport, burial, diffusion, bioirrigation, and boundary conditions
% should be handled outside this function in RTM_PDE_RHS.m.
% ---------- vector safety ----------
z   = Grid.z(:);
phi = Grid.poros(:);
ks  = Grid.k_sed(:);
n   = numel(z);
OM_lab = col(State.OM_lab, n, 'OM_lab');
OM_ref = col(State.OM_ref, n, 'OM_ref');
O2     = max(col(State.O2,    n, 'O2'),    0);
FeOOH  = max(col(State.FeOOH, n, 'FeOOH'), 0);
Fe2    = max(col(State.Fe2,   n, 'Fe2'),   0);
SO4    = max(col(State.SO4,   n, 'SO4'),   0);
HS     = max(col(State.HS,    n, 'HS'),    0);
CH4    = max(col(State.CH4,   n, 'CH4'),   0);
DIC    = max(col(State.DIC,   n, 'DIC'),   1e-12);
ALK    = max(col(State.ALK,   n, 'ALK'),   1e-12);
CaCO3  = max(col(State.CaCO3, n, 'CaCO3'), 0);
% ---------- scalar parameters ----------
rho = pick_scalar(Params, Config, Forcing, {'rho'}, 2.73);
k_O2       = pick_scalar(Params, Config, Forcing, {'k_O2'}, 2);
k_SO4      = pick_scalar(Params, Config, Forcing, {'k_SO4'}, 20);
KFEMonod   = pick_scalar(Params, Config, Forcing, {'KFEMonod'}, 1000);
K_HS       = pick_scalar(Params, Config, Forcing, {'K_HS'}, 7);
kFeOx      = pick_scalar(Params, Config, Forcing, {'kFeOx'}, 10);
kFeS       = pick_scalar(Params, Config, Forcing, {'kFeS'}, 1);
Kreox      = pick_scalar(Params, Config, Forcing, {'Kreox'}, 500);
K_CH4_SO4  = pick_scalar(Params, Config, Forcing, {'K_CH4_SO4'}, 100);
K_CH4_O2   = pick_scalar(Params, Config, Forcing, {'K_CH4_O2'}, 1);
k_AOM      = pick_scalar(Params, Config, Forcing, {'k_AOM'}, 0.2);
k_CH4_O2   = pick_scalar(Params, Config, Forcing, {'k_aerobic_CH4'}, 6);
Q10        = pick_scalar(Params, Config, Forcing, {'Q10'}, 2);
T_ref      = pick_scalar(Params, Config, Forcing, {'T_ref'}, 25);
T          = pick_scalar(Forcing, Config, Params, {'T', 'T_future'}, 25);
Salinity   = pick_scalar(Forcing, Config, Params, {'Salinity'}, 0.1);
F_FeOx     = pick_scalar(Forcing, Config, Params, {'F_FeOx'}, 0);
Ca         = pick_scalar(Forcing, Config, Params, {'Ca_top', 'Calcium'}, 1000);
k_ref_factor = pick_scalar(Config, Params, Forcing, {'k_ref_factor'}, 0);
% Carbonate kinetic parameters.
Ksp_ca     = pick_scalar(Params, Config, Forcing, {'Ksp_ca'}, 4.5e5);
k_calcite  = pick_scalar(Params, Config, Forcing, {'k_calcite'}, 1);
k_dis1     = pick_scalar(Params, Config, Forcing, {'k_calcite_dis1'}, 0.005);
k_dis2     = pick_scalar(Params, Config, Forcing, {'k_calcite_dis2'}, 10);
n_form     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO31'}, 1.76);
n_dis1     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO32'}, 0.11);
n_dis2     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO33'}, 4);
gamma_Ca   = pick_scalar(Params, Config, Forcing, {'Calcium_activity'}, 0.6);
gamma_CO3  = pick_scalar(Params, Config, Forcing, {'CO3_activity'}, 0.6);
% Numerical limiter parameters.
SO4_cut = pick_scalar(Config, Params, Forcing, {'SO4_cut'}, 1e-6);
Temp_factor = Q10.^((T - T_ref) ./ 10);
% ========================================================================
% 1. Carbonate speciation first, because HS free fraction needs pH.
% ========================================================================
Carb = RTM_Carbonate(ALK, DIC, T, Salinity);
pH  = Carb.pH;
CO3 = max(Carb.CO3, 1e-12);
sigma_carb = (gamma_Ca .* Ca .* gamma_CO3 .* CO3) ./ Ksp_ca - 1;
% ========================================================================
% 2. OM degradation and residual-carbon cascade.
% ========================================================================
% OM_lab is the redox-active pool.
% Unit follows current ODE convention:
%   RC_mol = mol C / cm3 bulk sediment / yr
%   RC_uM  = umol C / L bulk sediment / yr
RC_mol = Temp_factor .* ks .* OM_lab .* rho .* max(1 - phi, 1e-6) ./ 12;
RC_uM  = RC_mol .* 1e9;
% Solid OM reaction terms.
R_OM_lab = -Temp_factor .* ks .* OM_lab;
R_OM_ref = -Temp_factor .* k_ref_factor .* ks .* OM_ref;
% O2 branch.
R_respi = RC_uM .* O2 ./ max(O2 + k_O2, 1e-12);
R_respi = min(max(R_respi, 0), RC_uM);
RC_after_O2 = max(RC_uM - R_respi, 0);
% Fe branch.
Fe_gate = FeOOH ./ max(FeOOH + KFEMonod, 1e-12);
C_to_Fe_pot = RC_after_O2 .* Fe_gate;
R_FeRed_pot = 4 .* C_to_Fe_pot;  % umol Fe / L / yr
I_FeRed_pot = trapz(z, R_FeRed_pot) .* 1e-3;  % umol Fe / cm2 / yr
I_Fe_supply_ext = F_FeOx .* 36.5;             % mmol/m2/d -> umol/cm2/yr
Fe_supply_scale = min(1, I_Fe_supply_ext ./ max(I_FeRed_pot, 1e-12));
R_FeRed = R_FeRed_pot .* Fe_supply_scale;
C_to_Fe = R_FeRed ./ 4;
RC_after_Fe = max(RC_after_O2 - C_to_Fe, 0);
% SO4 branch with active depletion limiter.
SO4_pos = max(SO4, 0);
f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
f_SRR(SO4_pos <= SO4_cut) = 0;
R_SRR_pot = 0.5 .* RC_after_Fe;
R_SRR = R_SRR_pot .* f_SRR;
% AOM uses current CH4 and current SO4 in the coupled PDE RHS.
f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
f_AOM(SO4_pos <= SO4_cut) = 0;
R_AOM = k_AOM .* CH4 .* f_AOM;
% Methanogenesis branch.
RC_after_SO4 = max(RC_after_Fe - 2 .* R_SRR, 0);
R_Meth = 0.5 .* RC_after_SO4;
% Aerobic methane oxidation.
R_CH4Ox = k_CH4_O2 .* CH4 .* O2 ./ max(O2 + K_CH4_O2, 1e-12);
% ========================================================================
% 3. Fe/S secondary reactions.
% ========================================================================
HS_free = HS ./ (1 + (10.^(6 - pH)) ./ max(K_HS, 1e-12));
HS_free = max(HS_free, 0);
R_FeOx = kFeOx .* Fe2 .* O2;
R_FeS  = kFeS  .* Fe2 .* HS_free;
R_HSOx = Kreox .* HS_free .* O2;
% Convert FeOOH bulk-volume reaction rates to solid concentration rates:
% R_FeRed / R_FeOx: umol/L_bulk/yr
% FeOOH:            umol/gDW
% A = rho*(1-phi):  gDW/cm3_bulk
%
% umol/L_bulk/yr * 1e-3 = umol/cm3_bulk/yr
% solid rate = umol/cm3_bulk/yr / A
bulk_to_solid_Fe = 1e-3 ./ (rho .* max(1 - phi, 1e-6));
R_FeRed_solid = R_FeRed .* bulk_to_solid_Fe;
R_FeOx_solid  = R_FeOx  .* bulk_to_solid_Fe;
% ========================================================================
% 4. Carbonate precipitation/dissolution.
% ========================================================================
solid_to_uM_carb = 1e3 .* 1e6 .* 1e-2 .* rho .* max(1 - phi, 1e-6);
uM_to_solid_carb = 1 ./ solid_to_uM_carb;
R_carb_form_solid = zeros(n,1);
R_carb_dis_solid  = zeros(n,1);
idx_form = sigma_carb > 0;
R_carb_form_solid(idx_form) = ...
    abs(sigma_carb(idx_form)).^n_form .* k_calcite .* uM_to_solid_carb(idx_form);
idx_dis1 = sigma_carb < 0 & sigma_carb > -0.2;
R_carb_dis_solid(idx_dis1) = ...
    abs(sigma_carb(idx_dis1)).^n_dis1 .* k_dis1 .* CaCO3(idx_dis1);
idx_dis2 = sigma_carb <= -0.2;
R_carb_dis_solid(idx_dis2) = ...
    abs(sigma_carb(idx_dis2)).^n_dis2 .* k_dis2 .* CaCO3(idx_dis2);
% Positive means net CaCO3 precipitation into solid phase.
R_carb_net_solid = R_carb_form_solid - R_carb_dis_solid;
R_carb_net_uM    = R_carb_net_solid .* solid_to_uM_carb;
% ========================================================================
% 5. Final explicit reaction source/sink terms for PDE RHS.
% ========================================================================
Rates = struct();
% Solids.
Rates.OM_lab = R_OM_lab;
Rates.OM_ref = R_OM_ref;
Rates.FeOOH  = -R_FeRed_solid + R_FeOx_solid;
Rates.CaCO3  = R_carb_net_solid;
% Solutes.
% For first PDE benchmark, O2 follows the current ODE-compatible primary
% OM respiration sink. Secondary O2 sinks are stored in diagnostics below.
Rates.O2  = -R_respi;
Rates.Fe2 = +R_FeRed - R_FeOx - R_FeS;
% Include sulfide oxidation as sulfate regeneration in the transient sulfur
% ledger. This is more mass-consistent than the current steady FV SO4 solver.
Rates.SO4 = -R_SRR - R_AOM + R_HSOx;
Rates.HS  = +R_SRR + R_AOM - R_FeS - R_HSOx;
Rates.CH4 = +R_Meth - R_AOM - R_CH4Ox;
% DIC/ALK ledger follows explicit reaction signs:
% carbonate precipitation consumes DIC/ALK; dissolution releases them.
Rates.DIC = +R_respi ...
            +R_FeRed ./ 4 ...
            +2 .* R_SRR ...
            +R_Meth ...
            +R_CH4Ox ...
            +R_AOM ...
            -R_carb_net_uM;
Rates.ALK = +0.5 .* R_FeRed ...
            +2 .* R_SRR ...
            +2 .* R_AOM ...
            -2 .* R_FeS ...
            -R_HSOx ...
            -2 .* R_carb_net_uM;
% ========================================================================
% 6. Diagnostics.
% ========================================================================
Diag = struct();
Diag.pH = pH;
Diag.CO3 = CO3;
Diag.H2CO3 = Carb.H2CO3;
Diag.sigma_carb = sigma_carb;
Diag.RC_mol = RC_mol;
Diag.RC_uM = RC_uM;
Diag.R_respi = R_respi;
Diag.R_FeRed = R_FeRed;
Diag.R_SRR = R_SRR;
Diag.R_AOM = R_AOM;
Diag.R_Meth = R_Meth;
Diag.R_CH4Ox = R_CH4Ox;
Diag.R_FeOx = R_FeOx;
Diag.R_FeS = R_FeS;
Diag.R_HSOx = R_HSOx;
Diag.R_carb_form_solid = R_carb_form_solid;
Diag.R_carb_dis_solid = R_carb_dis_solid;
Diag.R_carb_net_solid = R_carb_net_solid;
Diag.R_carb_net_uM = R_carb_net_uM;
Diag.HS_free = HS_free;
Diag.f_SRR = f_SRR;
Diag.f_AOM = f_AOM;
Diag.R_SRR_pot = R_SRR_pot;
Diag.R_FeRed_pot = R_FeRed_pot;
Diag.Fe_supply_scale = Fe_supply_scale;
Diag.I_FeRed_pot = I_FeRed_pot;
Diag.I_Fe_supply_ext = I_Fe_supply_ext;
Diag.O2_secondary_sink = R_FeOx + R_HSOx + R_CH4Ox;
Diag.I_RC = trapz(z, RC_uM) .* 1e-3;
Diag.I_respi = trapz(z, R_respi) .* 1e-3;
Diag.I_FeRed_C = trapz(z, R_FeRed ./ 4) .* 1e-3;
Diag.I_SRR = trapz(z, R_SRR) .* 1e-3;
Diag.I_AOM = trapz(z, R_AOM) .* 1e-3;
Diag.I_Meth = trapz(z, R_Meth) .* 1e-3;
Diag.I_CaCO3_net = trapz(z, R_carb_net_uM) .* 1e-3;
end
% ========================================================================
% Local helpers.
% ========================================================================
function v = col(x, n, name)
    v = x(:);
    if numel(v) ~= n
        error('State.%s must have length %d, but got length %d.', name, n, numel(v));
    end
    v = real(v);
end
function val = pick_scalar(S1, S2, S3, names, default_val)
    val = default_val;
    structs = {S1, S2, S3};
    for is = 1:numel(structs)
        S = structs{is};
        if isempty(S) || ~isstruct(S)
            continue
        end
        for iname = 1:numel(names)
            nm = names{iname};
            if isfield(S, nm) && ~isempty(S.(nm))
                candidate = S.(nm);
                if isscalar(candidate)
                    val = candidate;
                    return
                end
            end
        end
    end
end
```

## File: Run_RTM_1D.m
```matlab
clear all
tic
% ----------------------------- INPUT PARAMETERS ---------------------------
global v_burial Mineral_Mass z_sed Oxygen Sulfate Corg_top
global k_sed k_O2 DSO4 DH2S DO2 DPO4 k_SO4 Kreox  Bioturb Calcium DHCO3 HCO3init
global O2init SO4init HSinit C_organic rho poros RC Alpha_Bioirrig
global R_respi R_SRR R_FeRed RC_after_Fe R_DIC_prod R_ALK_prod R_AOM R_CH4Ox Ksp_ca k_calcite DICinit R1_carb CO3_1 BE P_C_ratio Rviv1 R_FeS  R_FeOx Fe_3_init
global v_burial_Fluid CO3_activity Calcium_activity NPP kFeS FeooH Feinit R_HS_Ox kapatite
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2 CH4init Pinitial DCH4 kFeOx KFEMonod Sulfide Rapat CaCO3 F_CaCO3 O2_root
global C_HS C_Fe n_power_CaCO31 n_power_CaCO32 k_calcite_dis1 n_power_CaCO33 k_calcite_dis2 CaCO3_init Temp_factor
global T_future Rate_Meth Salinity pH K_HS R_AOM_lag R_AOM_actual R_AOM_pot SO4_diag F_FeOx
global F_lab_OM F_ref_OM C_organic_lab C_organic_ref C_organic_total k_ref_factor Fe_inventory_factor
%global KFe_HS Iron_conc R_iron Iron_C P_apaeq R1_carb_disso R1_carb_form
%     if nargin < 1 || isempty(Custom_Config)
%         Config = Config_Baseline();
%     else
%         Config = Custom_Config;
%     end
%     if nargin < 2 || isempty(Custom_Params)
%         Params = Params_Static();
%     else
%         Params = Custom_Params;
%     end
        Params = Params_Static();
    Config = Config_Baseline();
%     Hydro  = Hydro_Preprocessor(Config, Params);
    rho = Params.rho;
    Mineral_Mass = 215;   % keep as legacy until explicitly audited
    k_O2 = Params.k_O2;
    k_SO4 = Params.k_SO4;
    KFEMonod = Params.KFEMonod;
    DSO4 = Params.DSO4;
    DCH4 = Params.DCH4;
    DH2S = Params.DH2S;
    DO2  = Params.DO2;
    DHCO3 = Params.DHCO3;
    DPO4  = Params.DPO4;
    Kreox = Params.Kreox;
    kFeOx = Params.kFeOx;
    kFeS  = Params.kFeS;
    FeC_frac_max = Params.FeC_frac_max;
    Fe_inventory_factor = Params.Fe_inventory_factor;
    K_CH4_SO4 = Params.K_CH4_SO4;
    K_CH4_O2  = Params.K_CH4_O2;
    k_AOM = Params.k_AOM;
    k_aerobic_CH4 = Params.k_aerobic_CH4;
    Ksp_ca = Params.Ksp_ca;
    k_calcite = Params.k_calcite;
    k_calcite_dis1 = Params.k_calcite_dis1;
    k_calcite_dis2 = Params.k_calcite_dis2;
    n_power_CaCO31 = Params.n_power_CaCO31;
    n_power_CaCO32 = Params.n_power_CaCO32;
    n_power_CaCO33 = Params.n_power_CaCO33;
    Calcium_activity = Params.Calcium_activity;
    CO3_activity     = Params.CO3_activity;
    P_C_ratio = Params.P_C_ratio;
    kapatite = Params.kapatite;
%     KFeS = Params.KFeS;
    K_HS = Params.K_HS;
    Q10 = Params.Q10;
    T_ref = Params.T_ref;
    % Site inputs
    n = Config.n;
    Bioirrig_top    = Config.Bioirrig_top;
    Bioirrig_bottom = Config.Bioirrig_bottom;
    Bioirrig_scale  = Config.Bioirrig_scale;
    Lbottom         = Config.Lbottom;
    Bioturbtop      = Config.Bioturbtop;
    Bioturbbottom   = Config.Bioturbbottom;
    bioturbscale    = Config.bioturbscale;
    vbottom         = Config.vbottom;
    vbottom_fluid   = Config.vbottom_fluid;
    porostop        = Config.porostop;
    porosbottom     = Config.porosbottom;
    porosscale      = Config.porosscale;
    O2init          = Config.O2init;
    SO4init         = Config.SO4init;
    DICinit         = Config.DICinit;
    HCO3init        = Config.HCO3init;
    Calcium         = Config.Calcium;
    CH4init         = Config.CH4init;
    Feinit          = Config.Feinit;
    HSinit          = Config.HSinit;
    Pinitial        = Config.Pinitial;
    % Corg_top = Config.Corg_top;
    BE = Config.BE;
    NPP = Config.NPP;
    f_lab = Config.f_lab;
    k_ref_factor = Config.k_ref_factor;
    f_lab = max(0, min(1, f_lab));
    F_OM_total = BE * NPP * 1E-4;   % current transitional OM input
    F_lab_OM   = f_lab * F_OM_total;
    F_ref_OM   = (1 - f_lab) * F_OM_total;
    if Config.use_hydro_npp_multiplier
        NPP = NPP * Hydro.NPP_multiplier;
    end
    F_FeOx = Config.F_FeOx;
    F_CaCO3 = Config.F_CaCO3;
    T_future = Config.T_future;
    Salinity = Config.Salinity;
    ageinit  = Config.ageinit;
    age_root = Config.age_root;
    % Optional weak hydro injection to diffusion coefficients
    if Config.use_hydro_diffusion_multiplier
        DSO4 = DSO4 * Hydro.diffusion_multiplier;
        DCH4 = DCH4 * Hydro.diffusion_multiplier;
        DH2S = DH2S * Hydro.diffusion_multiplier;
        DO2  = DO2  * Hydro.diffusion_multiplier;
        DHCO3 = DHCO3 * Hydro.diffusion_multiplier;
    end
    DOC_root_1      = Config.DOC_root_1;
    O2_root_1       = Config.O2_root_1;
    POC_root_1      = Config.POC_root_1;
% --------------- Calculating initial depth profiles for input parameters ----------------
% k_sed = k_sed / 100;
% k_AOM = k_AOM / 10;
% k_aerobic_CH4 = k_aerobic_CH4 / 5;
% k_calcite = k_calcite / 10;
% kFeOx = kFeOx / 5;
    z_sed = linspace(0, Lbottom, n);
    z_biodiff = linspace(0, Lbottom, 10001);
    dz_sed = Lbottom / (n - 1);
    if Config.use_constant_porosity
        poros = Config.constant_porosity .* ones(1,n);
    elseif Config.use_hydro_phi
        poros = Hydro.phi_mix .* ones(1,n);
    else
        poros = Config.porosbottom + (Config.porostop - Config.porosbottom) .* exp(-z_sed / Config.porosscale);
    end
    Bioturb_1 = Config.Bioturbbottom + (Config.Bioturbtop - Config.Bioturbbottom) .* exp(-z_biodiff / Config.bioturbscale);
    Bioturb = interp1(z_biodiff, Bioturb_1, z_sed);
    Alpha_Bioirrig = Config.Bioirrig_bottom + (Config.Bioirrig_top - Config.Bioirrig_bottom) .* exp(-z_sed / Config.Bioirrig_scale);
    v_burial = Config.vbottom .* (1 - Config.porosbottom) ./ (1 - poros);
    v_burial_Fluid = Config.vbottom_fluid .* (1 + Config.porosbottom) ./ (1 + poros);
    age = Config.ageinit + cumsum(dz_sed ./ v_burial);
    k_sed = 10.^(-0.95 .* log10(age) - 0.81);
    Temp_factor = Q10.^((T_future - T_ref) / 10);
% ----------------------- Initial Carbonate concentration -----------------
[pH_top, CO3_top, ~] = River_Carbonate(HCO3init, DICinit, T_future, Salinity, 1);
% CO3_top = Carb_CO3(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
CO3_1 = CO3_top*ones(1,n);
C_HS  = zeros(1,n); %intial value for sulfide
R_HS_Ox = zeros(1,n);
CaCO3 = 1E5.*zeros(1,n);
Rapat = zeros(1,n);
Rviv1 = zeros(1,n); %umol/Lsed/yr
R_FeS = ones(1,n); %umol/Lsed/yr
pH    = pH_top.*ones(1,n);
R_AOM_lag = zeros(1,n);
R_AOM_pot    = zeros(1,n);   % potential AOM from previous CH4 profile, before SO4 limitation
R_AOM_actual = zeros(1,n);   % actual SO4-supported AOM used by sulfur/carbonate ledger
R_FeOx = zeros(1,n);   % lagged Fe(II) reoxidation used to cap Fe reduction supply
% ----------------------- Correcting organic matter reactivity based on oxygen penetration depth ----------------
% After calculating the oxygen penetration depth, reactivity profiles would
% be corrected using oxic and anoxic power law by Katsev & Crowe (2015).
% ------------- ORGANIC MATTER DEGRADATION --------------------------------
hold on
% Solving ODE
if Bioturbtop == 0
        x = linspace(0,Lbottom,n);
        CorgInit  = (BE * NPP * 1E-4)./(v_burial(1) * rho * (1-poros(1)));
        C_organic = CorgInit*exp(-cumsum(k_sed./v_burial.*dz_sed));
        BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic
else
    nmesh=1000;
    x=linspace(0,Lbottom,nmesh);
    solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
    sol = bvp4c(@organicODE,@organicbc,solinit);
    x = linspace(0,Lbottom,n);
    y = deval(sol,x);
    if min(y) < 0
        fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
    end
    y = max(y, 1e-12);
    C_organic = y(1,:);
    C_organic_lab = C_organic;
    solid_flux = rho .* max(1 - poros, 1e-6) .* max(v_burial, 1e-12);
    C_organic_ref = F_ref_OM ./ solid_flux;
    C_organic_total = C_organic_lab + C_organic_ref;
    C_organic = C_organic_total;  % only for plotting/output compatibility
    BEsed_org = C_organic_total ./ max(C_organic_total(1), 1e-12);
end
% if Bioturbtop == 0
%     x = linspace(0,Lbottom,n);
%     C_organic = Corg_top .* exp(-cumsum((Temp_factor .* k_sed) ./ v_burial .* dz_sed));
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% else
%     nmesh = 1000;
%     x = linspace(0,Lbottom,nmesh);
%     solinit = bvpinit(linspace(0,Lbottom,nmesh), [Corg_top 0]);
%     sol = bvp4c(@organicODE, @organicbc, solinit);
%     x = linspace(0,Lbottom,n);
%     y = deval(sol,x);
%     if min(y) < 0
%         fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
%     end
%     y = max(y, 1e-12);
%     C_organic = y(1,:);
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% end
% ---------------------------- OXYGEN -------------------------------------
RC_lab = Temp_factor .* k_sed .* C_organic_lab .* rho .* ((1-poros)./12);
RC = RC_lab; % molCorg/cm3/yr mineralization rate
O2_root = zeros(1,n);
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@O2_ODE,@O2_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('Initial O2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_O2 = y(1,:);
Oxygen = C_O2;
% ---------------------------- OXYGEN PENTRATION DEPTS_ --------------------
count_OPD = 0;
for i=1:n
    if Oxygen(1,i) < 1
        count_OPD = count_OPD + 1;
        OPD_1(1,count_OPD) = z_sed(1,i);
        num_OPD1(1,count_OPD) = i;
    end
end
if min(Oxygen) > 1
    OPD_1 = Lbottom;
    num_OPD1 = n;
end
OPD = min(OPD_1);
num_OPD = min(num_OPD1);
mm_count = 0;
for i=1:n
    FeOxInit_raw = (F_FeOx .* 36.5) ./ ...
    (v_burial(1) .* rho .* max(1 - poros(1), 1e-6));
    FeOxInit = Fe_inventory_factor .* FeOxInit_raw;
    Ironoxy(1,i) = 0.01.*FeOxInit.*((1/(1+exp(z_sed(1,i)-OPD)))+2*exp(-((z_sed(1,i)-OPD)^2)/2));
    mm_count=mm_count+1;
    FeOx(1,mm_count)=Ironoxy(1,i);
end
% FeooH =  FeOx;wor
% FeooH =  zeros(1,n);
FeooH = max(FeOx, 1e-12);
% ---------------------------- CORRECTING ACTIVITY PROFILES ---------------
b_oxic = 0.95;%0.977;
a_oxic = 0.81;%0.312;
b_anoxic = 0.95;%0.857;
a_anoxic = 0.81;%1.1;
% % slower reactivity decay below OPD
% b_anoxic = 0.75;
% a_anoxic = 0.60;
% k_sed(1,1:num_OPD) = 10.^(-b_oxic*log10(age(1,1:num_OPD)) - a_oxic);  % Oxic
% k_sed(1,(num_OPD+1):n) = 10.^(-b_anoxic*log10(age(1,(num_OPD+1):n)) - a_anoxic);  % Anoxic
k_oxic   = 10.^(-b_oxic  .* log10(age) - a_oxic);
k_anoxic = 10.^(-b_anoxic .* log10(age) - a_anoxic);
transition_width = 0.5;   % cm
w_oxic = 1 ./ (1 + exp((z_sed - OPD) ./ transition_width));
k_sed = w_oxic .* k_oxic + (1 - w_oxic) .* k_anoxic;
% ------------------- Blue Carbon ------------------------------------
% DOC release
depth_rootzone = 10; % seagrass root length (cm)
z_root = 0:0.1:depth_rootzone;
mu_root = 4; % value for the center of rootzone in normal distribution
sigma_root = 6; % sigma for normal distribution of flux in the rootzone
DOC_root_2 = DOC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/year
DOC_root_22 = interp1(z_root,DOC_root_2,z_sed);
% poros_root = interp1(z_sed,poros,z_sed);
for i_root = 1:n
    if depth_rootzone >= z_sed(1,i_root)
DOC_root(1,i_root) = 1E6 * (DOC_root_22(1,i_root)./(z_sed(1,2)-z_sed(1,1))); % (umol/l/year)
    else
DOC_root(1,i_root) = 0;
    end
end
% O2 release
muO2_root = 4; % value for the center of rootzone in normal distribution
sigmaO2_root = 6; % sigma for normal distribution of flux in the rootzone
O2_root_2 = O2_root_1 * 1E-4 * normpdf(z_root,muO2_root,sigmaO2_root);  % mmol/cm2/year
O2_root_22 = interp1(z_root,O2_root_2,z_sed);
for i_root = 1:n
    if depth_rootzone >= z_sed(1,i_root)
O2_root(1,i_root) = 1E6 * (O2_root_22(1,i_root)./(z_sed(1,2)-z_sed(1,1))); % (umol/l/year)
    else
O2_root(1,i_root) = 0;
    end
end
% Seagrass POC release
muPOC_root = 4; % value for the center of rootzone in normal distribution
sigmaPOC_root = 6; % sigma for normal distribution of flux in the rootzone
POC_root_2 = POC_root_1 * normpdf(z_root,muPOC_root,sigmaPOC_root);  % mmol/cm2/year
POC_root = interp1(z_root,POC_root_2,z_sed);
k_sed_root = 10.^(-0.95*log10(age_root) - 0.8); % more reactive
poros2 = interp1(z_sed,poros,z_sed);
for i_root = 1:n
    if depth_rootzone >= z_sed(1,i_root)
% RC_root(1,i_root) = Temp_factor.*k_sed_root .*POC_root(1,i_root).*rho.*((1-poros2(1,i_root))./(12)); % (umol/l/year)
RC_root(1,i_root) = 1.*k_sed_root .*POC_root(1,i_root).*rho.*((1-poros2(1,i_root))./(12)); % (umol/l/year)
    else
RC_root(1,i_root) = 0;
    end
end
% -------------------------------------------------------------------------
iteration = 1;
count_loop = 1;
max_outer_iter = 12;
min_outer_iter = 5;
conv_tol = 0.05;     % 5% profile/rate convergence
C_Fe_prev   = [];
FeooH_prev  = [];
Sulfate_prev = [];
CH4_prev    = [];
ALK_prev    = [];
pH_prev     = [];
C_HS_prev = [];
R_FeS_prev = [];
while iteration <= max_outer_iter
% ------------- ORGANIC MATTER DEGRADATION --------------------------------
hold on
% Solving ODE
if Bioturbtop == 0
x = linspace(0,Lbottom,n);
CorgInit  = (BE * NPP * 1E-4)./(v_burial(1) * rho * (1-poros(1)));
C_organic = CorgInit.*exp(-cumsum((Temp_factor.*k_sed)./v_burial.*dz_sed));
BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic
else
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@organicODE,@organicbc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_organic = y(1,:);
C_organic_lab = C_organic;
solid_flux = rho .* max(1 - poros, 1e-6) .* max(v_burial, 1e-12);
C_organic_ref = F_ref_OM ./ solid_flux;
C_organic_total = C_organic_lab + C_organic_ref;
C_organic = C_organic_total;  % only for plotting/output compatibility
BEsed_org = C_organic_total ./ max(C_organic_total(1), 1e-12);
end
% if Bioturbtop == 0
%     x = linspace(0,Lbottom,n);
%     C_organic = Corg_top .* exp(-cumsum((Temp_factor .* k_sed) ./ v_burial .* dz_sed));
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% else
%     nmesh = 1000;
%     x = linspace(0,Lbottom,nmesh);
%     solinit = bvpinit(linspace(0,Lbottom,nmesh), [Corg_top 0]);
%     sol = bvp4c(@organicODE, @organicbc, solinit);
%     x = linspace(0,Lbottom,n);
%     y = deval(sol,x);
%     if min(y) < 0
%         fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
%     end
%     y = max(y, 1e-12);
%     C_organic = y(1,:);
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% end
% ---------------------------- OXYGEN -------------------------------------
RC_lab = Temp_factor .* k_sed .* C_organic_lab .* rho .* ((1-poros)./12);
RC = RC_lab + RC_root; % molCorg/cm3/yr mineralization rate
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@O2_ODE,@O2_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('O2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_O2 = y(1,:);
Oxygen = C_O2;
% -------- update OPD-based k_sed every iteration --------
OPD_1 = [];
num_OPD1 = [];
for i = 1:n
    if Oxygen(i) < 1
        OPD_1(end+1) = z_sed(i);
        num_OPD1(end+1) = i;
    end
end
if isempty(OPD_1)
    OPD = Lbottom;
    num_OPD = n;
else
    OPD = OPD_1(1);
    num_OPD = num_OPD1(1);
end
R_respi = RC .* (Oxygen ./ (Oxygen + k_O2)) .* 1E9;   % umol/L/yr
% residual-carbon cascade: O2 first, then Fe, then SO4, then CH4
RC_total_uM = RC .* 1E9;                              % umol C / L / yr
RC_after_O2 = max(RC_total_uM - R_respi, 0);
% FeooH_old = FeooH;
% % cap Fe reduction so Fe branch does not consume nearly all residual carbon
% Fe_gate = FeooH_old ./ (FeooH_old + KFEMonod);
% C_to_Fe = min(FeC_frac_max .* RC_after_O2, RC_after_O2 .* Fe_gate);
% R_FeRed = 4 .* C_to_Fe;
FeooH_old = FeooH;
% Potential Fe reduction demand from carbon and local FeOOH availability.
% FeC_frac_max is no longer active.
Fe_gate = FeooH_old ./ max(FeooH_old + KFEMonod, 1e-12);
C_to_Fe_pot = RC_after_O2 .* Fe_gate;
R_FeRed_pot = 4 .* C_to_Fe_pot;
I_FeRed_pot = trapz(z_sed, R_FeRed_pot) .* 1e-3;
solid1 = rho .* max(1 - poros(1), 1e-6);
J_Fe3_burial_in = solid1 .* v_burial(1) .* FeooH_old(1);
J_Fe3_mix_in = -solid1 .* Bioturb(1) .* ...
    ((FeooH_old(2) - FeooH_old(1)) ./ dz_sed);
I_Fe_transport_supply = max(0, J_Fe3_burial_in + J_Fe3_mix_in);
I_Fe_recycle_lag = trapz(z_sed, max(R_FeOx,0)) .* 1e-3;
I_Fe_supply_ext = F_FeOx .* 36.5;  % mmol/m2/d -> umol/cm2/yr
Fe_recycle_cap_frac = 0.25;
I_Fe_recycle_eff = min(I_Fe_recycle_lag, ...
                       Fe_recycle_cap_frac * I_Fe_supply_ext);
% I_Fe_supply_cap = I_Fe_supply_ext + I_Fe_recycle_eff;
I_Fe_supply_cap = I_Fe_supply_ext ;%...
                % + I_Fe_transport_supply ...
                % + I_Fe_recycle_lag;
Fe_supply_scale = min(1, I_Fe_supply_cap ./ max(I_FeRed_pot, 1e-12));
R_FeRed = R_FeRed_pot .* Fe_supply_scale;
C_to_Fe = R_FeRed ./ 4;
C_HS_used_for_Fe2 = C_HS;
pH_used_for_Fe2 = pH;
% ------------------------ IRON(II) ---------------------------------------
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@Fe_ODE,@Fe_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('Fe2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_Fe = y(1,:);
% Iron_C(iteration,:) = C_Fe;
%
% for i=1:n
%   if Iron_C(iteration,i) < 0
%       Iron_C(iteration,i) = 0;
%   end
% end
% F_diff_Fe(1,count_loop) = DH2S.*((C_Fe(1,2) - C_Fe(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr
% ------------------------ IRON(III) ---------------------------------------
% Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
% R_iron(count_loop,:) = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9; %rate of iron reduction umol/l/year
% R_FeOx = (kFeOx.*C_Fe.*Oxygen);
% % R_FeOx_1(count_loop,:) = R_FeOx;
% % Fe_3_init  = 365.*1E2.*(F_FeOx)./(v_burial(1));  %umol/l
% Fe_3_init = 36.5.*(F_FeOx).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %umol/l
%
% % KFEMonod = 2000;
%
% % % Solving ODE
%
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[Fe_3_init 0]);
% sol = bvp4c(@Fe3_ODE,@Fe3_bc,solinit);
%
% x = linspace(0,Lbottom,n);
%
% y = deval(sol,x);
%
% if min(y) < 0
%     fprintf('Fe3：%.2e\n', min(y));
% end
% y = max(y, 1e-12);
%
% C_Fe_3 = y(1,:);
% FeooH = C_Fe_3;
% % FeooH = max(C_Fe_3, 0.2 .* FeOx);
% RC_after_Fe = max(RC_after_O2 - R_FeRed ./ 4, 0);   % umol C / L / yr
R_FeOx = kFeOx .* C_Fe .* Oxygen;
Fe_3_init_raw = 36.5 .* F_FeOx .* (poros(1) ./ max(1 - poros(1), 1e-6)) ./ (v_burial(1) .* rho);  % umol/g
Fe_3_init = Fe_inventory_factor .* Fe_3_init_raw;
% Fe3 top boundary is a fixed FeOOH inventory derived from Fe input using a global Fe_inventory_factor.
% Do not convert F_FeOx into a fixed FeOOH surface concentration here.
% Solve Fe(III) with the current Fe3_ODE instead of explicit forward update
nmesh = 1000;
x = linspace(0, Lbottom, nmesh);
opts_Fe3 = bvpset('NMax', 5000, 'RelTol', 1e-4);
guess_fun = @(xq) [max(interp1(z_sed, FeooH, xq, 'linear', 'extrap'), 1e-12); 0];
solinit = bvpinit(linspace(0, Lbottom, nmesh), guess_fun);
sol = bvp4c(@Fe3_ODE, @Fe3_bc, solinit, opts_Fe3);
x = linspace(0, Lbottom, n);
y = deval(sol, x);
if min(y(1,:)) < 0
    fprintf('Fe3 is negative in the current iteration! Minimum：%.2e\n', min(y(1,:)));
end
C_Fe_3_new = max(y(1,:), 1e-12);
% Under-relaxed FeOOH update for the NEXT outer iteration.
% Do not recompute R_FeRed inside the same iteration.
relax_Fe3 = 0.5;
FeooH = (1 - relax_Fe3) .* FeooH_old + relax_Fe3 .* C_Fe_3_new;
% Use the same C_to_Fe / R_FeRed that was used to solve Fe2 in this iteration.
RC_after_Fe = max(RC_after_O2 - C_to_Fe, 0);
% % ------------------------ SO4 budget pre-diagnostic ----------------------
% R_SRR_pot_max = 0.5 .* RC_after_Fe;  % maximum sulfate demand, umol SO4/L/yr
%
% I_SO4_demand_pot = trapz(z_sed, R_SRR_pot_max) .* 1e-3;  % umol SO4/cm2/yr
%
% % crude lower-bound/top-down supply estimates
% J_SO4_diff_ref = DSO4 .* SO4init ./ Lbottom .* 1e-3;      % umol SO4/cm2/yr
% I_SO4_irrig_max = trapz(z_sed, Alpha_Bioirrig .* SO4init) .* 1e-3;
%
% SO4_supply_ref = J_SO4_diff_ref + I_SO4_irrig_max;
%
% fprintf('\n--- SO4 pre-diagnostic ---\n');
% fprintf('Potential SO4 demand from RC_after_Fe: %.3f umol/cm2/yr\n', I_SO4_demand_pot);
% fprintf('Reference diffusive SO4 supply:        %.3f umol/cm2/yr\n', J_SO4_diff_ref);
% fprintf('Max irrigation SO4 supply:             %.3f umol/cm2/yr\n', I_SO4_irrig_max);
% fprintf('Supply reference total:                %.3f umol/cm2/yr\n', SO4_supply_ref);
% fprintf('Demand / supply_ref:                   %.2f\n', I_SO4_demand_pot ./ max(SO4_supply_ref, 1e-12));
% fprintf('--------------------------\n\n');
% ------------------------ SULFATE ---------------------------------------
% Use positivity-preserving finite-volume depletion-front solver.
[Sulfate, R_SRR, R_AOM_actual, SO4_diag] = Solve_SO4_FV_Front();
% fprintf('SO4 FV front: min=%.3e uM, front=%.2f cm, I_SRR=%.3f, I_AOM=%.3f, I_pot=%.3f umol/cm2/yr\n', ...
%     SO4_diag.min_SO4, SO4_diag.front_depth, SO4_diag.I_SRR, SO4_diag.I_AOM, SO4_diag.I_demand_pot);
fprintf('SO4 I_demand_pot = %.3f\n', SO4_diag.I_demand_pot);
fprintf('SO4 I_SRR        = %.3f\n', SO4_diag.I_SRR);
fprintf('SO4 I_AOM        = %.3f\n', SO4_diag.I_AOM);
fprintf('SO4 I_irrig      = %.3f\n', SO4_diag.I_irrig_source);
fprintf('SO4 J_top_down   = %.3f\n', SO4_diag.J_top_down);
fprintf('SO4 supply/demand = %.3f\n', ...
    (SO4_diag.I_irrig_source + SO4_diag.J_top_down) ./ ...
    max(SO4_diag.I_SRR + SO4_diag.I_AOM, 1e-12));
% % ------------------------ SULFATE ---------------------------------------
%
% % Solving ODE
%
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% % solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[SO4init 0]);
% sol = bvp4c(@SO4_ODE,@SO4_bc,solinit);
%
% x = linspace(0,Lbottom,n);
%
% y = deval(sol,x);
%
% if min(y) < 0
%     fprintf('SO4 is negative in the current iteration! Minimum：%.2e\n', min(y));
% end
% y = max(y, 1e-12);
%
% % C_SO4 = y(1,:);
% % Sulfate = C_SO4;
% %
% % %R_SRR = RC.*Inhib.* (Sulfate./(Sulfate+k_SO4)).*1E9; %rate of sulfate reduction umol/l/year
% %
% % R_SRR = 0.5 .* RC_after_Fe .* (Sulfate ./ (Sulfate + k_SO4));   % umol SO4/L/yr
%
% SO4_raw = real(y(1,:));
%
% if min(SO4_raw) < -1
%     fprintf('SO4 is substantially negative before clipping! Minimum：%.2e\n', min(SO4_raw));
% elseif min(SO4_raw) < 0
%     fprintf('SO4 has small negative numerical values before clipping. Minimum：%.2e\n', min(SO4_raw));
% end
%
% Sulfate = max(SO4_raw, 0);
%
% SO4_cut = 1e-6;   % uM; below this sulfate is treated as depleted
%
% SO4_pos = max(Sulfate, 0);
%
% f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
% f_SRR(SO4_pos <= SO4_cut) = 0;
%
% R_SRR = 0.5 .* RC_after_Fe .* f_SRR;   % actual organoclastic SRR, umol SO4/L/yr
%
% % Actual AOM supported by current sulfate, using previous CH4-derived potential.
% f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
% f_AOM(SO4_pos <= SO4_cut) = 0;
%
% R_AOM_actual = R_AOM_pot .* f_AOM;     % actual AOM used by sulfur/carbonate ledger
% ------------------------ SULFIDE ---------------------------------------
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@HS_ODE,@HS_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('HS is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_HS = max(real(y(1,:)), 0);
HS_conc = C_HS ./ (1 + ((10.^(6-pH))./K_HS));
Sulfide(iteration,:) = C_HS;
% sigma_FeS_1 =  (Iron_C(iteration,:).*HS_conc)./((10.^(6-pH)).*KFeS);
% delta_FeS  = (sigma_FeS_1 - 1);
%
% for i=1:n
%   if delta_FeS (1,i) > 0
%       delta_FeS1(1,i) = 1;
%   else
%       delta_FeS1(1,i) = 0;
%   end
% end
% R_FeS = kFeS.*C_Fe.*C_HS;
R_FeS = kFeS .* C_Fe .* HS_conc;
% R_FeS_1(count_loop,:) = R_FeS;
% R_FeS_store(iteration,:) = R_FeS;
for i=1:n
  if R_FeS(1,i) < 0
      R_FeS(1,i) = 0;
  end
end
R_HS_Ox = max(Kreox .* HS_conc .* Oxygen, 0);   % umol/L/yr
% ---------------- Fe diagnostics ----------------
% I_FeRed = trapz(z_sed, R_FeRed) .* 1e-3;
% I_FeOx  = trapz(z_sed, R_FeOx)  .* 1e-3;
% I_FeS   = trapz(z_sed, R_FeS)   .* 1e-3;
%
% I_Fe_irrig = trapz(z_sed, Alpha_Bioirrig .* max(C_Fe - Feinit, 0)) .* 1e-3;
%
% J_Fe2_top_up = DH2S .* ((C_Fe(2) - C_Fe(1)) ./ dz_sed) .* 1e-3;
%
% Fe_sink_total = I_FeOx + I_FeS + I_Fe_irrig + J_Fe2_top_up;
% fprintf('\n--- Fe diagnostics ---\n');
% fprintf('I_FeRed_pot          = %.3f umol Fe/cm2/yr\n', I_FeRed_pot);
% fprintf('I_Fe external supp = %.3f umol Fe/cm2/yr\n', I_Fe_supply_ext);
% fprintf('I_Fe transport supply= %.3f umol Fe/cm2/yr\n', I_Fe_transport_supply);
% fprintf('I_Fe recycle lag     = %.3f umol Fe/cm2/yr\n', I_Fe_recycle_lag);
% fprintf('I_Fe supply cap      = %.3f umol Fe/cm2/yr\n', I_Fe_supply_cap);
% fprintf('Fe supply scale      = %.3f\n', Fe_supply_scale);
% fprintf('I_FeRed          = %.3f umol Fe/cm2/yr\n', I_FeRed);
% fprintf('I_FeOx           = %.3f umol Fe/cm2/yr\n', I_FeOx);
% fprintf('I_FeS            = %.3f umol Fe/cm2/yr\n', I_FeS);
% fprintf('I_Fe irrigation  = %.3f umol Fe/cm2/yr\n', I_Fe_irrig);
% fprintf('J_Fe2 top up     = %.3f umol Fe/cm2/yr\n', J_Fe2_top_up);
% fprintf('Fe sink total    = %.3f umol Fe/cm2/yr\n', Fe_sink_total);
% fprintf('max Fe2          = %.1f uM\n', max(C_Fe));
% fprintf('max Fe3          = %.1f umol/g\n', max(FeooH));
% fprintf('Fe sink/source   = %.3f\n', Fe_sink_total / max(I_FeRed,1e-12));
%
% Fe_budget_residual = Fe_sink_total - I_FeRed;
%
% fprintf('Fe residual       = %.3f umol Fe/cm2/yr\n', Fe_budget_residual);
% fprintf('Fe residual/source= %.3f\n', Fe_budget_residual ./ max(I_FeRed,1e-12));
% Fe_sink_no_top = I_FeOx + I_FeS + I_Fe_irrig;
%
%
% fprintf('Fe sink no top    = %.3f umol Fe/cm2/yr\n', Fe_sink_no_top);
% fprintf('Fe no-top/source  = %.3f\n', Fe_sink_no_top ./ max(I_FeRed,1e-12));
%
% HS_used_for_Fe2 = C_HS_used_for_Fe2 ./ (1 + ((10.^(6 - pH_used_for_Fe2)) ./ K_HS));
% R_FeS_used_for_Fe2 = kFeS .* C_Fe .* HS_used_for_Fe2;
% I_FeS_used_for_Fe2 = trapz(z_sed, R_FeS_used_for_Fe2) .* 1e-3;
%
% fprintf('I_FeS used Fe2   = %.3f umol Fe/cm2/yr\n', I_FeS_used_for_Fe2);
% fprintf('Fe used sink/source = %.3f\n', ...
%     (I_FeOx + I_FeS_used_for_Fe2 + I_Fe_irrig + J_Fe2_top_up) ./ max(I_FeRed,1e-12));
fprintf('SO4 min = %.2f uM\n', min(Sulfate));
fprintf('SO4 bottom = %.2f uM\n', Sulfate(end));
idx_SO4_10 = find(Sulfate < 10, 1);
if isempty(idx_SO4_10)
    fprintf('SO4 depletion front <10 uM: none\n');
else
    fprintf('SO4 depletion front <10 uM: %.2f cm\n', z_sed(idx_SO4_10));
end
I_FeRed = trapz(z_sed, R_FeRed) .* 1e-3;
I_FeS   = trapz(z_sed, R_FeS) .* 1e-3;
I_FeOx_int = trapz(z_sed, R_FeOx) .* 1e-3;
fprintf('I_FeRed = %.3f umol Fe/cm2/yr\n', I_FeRed);
fprintf('I_FeS   = %.3f umol Fe/cm2/yr\n', I_FeS);
fprintf('I_FeOx  = %.3f umol Fe/cm2/yr\n', I_FeOx_int);
fprintf('FeS/FeRed = %.3f\n', I_FeS ./ max(I_FeRed,1e-12));
fprintf('max Fe2 = %.2f uM\n', max(C_Fe));
fprintf('max FeOOH = %.2f umol/g\n', max(FeooH));
fprintf('----------------------\n\n');
% ------------------------ METHANE ---------------------------------------
RC_after_SO4 = max(RC_after_Fe - 2 .* R_SRR, 0);   % umol C/L/yr
Rate_Meth = 0.5 .* RC_after_SO4;                   % umol CH4/L/yr
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@CH4_ODE,@CH4_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('CH4 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_CH4 = y(1,:);
CH4 = C_CH4;
SO4_pos = max(Sulfate, 0);
f_AOM_post = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
f_AOM_post(SO4_pos <= 1e-6) = 0;
R_AOM   = k_AOM .* CH4 .* f_AOM_post;  % diagnostic current AOM based on current CH4/SO4
R_CH4Ox = k_aerobic_CH4 .* CH4 .* (Oxygen ./ max(Oxygen + K_CH4_O2, 1e-12));
% Use the same actual AOM that SO4/HS used this iteration.
% This keeps sulfur and carbonate ledger internally consistent.
R_DIC_prod = R_respi ...
           + R_FeRed ./ 4 ...
           + 2 .* R_SRR ...
           + Rate_Meth ...
           + R_CH4Ox ...
           + R_AOM_actual;
R_ALK_prod = 0.5 .* R_FeRed ...
           + 2 .* R_SRR ...
           + 2 .* R_AOM_actual ...
           - 2 .* R_FeS ...
           - R_HS_Ox;
% Store potential AOM for the next iteration's SO4 solve.
R_AOM_pot = k_AOM .* CH4;
R_AOM_lag = R_AOM;
% ---------------- Redox partition diagnostics ----------------
% I_RC     = trapz(z_sed, RC .* 1E9) .* 1e-3;        % umol C/cm2/yr
% I_respi  = trapz(z_sed, R_respi) .* 1e-3;          % umol C/cm2/yr
% I_FeC    = trapz(z_sed, R_FeRed ./ 4) .* 1e-3;     % umol C/cm2/yr
% I_SO4C   = trapz(z_sed, 2 .* R_SRR) .* 1e-3;       % umol C/cm2/yr
% I_methC  = trapz(z_sed, 2 .* Rate_Meth) .* 1e-3;   % umol C/cm2/yr
%
% fprintf('\n--- Carbon redox partition ---\n');
% fprintf('I_RC total       = %.3f umol C/cm2/yr\n', I_RC);
% fprintf('I_O2 respiration = %.3f (%.1f%%)\n', I_respi, 100*I_respi/max(I_RC,1e-12));
% fprintf('I_Fe reduction   = %.3f (%.1f%%)\n', I_FeC,   100*I_FeC/max(I_RC,1e-12));
% fprintf('I_SO4 reduction  = %.3f (%.1f%%)\n', I_SO4C,  100*I_SO4C/max(I_RC,1e-12));
% fprintf('I_methanogenesis = %.3f (%.1f%%)\n', I_methC, 100*I_methC/max(I_RC,1e-12));
% fprintf('Closure ratio     = %.3f\n', (I_respi + I_FeC + I_SO4C + I_methC)/max(I_RC,1e-12));
% fprintf('------------------------------\n\n');
% ------------------------------- Coupled Carbonate -----------------------------------
CaCO3_init = 1E-4 .* (F_CaCO3) .* (poros(1) / (1 - poros(1))) / (v_burial(1) * rho);  % gr/grDw
solinit_coupled = bvpinit(linspace(0, Lbottom, nmesh), ...
    [DICinit, 0, HCO3init, 0, CaCO3_init, 0]);
sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled);
y_coupled = deval(sol_coupled, x);
C_DIC  = max(real(y_coupled(1,:)), 1e-12);
C_alka = max(real(y_coupled(3,:)), 1e-12);
CaCO3  = max(real(y_coupled(5,:)), 0);
ALK = C_alka;
F_diff_DIC = DHCO3 .* ((C_DIC(1,2)  - C_DIC(1,1))  ./ (x(1,2) - x(1,1))) * 1E-3;
F_diff     = DHCO3 .* ((C_alka(1,2) - C_alka(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;
% -------------- pH / CO3 / H2CO3 postprocessing -----------------------
pH_1    = zeros(1,n);
CO3_1   = zeros(1,n);
C_H2CO3 = zeros(1,n);
for i = 1:n
    [pH_1(i), CO3_1(i), C_H2CO3(i)] = River_Carbonate(ALK(i), C_DIC(i), T_future, Salinity, 1);
end
pH = pH_1;
sigma_carb = (Calcium_activity .* Calcium .* CO3_activity .* CO3_1) ./ Ksp_ca - 1;
unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho .* max(1 - poros, 1e-6));
R_carb_form = zeros(1,n);
R_carb_disso = zeros(1,n);
for i = 1:n
    if sigma_carb(i) > 0
        R_carb_form(i) = abs(sigma_carb(i))^n_power_CaCO31 .* k_calcite .* unit_conversion(i);
    else
        R_carb_form(i) = 0;
    end
    if (sigma_carb(i) < 0) && (sigma_carb(i) > -0.2)
        R_carb_disso(i) = abs(sigma_carb(i))^n_power_CaCO32 .* k_calcite_dis1 .* CaCO3(i);
    elseif sigma_carb(i) <= -0.2
        R_carb_disso(i) = abs(sigma_carb(i))^n_power_CaCO33 .* k_calcite_dis2 .* CaCO3(i);
    else
        R_carb_disso(i) = 0;
    end
end
R1_carb = R_carb_form - R_carb_disso .* (1E3 * 1E6 * 1E-2 * rho .* max(1 - poros, 1e-6));
% ------------------------ PHOSPHOROUS ------------------------------------
% % Solving ODE
%
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
% sol = bvp4c(@Phos_ODE,@Phos_bc,solinit);
%
% x = linspace(0,Lbottom,n);
%
% y = deval(sol,x);
% C_Phos = y(1,:);
% PO4 = C_Phos;
%
% P_apaeq = (305.831-74.654.*pH + 4.583.*pH.^2).*(1-(T_apa-10).*X_apa);
% FeOxide=interp1(z_sed,FeooH./1000,x);
% Kadsp = (Kp_ads.*FeOxide.*S_ads_coef.*((1-poros)./poros).*rho.*1E3)./((10.^-pH.*1E6)+Kp_ads.*C_Phos);
% ads_P = Kadsp.*C_Phos; %Fe-bound P
%
% PO4_diss = C_Phos - ads_P; %PO4 dissolved
% delta_apat2 = (PO4_diss-P_apaeq);
%
% for i=1:n
%   if delta_apat2 (1,i) > 0
%       delta_apat1(1,i) = 1;
%   else
%       delta_apat1(1,i) = 0;
%   end
% end
%
% sigma1_viv = ((Iron_C.^3).*(PO4_diss.^2))./Kviv; %umol/Lsed/yr
% delta_viv1 = (sigma1_viv-1);
%
% for i=1:n
%   if delta_viv1 (1,i) > 0
%       delta_viv1(1,i) = 1;
%   else
%       delta_viv1(1,i) = 0;
%   end
% end
%
% Rapat = kapatite.*delta_apat1.*(delta_apat2.^n_apa);
% Rviv1 = kviv.*(sigma1_viv.^alpha_viv-1).*delta_viv1.*(1-poros).*rho.*1E3; %umol/Lsed/yr
% -------------------------- Convergence coefficient ----------------------
if iteration > 1
    abs_Fe2 = max(abs(C_Fe(:) - C_Fe_prev(:)));
    rel_Fe2 = abs_Fe2 ./ max(max(abs(C_Fe_prev(:))), 1);
    % If Fe2 is low, absolute uM-scale changes should not dominate convergence.
    conv_Fe2 = min(rel_Fe2, abs_Fe2 ./ 25);
    conv_Fe3 = max(abs(FeooH(:) - FeooH_prev(:))) ./ max(max(abs(FeooH_prev(:))), 1);
    conv_SO4 = max(abs(Sulfate(:) - Sulfate_prev(:))) ./ max(max(abs(Sulfate_prev(:))), 1);
    conv_CH4 = max(abs(CH4(:) - CH4_prev(:))) ./ max(max(abs(CH4_prev(:))), 1);
    conv_ALK = max(abs(ALK(:) - ALK_prev(:))) ./ max(max(abs(ALK_prev(:))), 1);
    conv_pH  = max(abs(pH(:) - pH_prev(:))) ./ max(max(abs(pH_prev(:))), 1);
conv_FeS = max(abs(R_FeS(:) - R_FeS_prev(:))) ./ max(max(abs(R_FeS_prev(:))), 1);
abs_HS = max(abs(C_HS(:) - C_HS_prev(:)));
conv_HS_abs = abs_HS ./ 200;   % diagnostic only, not hard convergence
conv_all = max([conv_Fe2, conv_Fe3, conv_SO4, conv_CH4, conv_ALK, conv_pH, conv_FeS]);
fprintf('Outer iteration %d convergence: Fe2 %.3f, Fe3 %.3f, SO4 %.3f, HSabs %.3f, FeS %.3f, CH4 %.3f, ALK %.3f, pH %.3f, max %.3f\n', ...
    iteration, conv_Fe2, conv_Fe3, conv_SO4, conv_HS_abs, conv_FeS, conv_CH4, conv_ALK, conv_pH, conv_all);
else
    conv_all = Inf;
end
C_Fe_prev    = C_Fe;
FeooH_prev   = FeooH;
Sulfate_prev = Sulfate;
CH4_prev     = CH4;
ALK_prev     = ALK;
pH_prev      = pH;
C_HS_prev = C_HS;
R_FeS_prev = R_FeS;
if iteration >= min_outer_iter && conv_all < conv_tol
    break
end
iteration = iteration + 1;
count_loop = count_loop + 1;
end
%
% % -------------------------- Sulfur budget diagnostics --------------------------
% dz = z_sed(2) - z_sed(1);
%
% % Top boundary diffusive fluxes
% % Sign convention used here:
% %   J_SO4_top_down > 0 : downward sulfate flux INTO sediment
% %   J_HS_top_up   > 0 : upward sulfide flux OUT OF sediment
% J_SO4_top_down = -DSO4 .* ((Sulfate(2) - Sulfate(1)) ./ dz) .* 1e-3;   % umol/cm2/yr
% J_HS_top_up    =  DH2S .* ((C_HS(2)   - C_HS(1))   ./ dz) .* 1e-3;     % umol/cm2/yr
%
% % Depth-integrated rates (convert umol/L/yr over cm to umol/cm2/yr)
% I_SRR   = trapz(z_sed, R_SRR)     .* 1e-3;
% I_AOM   = trapz(z_sed, R_AOM_lag) .* 1e-3;
% I_FeS   = trapz(z_sed, R_FeS)     .* 1e-3;
% I_HSOx  = trapz(z_sed, R_HS_Ox)   .* 1e-3;
%
% % Irrigation terms as net source/sink diagnostics
% I_SO4_irrig = trapz(z_sed, Alpha_Bioirrig .* (SO4init - Sulfate)) .* 1e-3;   % >0 means irrigation supplies SO4
% I_HS_irrig  = trapz(z_sed, Alpha_Bioirrig .* (C_HS - HSinit))    .* 1e-3;    % >0 means irrigation removes HS
%
% % Net sulfur source/sink summaries
% SO4_total_sink   = I_SRR + I_AOM;
% SO4_total_supply = J_SO4_top_down + I_SO4_irrig;
%
% HS_total_source  = I_SRR + I_AOM;
% HS_total_sink    = I_FeS + I_HSOx + I_HS_irrig + J_HS_top_up;
%
% % Where is SRR happening?
% idx_deep10 = z_sed >= 10;
% idx_deep15 = z_sed >= 15;
%
% I_SRR_deep10 = trapz(z_sed(idx_deep10), R_SRR(idx_deep10)) .* 1e-3;
% I_SRR_deep15 = trapz(z_sed(idx_deep15), R_SRR(idx_deep15)) .* 1e-3;
%
% frac_SRR_deep10 = I_SRR_deep10 ./ max(I_SRR, 1e-12);
% frac_SRR_deep15 = I_SRR_deep15 ./ max(I_SRR, 1e-12);
%
% % Sulfide sink partition
% frac_FeS   = I_FeS   ./ max(HS_total_source, 1e-12);
% frac_HSOx  = I_HSOx  ./ max(HS_total_source, 1e-12);
% frac_irrig = I_HS_irrig ./ max(HS_total_source, 1e-12);
% frac_flux  = J_HS_top_up ./ max(HS_total_source, 1e-12);
%
% fprintf('\n================ Sulfur Budget Diagnostics ================\n');
% fprintf('Top SO4 diffusive influx (downward +):   %.3f umol/cm2/yr\n', J_SO4_top_down);
% fprintf('Integrated SO4 irrigation source:        %.3f umol/cm2/yr\n', I_SO4_irrig);
% fprintf('Integrated organoclastic SRR sink:       %.3f umol/cm2/yr\n', I_SRR);
% fprintf('Integrated AOM sulfate sink:             %.3f umol/cm2/yr\n', I_AOM);
% fprintf('SO4 total supply (top+irrigation):       %.3f umol/cm2/yr\n', SO4_total_supply);
% fprintf('SO4 total sink   (SRR+AOM):              %.3f umol/cm2/yr\n', SO4_total_sink);
%
% fprintf('\nIntegrated HS source (SRR+AOM):          %.3f umol/cm2/yr\n', HS_total_source);
% fprintf('Integrated FeS sink:                     %.3f umol/cm2/yr (%.2f%%)\n', I_FeS,   100*frac_FeS);
% fprintf('Integrated HS oxidation sink:            %.3f umol/cm2/yr (%.2f%%)\n', I_HSOx,  100*frac_HSOx);
% fprintf('Integrated HS irrigation sink:           %.3f umol/cm2/yr (%.2f%%)\n', I_HS_irrig, 100*frac_irrig);
% fprintf('Top HS diffusive efflux (upward +):      %.3f umol/cm2/yr (%.2f%%)\n', J_HS_top_up, 100*frac_flux);
% fprintf('HS total sink (FeS+Ox+irrig+flux):       %.3f umol/cm2/yr\n', HS_total_sink);
%
% fprintf('\nDeep SRR fraction (>10 cm):              %.2f%%\n', 100*frac_SRR_deep10);
% fprintf('Deep SRR fraction (>15 cm):              %.2f%%\n', 100*frac_SRR_deep15);
% fprintf('===========================================================\n\n');
% Storing steady-state solutions as initial conditions for PDEs
Diff_fluxes = [F_diff_DIC F_diff F_diff_DIC./F_diff]';
R_ALK_DIC = F_diff./F_diff_DIC;
F_diff_CH4 = DCH4.*((CH4(1,2) - CH4(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr
% % ----------------------------- OUTPUT PACKAGING ---------------------------
%     Outputs.z_sed = z_sed;
%     Outputs.pH_profile = pH;
%     Outputs.CH4_profile = CH4;
%     Outputs.O2_profile = Oxygen;
% %     Outputs.SO4_profile = Sulfate;
% %     Outputs.DIC_profile = C_DIC;
%
%     % Core Diagnostics
%     Outputs.Max_CH4 = max(CH4);
%
%
%     Outputs.Org_Bottom = C_organic(end) * 100; % %gDw
%     Outputs.ALK_Bottom = ALK(end);             % uM
%     Outputs.pH_Bottom  = pH(end);              %
%     Outputs.CH4_Bottom = CH4(end);             % uM
% %     Outputs.Org_Top    = C_organic(1) * 100; % %gDw
%
%     % OPD: O2 < 1 uM
%     idx_O2 = find(Oxygen < 1, 1);
%     if isempty(idx_O2), Outputs.OPD = z_sed(end); else, Outputs.OPD = z_sed(idx_O2); end
%
%     % SO4_Depth: SO4 降至 < 10 uM
%     idx_SO4 = find(Sulfate < 10, 1);
%     if isempty(idx_SO4), Outputs.SO4_Depth = z_sed(end); else, Outputs.SO4_Depth = z_sed(idx_SO4); end
%
%
%     idx_top5 = (z_sed <= 5);                   % 圈定 0-5 cm 网格
%     idx_bot5 = (z_sed >= (Lbottom - 5));       % 圈定底部 5 cm 网格
%
%     Outputs.ALK_Bot5   = mean(ALK(idx_bot5));             % 底层 5cm 平均碱度
%     Outputs.Sigma_Top5 = mean(sigma_carb(idx_top5));      % 表层 5cm 平均饱和度 (Omega-1)
%     Outputs.CaCO3_Top5 = mean(CaCO3(idx_top5)) * 100;     % 表层 5cm 平均 CaCO3 (%gDw)
%     Outputs.Integ_Meth = trapz(z_sed, Rate_Meth);         %integrated Rate_Meth
%
% %     % CH4_Onset_Depth: CH4 超过 10 uM 的深度
% %     idx_CH4 = find(CH4 > 10, 1);
% %     if isempty(idx_CH4), Outputs.CH4_Onset = z_sed(end); else, Outputs.CH4_Onset = z_sed(idx_CH4); end
% %
% %     % Sigma0_Depth: 碳酸钙饱和度 Omega-1 穿过 0 的深度 (>= 0)
% %     idx_sigma = find(sigma_carb >= 0, 1);
% %     if isempty(idx_sigma), Outputs.Sigma0_Depth = z_sed(end); else, Outputs.Sigma0_Depth = z_sed(idx_sigma); end
% %
% %     % CaCO3_Front_Depth: 碳酸钙开始显著积累的深度 (设定阈值为 1e-4，即脱离初始极小值)
% %     idx_CaCO3 = find(CaCO3 > 1e-4, 1);
% %     if isempty(idx_CaCO3), Outputs.CaCO3_Front = z_sed(end); else, Outputs.CaCO3_Front = z_sed(idx_CaCO3); end
%
%
% %     % Methane Appearance Depth (Depth where CH4 > 5 uM)
% %     ch4_idx = find(CH4 > 5, 1);
% %     if isempty(ch4_idx)
% %         Outputs.CH4_Depth = Lbottom; % No significant methane
% %     else
% %         Outputs.CH4_Depth = z_sed(ch4_idx);
% %     end
% %
% %     Outputs.Convergence_Status = K_converge;
%
% end % End of Function
% -------------------------------------------------------------------------
% ----------------------------- PLOTS -------------------------------------
clf;
n_plot = 6; % number of plots in each row
m_plot = 3; % number of total rows
% Organic
subplot(m_plot,n_plot,1);
% plot((C_organic + POC_root).*100,z_sed,'lineWidth',2); axis ij
plot((C_organic ).*100,z_sed,'lineWidth',2); axis ij
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
plot(Sulfate,z_sed,'lineWidth',2);
% xlim([0 SO4init]);
axis ij
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
% plot(CaCO3.*(1E-5.*(poros2./(1-poros2)).*(1./rho)),z_sed,'lineWidth',2); axis ij
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
title('OM Burial Efficiency')
ylabel('Depth (cm)');
box on
grid on
ax.LineWidth = 2;
% Aerobic respiration and sulfate reduction rates
% subplot(m_plot,n_plot,14);
% plot(R_SRR,z_sed,R_respi,z_sed,'lineWidth',2); axis ij %umol/l/year
% title('Rate (\mumol/l/year)')
% legend('Sulfate Red','Aerobic Resp');
subplot(m_plot,n_plot,15);
plot((0.5.*R_SRR)./365,z_sed,'lineWidth',2); axis ij %umol/l/year
title('Sulfate Reduction Rate (nmol/cm3/d)')
% legend('Sulfate Red');
box on
grid on
ax.LineWidth = 2;
% Carbonate saturation index
% subplot(m_plot,n_plot,15);
%
% plot(sigma_carb(iteration-1,:),z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Caclite saturation (\Omega - 1)')
%
% box on
subplot(m_plot,n_plot,17);
plot(sigma_carb(end,:),z_sed,'lineWidth',2); axis ij
title('Calciite saturation (\Omega - 1)')
box on
grid on
ax.LineWidth = 2;
subplot(m_plot,n_plot,18);
plot(FeooH(1,:),z_sed,'lineWidth',2); axis ij
title('Fe(III) (\mumol/gr)')
box on
grid on
ax.LineWidth = 2;
toc
% Mineral saturation indices
%
% subplot(m_plot,n_plot,16);
% plot(delta_viv1,z_sed,delta_apat2,z_sed,delta_FeS,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Mineral saturation (\Omega - 1)')
% legend('Vivianite','Apatite','FeS');
%
% box on
% grid on
%
% ax.LineWidth = 2;
% AAA_time = toc;
% AAA_data = [Oxygen' C_DIC' C_alka' pH' sigma_carb' z_sed'];
%
% R_iron_unit = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9.*...
%                        (poros./(1-poros)).*1E-3.*(1./rho); %rate of iron reduction umol/g/year
% R_FeOx_unit = (kFeOx.*C_Fe.*Oxygen).*...
%               (poros./(1-poros)).*1E-3.*(1./rho);
%
% R_SRR_integ = cumsum(0.5.*R_SRR.*dz_sed.*1E-3);
% R_RC_integ = cumsum(RC.*dz_sed);
% R_iron_integ = cumsum(R_iron(count_loop-1,:).*dz_sed.*1E-3);
% R_iron_integ_unit = cumsum(R_iron_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
% R_FeOx_integ_unit = cumsum(R_FeOx_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
% R_FeOx_integ = cumsum(R_FeOx_1(count_loop-1,:).*dz_sed.*1E-3);
% R_FeS_integ = cumsum(R_FeS.*dz_sed.*1E-3);
% R_FeS_1_integ = cumsum(R_FeS_1(count_loop-1,:).*dz_sed.*1E-3);
% R_HSOX_integ = cumsum(R_Ox.*dz_sed.*1E-3);
% R_biorrig_integ = cumsum((Alpha_Bioirrig.*(HSinit-C_HS)).*dz_sed.*1E-3);
% R_biorrig_integ_iron = cumsum((Alpha_Bioirrig.*(Feinit-C_Fe)).*dz_sed.*1E-3);
% R_biorrigALK_integ = cumsum((Alpha_Bioirrig.*(HCO3init-ALK)).*dz_sed.*1E-3);
% R_ALK = RC.*1E9 - R_respi;
% R_carb_integ = cumsum(R1_carb.*dz_sed.*1E-3);
% R_ALK_integ_WITHOUT = 2.*R_FeS_integ;
% R_ALK_integ_WITH = 2.*R_FeS_integ - (R_carb_integ);
% R_CH4O2_integ = cumsum((k_aerobic_CH4.* CH4.* (Oxygen./(Oxygen+K_CH4_O2))).*dz_sed.*1E-3);
% R_CH4SO4_integ = cumsum((k_AOM.* CH4.* (Sulfate./(Sulfate+K_CH4_SO4))).*dz_sed.*1E-3);
%
% R_net = R_SRR_integ - R_FeS_integ - R_HSOX_integ + R_biorrig_integ - F_diff_HS;
%
% R_net_iron = R_iron_integ - R_FeS_integ - R_FeOx_integ - F_diff_Fe(1,count_loop-1) + R_biorrig_integ_iron;
%
% R_net_Fe3  =  (v_burial(end).*C_Fe_3(end).*rho.*((1-poros(end))./poros(end))) - ...
%               (v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1))) + R_iron_integ_unit - R_FeOx_integ_unit;
%
% R_net_Fe3_percent  =  100.*R_net_Fe3(end)./(v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1)));
%
% R_net_org  =  (v_burial(end).*C_organic(end).*rho.*((1-poros(end))./12)) - (v_burial(1).*C_organic(1).*rho.*((1-poros(1))./12)) + ...
%               R_RC_integ;
%
% R_net_org_percent  =  100.*R_net_org(end)./((rho.*((1-poros(1))./12)).*v_burial(1).*C_organic(1));
%
% F_HS_tot = (R_SRR_integ - R_FeS_integ).*0.0274;
%
% F_S_out = (R_SRR_integ - R_FeS_integ - R_HSOX_integ).*0.0274;
%
% R_SRR_integ_store = R_SRR_integ(end).*0.0274;  %mmol/m2/d
% R_ALK_integ_WITH_store = R_ALK_integ_WITH(end).*0.0274; %mmol/m2/d
%
% F_ox_py = R_FeS_integ./R_SRR_integ;
%
% AAA_Store_1 = [F_FeOx R_SRR_integ_store R_ALK_integ_WITH_store];
% AAA_Store = [NPP.*BE R_ALK_integ_WITH(end) R_ALK_integ_WITHOUT(end) 2.*R_carb_integ(end) F_diff];
```

## File: Run_RTM_1D_PDE.m
```matlab
function Result = Run_RTM_1D_PDE()
% RUN_RTM_1D_PDE
% First coupled FV-MOL transient RTM driver.
%
% Purpose:
%   1. Run constant-forcing spin-up.
%   2. Test coupled RHS stability.
%   3. Generate profiles comparable to the stable ODE steady-state model.
clearvars -except Result
tic
Params = Params_Static();
Config = Config_Baseline();
% First PDE benchmark settings.
if ~isfield(Config, 't_spinup')
    Config.t_spinup = 300;  % yr
end
if ~isfield(Config, 'dt_out_spinup')
    Config.dt_out_spinup = 2;  % yr
end
Grid = Build_Grid_1D(Config, Params);
Forcing = Build_Forcing_PDE(Config);
State0 = Initialize_PDE_State(Grid, Forcing, Config, Params);
Y0 = Pack_State(State0);
tspan = 0:Config.dt_out_spinup:Config.t_spinup;
nonnegative_idx = 1:numel(Y0);
opts = odeset( ...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-8, ...
    'NonNegative', nonnegative_idx, ...
    'MaxStep', 1);
rhs = @(t,Y) RTM_PDE_RHS(t, Y, Grid, Forcing, Params, Config);
fprintf('Starting coupled FV-MOL PDE spin-up...\n');
[tout, Yout] = ode15s(rhs, tspan, Y0, opts);
fprintf('PDE spin-up finished.\n');
State_final = Unpack_State(Yout(end,:).', Grid);
State_final = floor_output_state(State_final);
[Rates_final, Diag_final] = RTM_Reaction_Rates(State_final, Grid, Forcing, Params, Config);
Result = struct();
Result.t = tout;
Result.Y = Yout;
Result.Grid = Grid;
Result.Forcing = Forcing;
Result.Params = Params;
Result.Config = Config;
Result.State_final = State_final;
Result.Rates_final = Rates_final;
Result.Diag_final = Diag_final;
Result.Summary = Summarize_PDE_Result(Result);
Result.Budget = RTM_Budget(Result);
% Plot_PDE_Result(Result);
toc
end
function State = floor_output_state(State)
names = fieldnames(State);
for i = 1:numel(names)
    name = names{i};
    if strcmp(name, 'DIC') || strcmp(name, 'ALK')
        State.(name) = max(real(State.(name)), 1e-12);
    else
        State.(name) = max(real(State.(name)), 0);
    end
end
end
```

## File: Sensitivity.m
```matlab

clc; clear; close all;
% 1. Define Parameter Space: {Name, BaseValue, MinVal, MaxVal}
% Grounded in realistic estuarine/riverine bounds
Param_Space = {
    'NPP',        200,   50,    600;   % Primary Production (g/m2/yr)
    'BE',         0.1,  0.03,  0.15;  % Burial Efficiency (fraction)
    'vbottom',    0.5,   0.1,   2.0;   % Sedimentation Rate (cm/yr)
    'F_FeOx',     2,   0.1,   5.0;   % Fe(III) Flux (mmol/m2/d)
    'F_CaCO3',        10,    1,     30;
    'SO4init',    200,   50,    500;  % Boundary SO4 (uM)
    'HCO3init',       950,   500,   2000;
    'DICinit',        1000,  700,   2500;
    'Calcium',        1000,  200,   2000;
    'Bioturbtop', 10,    1,     30;     % Bioturbation (cm2/yr)
    'k_SO4',      20,    5,     30;     % uM
    'DSO4',      310,   100,   500;    % cm2 / yr
    'DCH4',      300,   100,   500;    % cm2 / yr
    'DH2S',      300,   429,   650;    % cm2 / yr
    'DO2',      300,   370,   730;    % cm2 / yr
    'Kreox',      500,   10,   1000;    % 1 / umol / L / yr
    'kFeOx',      10,   1,   100;    % 1 / umol / L / yr
    'kFeS',       10,   1,   100;    % 1 / umol / L / yr
    'K_CH4_SO4',  100,   10,   500;    % uM
    'k_AOM',      1,   0.1,   10;    % 1 / yr
    'k_aerobic_CH4',      6,   1,   10 ;   % 1 / yr
    % 'k_calcite',      1,     0.2,   5;
    % 'k_calcite_dis1', 0.005, 0.001, 0.05;
    'k_calcite_dis2', 10,    1,     50;
};
num_params = size(Param_Space, 1);
range_fraction = 0.10; % Perturb by 10% of the total physical range
% 2. Execute Baseline Run
Base_Config = Config_Baseline();
Base_Params = Params_Static();
fprintf('Executing Baseline Run...\n');
try
    Base_Outputs = Run_RTM_1D(Base_Config, Base_Params);
catch
    error('Baseline run failed.');
end
% % Extract Baseline Targets
% Base_Max_CH4   = Base_Outputs.Max_CH4;
% Base_Bottom_pH = Base_Outputs.Bottom_pH;
%
% % Preallocate
% S_Max_CH4   = zeros(num_params, 1);
% S_Bottom_pH = zeros(num_params, 1);
% S_Org_Burial  = zeros(num_params, 1);
% S_ALK_Bottom  = zeros(num_params, 1);
% S_Carb_Burial = zeros(num_params, 1);
% S_CH4_Flux    = zeros(num_params, 1);
% Param_Names = cell(num_params, 1);
%
% % 3. Execute Range-Scaled Perturbation Loop
% fprintf('Starting Range-Scaled Scan (%.0f%% of feasible range)...\n', range_fraction * 100);
% tic;
%
% for i = 1:num_params
%     Param_Names{i} = Param_Space{i, 1};
%     base_val = Param_Space{i, 2};
%     min_val  = Param_Space{i, 3};
%     max_val  = Param_Space{i, 4};
%
%     % Calculate delta based on RANGE, not baseline
%     delta_X = range_fraction * (max_val - min_val);
%     perturb_val = base_val + delta_X;
%
%     Run_Config = Base_Config;
%     Run_Config.(Param_Names{i}) = perturb_val;
%
%     fprintf('Testing %s: %.3f -> %.3f ... ', Param_Names{i}, base_val, perturb_val);
%
%     try
%         Outputs = Run_RTM_1D(Run_Config);
%
%         % Calculate Range-Scaled Sensitivity (S)
%         delta_CH4 = Outputs.Max_CH4 - Base_Max_CH4;
%         % Mathematical safeguard for zero baseline methane (avoids Inf)
%         denom_CH4 = max(Base_Max_CH4, 1e-6);
%         S_Max_CH4(i) = (delta_CH4 / denom_CH4) / range_fraction;
%
%         delta_pH = Outputs.Bottom_pH - Base_Bottom_pH;
%         S_Bottom_pH(i) = (delta_pH / Base_Bottom_pH) / range_fraction;
%
%         fprintf('Done.\n');
%     catch ME
%         fprintf('FAILED (Stiff ODE). S assigned as NaN.\n');
%         S_Max_CH4(i)   = NaN;
%         S_Bottom_pH(i) = NaN;
%     end
% end
% exec_time = toc;
% fprintf('Scan Complete in %.2f seconds.\n', exec_time);
%
% % 4. Visualization (Tornado Plots)
% figure('Name', 'Range-Scaled Sensitivity Analysis', 'Color', 'w', 'Position', [150, 150, 1000, 450]);
%
% % Subplot 1: Sensitivity of Max CH4
% subplot(1,2,1);
% [sorted_S_CH4, idx_CH4] = sort(S_Max_CH4, 'ascend');
% sorted_Names_CH4 = Param_Names(idx_CH4);
%
% barh(sorted_S_CH4, 'FaceColor', [0.85 0.32 0.09], 'EdgeColor', 'k');
% set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_CH4, 'TickLabelInterpreter', 'none');
% xlabel('Range-Scaled Sensitivity (S_{range})');
% title('Sensitivity of Max CH_4');
% xline(0, 'k--', 'LineWidth', 1.5);
% grid on;
%
% % Subplot 2: Sensitivity of Bottom pH
% subplot(1,2,2);
% [sorted_S_pH, idx_pH] = sort(S_Bottom_pH, 'ascend');
% sorted_Names_pH = Param_Names(idx_pH);
%
% barh(sorted_S_pH, 'FaceColor', [0.0 0.44 0.74], 'EdgeColor', 'k');
% set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_pH, 'TickLabelInterpreter', 'none');
% xlabel('Range-Scaled Sensitivity (S_{range})');
% title('Sensitivity of Bottom pH');
% xline(0, 'k--', 'LineWidth', 1.5);
% grid on;
%
% sgtitle(sprintf('Range-Scaled Sensitivity (+%.0f%% of Physical Bound)', range_fraction*100), 'FontWeight', 'bold');
% --- 2. 预分配空间与基线提取 ---
Base_Org_Bottom = Base_Outputs.Org_Bottom;
% Base_Org_Top  = Base_Outputs.Org_Top;
Base_ALK = Base_Outputs.ALK_Bottom;
Base_pH  = Base_Outputs.pH_Bottom;
Base_CH4 = Base_Outputs.CH4_Bottom;
Base_OPD   = Base_Outputs.OPD;
Base_SO4D  = Base_Outputs.SO4_Depth;
% Base_CH4D  = Base_Outputs.CH4_Onset;
% Base_SigD  = Base_Outputs.Sigma0_Depth;
% Base_CaD   = Base_Outputs.CaCO3_Front;
B_Meth  = Base_Outputs.Integ_Meth;
B_ALK5  = Base_Outputs.ALK_Bot5;
B_Sig5  = Base_Outputs.Sigma_Top5;
B_Ca5   = Base_Outputs.CaCO3_Top5;
% 无论你有多少个 Output，都先建好全零数组
S_OrgB = zeros(num_params, 1);
% S_OrgT = zeros(num_params,1);
S_ALK = zeros(num_params, 1);
S_pH  = zeros(num_params, 1);
S_CH4 = zeros(num_params, 1);
S_OPD  = zeros(num_params,1);
S_SO4D = zeros(num_params,1);
% S_CH4D = zeros(num_params,1);
% S_SigD = zeros(num_params,1);
% S_CaD  = zeros(num_params,1);
S_ALK5 = zeros(num_params, 1);
S_Sig5 = zeros(num_params,1);
S_Ca5  = zeros(num_params,1);
S_Meth = zeros(num_params,1);
Param_Names = cell(num_params, 1);
% --- 3. 执行扰动循环 ---
fprintf('Starting Range-Scaled Scan...\n');
tic;
for i = 1:num_params
    Param_Names{i} = Param_Space{i, 1};
    base_val = Param_Space{i, 2};
    min_val  = Param_Space{i, 3};
    max_val  = Param_Space{i, 4};
    delta_X = range_fraction * (max_val - min_val);
    perturb_val = base_val + delta_X;
    Run_Config = Base_Config;
    Run_Params = Base_Params;
    % Auto-Route the parameter to the correct struct
    if isfield(Run_Config, Param_Names{i})
        Run_Config.(Param_Names{i}) = perturb_val;
    elseif isfield(Run_Params, Param_Names{i})
        Run_Params.(Param_Names{i}) = perturb_val;
    else
        error(['Parameter ', Param_Names{i}, ' not found in Config or Params.']);
    end
    fprintf('Testing %s: %.3f -> %.3f ... \n', Param_Names{i}, base_val, perturb_val);
    try
        Outputs = Run_RTM_1D(Run_Config, Run_Params); % Pass both!
        % 计算各个指标的归一化敏感度
        S_OrgB(i) = ((Outputs.Org_Bottom - Base_Org_Bottom) / max(Base_Org_Bottom, 1e-6)) / range_fraction;
%         S_OrgT(i) = ((Outputs.Org_Top    - Base_Org_Top) / max(Base_Org_Top, 1e-6)) / range_fraction;
        S_ALK(i) = ((Outputs.ALK_Bottom - Base_ALK) / max(Base_ALK, 1e-6)) / range_fraction;
        S_pH(i)  = ((Outputs.pH_Bottom  - Base_pH)  / max(Base_pH,  1e-6)) / range_fraction;
        S_CH4(i) = ((Outputs.CH4_Bottom - Base_CH4) / max(Base_CH4, 1e-6)) / range_fraction;
        S_ALK5(i) = ((Outputs.ALK_Bot5   - B_ALK5) / max(B_ALK5, 1e-6)) / range_fraction;
        S_Sig5(i) = ((Outputs.Sigma_Top5 - B_Sig5) / max(abs(B_Sig5), 1e-4)) / range_fraction;
        S_Ca5(i)  = ((Outputs.CaCO3_Top5 - B_Ca5)  / max(B_Ca5,  1e-6)) / range_fraction;
        S_Meth(i) = ((Outputs.Integ_Meth - B_Meth) / max(B_Meth, 1e-6)) / range_fraction;
        S_OPD(i)  = ((Outputs.OPD       - Base_OPD)  / max(Base_OPD,  0.1)) / range_fraction;
        S_SO4D(i) = ((Outputs.SO4_Depth - Base_SO4D) / max(Base_SO4D, 0.1)) / range_fraction;
%         S_CH4D(i) = ((Outputs.CH4_Onset - Base_CH4D) / max(Base_CH4D, 0.1)) / range_fraction;
%         S_SigD(i) = ((Outputs.Sigma0_Depth- Base_SigD) / max(Base_SigD, 0.1)) / range_fraction;
%         S_CaD(i)  = ((Outputs.CaCO3_Front- Base_CaD)  / max(Base_CaD,  0.1)) / range_fraction;
    catch ME
        S_OrgB(i) = NaN; S_ALK(i) = NaN; S_pH(i) = NaN; S_CH4(i) = NaN;
        fprintf('报错信息为: %s\n', ME.message);
    end
end
fprintf('Scan Complete in %.2f seconds.\n', toc);
% --- 4. 可扩展动态绘图模块 ---
% 【扩展指南】未来若要增加输出，只需在这两个 Cell Array 中添加新变量和标题即可
Targets = {S_OrgB, S_ALK, S_pH, S_CH4, S_OPD, S_SO4D, S_Meth, S_Sig5, S_Ca5};
Titles  = {'Bottom Organic (%)', 'Bottom ALK (\muM)', 'Bottom pH', 'Bottom CH_4 (\muM)', ...
        'OPD (O_2 < 1\muM)', 'SO_4 Depletion Depth','Integrated Methanogenesis',...
         'Mean \Omega-1 (Top 5cm)', 'Mean CaCO_3 (Top 5cm)'};
n_targets = length(Targets);
% 自动计算子图的行列数
cols = ceil(sqrt(n_targets));
rows = ceil(n_targets / cols);
% 自动调整窗口大小以适应子图数量
figure('Name', 'Multi-Output Sensitivity', 'Color', 'w', 'Position', [100, 100, cols*400, rows*350]);
for k = 1:n_targets
    subplot(rows, cols, k);
%     [sorted_S, idx] = sort(Targets{k}, 'ascend');
%     sorted_Names = Param_Names(idx);
    % 使用统一的配色
%      barh(sorted_S, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
barh(Targets{k}, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', Param_Names, 'TickLabelInterpreter', 'none');
set(gca, 'YDir', 'reverse');
%      set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names, 'TickLabelInterpreter', 'none');
%     xlabel('Sensitivity Index (S_{range})');
    title(['Sensitivity of: ', Titles{k}]);
    xline(0, 'k-', 'LineWidth', 1.2);
    grid on;
    max_abs_val = max(abs(Targets{k}(~isnan(Targets{k}))));
    % 2. 安全保护：如果算出全是 0 或报错，给一个默认基准 1
    if isempty(max_abs_val) || max_abs_val == 0
        max_abs_val = 1;
    end
    % 3. 强制锁定坐标轴边界，乘 1.1 是为了左右两边留出 10% 的空白余量好看
    xlim([-max_abs_val * 1.1, max_abs_val * 1.1]);
end
```

## File: SO4_bc.m
```matlab

function res = SO4_bc(SO4a,SO4b)
global SO4init
  res = [ SO4a(1)-SO4init
          SO4b(2) ];
end
```

## File: SO4_ODE.m
```matlab

function dydx = SO4_ODE(x, SO4)
global k_SO4 K_CH4_SO4 RC_after_Fe DSO4 v_burial_Fluid SO4init Alpha_Bioirrig z_sed poros
global R_AOM_pot
v_burial_f = interp1(z_sed, v_burial_Fluid, x, 'linear', 'extrap');
Alpha      = interp1(z_sed, Alpha_Bioirrig, x, 'linear', 'extrap');
fi         = interp1(z_sed, poros, x, 'linear', 'extrap');
RC1        = interp1(z_sed, RC_after_Fe, x, 'linear', 'extrap');
% potential AOM from previous CH4 profile
R_AOM_pot_1 = interp1(z_sed, R_AOM_pot, x, 'linear', 'extrap');
% True nonnegative sulfate for reaction limitation.
% Important: no softplus. If SO4 <= 0, reactions using SO4 shut down.
SO4_pos = max(real(SO4(1)), 0);
f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
R_SRR_local = 0.5 .* RC1 .* f_SRR;
R_AOM_local = R_AOM_pot_1 .* f_AOM;
NR = + v_burial_f .* (SO4(2) ./ (fi .* DSO4)) ...
     + R_SRR_local ...
     + R_AOM_local ...
     - Alpha .* (SO4init - SO4(1));
dydx = [ SO4(2) ./ (fi .* DSO4)
         NR ];
end
```

## File: Solve_SO4_FV_Front.m
```matlab
function [Sulfate, R_SRR, R_AOM_actual, SO4_diag] = Solve_SO4_FV_Front()
% Positivity-preserving finite-volume SO4 solver with depletion front behavior.
%
% This replaces the ordinary SO4 bvp4c solve.
% It keeps SO4 >= 0 and makes SRR/AOM vanish in sulfate-depleted cells.
global z_sed poros DSO4 SO4init Alpha_Bioirrig
global k_SO4 K_CH4_SO4 RC_after_Fe R_AOM_pot
z = z_sed(:);
n = numel(z);
if n < 3
    error('Solve_SO4_FV_Front requires at least 3 grid cells.');
end
dz_all = diff(z);
if max(abs(dz_all - dz_all(1))) > 1e-10
    error('Solve_SO4_FV_Front currently assumes uniform z_sed.');
end
dz = dz_all(1);
phi = poros(:);
Alpha = Alpha_Bioirrig(:);
Knode = max(phi .* DSO4, 1e-12);
Kface = 0.5 .* (Knode(1:end-1) + Knode(2:end));
R_SRR_pot = max(0, 0.5 .* RC_after_Fe(:));  % umol SO4/L/yr
R_AOM_pot_col = max(0, R_AOM_pot(:));       % umol SO4/L/yr potential
S_floor = 0;
S_cut = 1e-6;
tol = 1e-8;
max_iter = 200;
relax = 0.7;
% Initial guess: positive, decreasing sulfate.
S = SO4init .* exp(-z ./ max(0.25*z(end), dz));
S(1) = SO4init;
for it = 1:max_iter
    S_old = max(S, S_cut);
    % Linearized Monod sinks:
    % R ~= lambda * S
    lambda_srr = R_SRR_pot ./ max(S_old + k_SO4, 1e-12);
    lambda_aom = R_AOM_pot_col ./ max(S_old + K_CH4_SO4, 1e-12);
    lambda = lambda_srr + lambda_aom;
    % Depleted cells remain non-reactive.
    depleted_old = S <= S_cut;
    lambda(depleted_old) = 0;
    A = spalloc(n, n, 3*n);
    b = zeros(n,1);
    % Top Dirichlet: bottom-water sulfate.
    A(1,1) = 1;
    b(1) = SO4init;
    for i = 2:n
        Km = Kface(i-1);
        if i < n
            Kp = Kface(i);
        else
            Kp = 0;  % bottom no-flux
        end
        % Equation:
        % -div(K grad S) + (Alpha + lambda)*S = Alpha*SO4init
        A(i,i-1) = -Km / dz^2;
        A(i,i)   = (Km + Kp) / dz^2 + Alpha(i) + lambda(i);
        if i < n
            A(i,i+1) = -Kp / dz^2;
        end
        b(i) = Alpha(i) .* SO4init;
    end
    S_lin = A \ b;
    S_new = max(S_floor, relax .* S_lin + (1 - relax) .* S);
    S_new(1) = SO4init;
    err = norm(S_new - S, inf) ./ max(SO4init, 1);
    S = S_new;
    if err < tol
        break
    end
end
% Active-set depletion projection.
S(S <= S_cut) = 0;
S(1) = SO4init;
% Actual rates from final nonnegative sulfate.
f_SRR = S ./ max(S + k_SO4, 1e-12);
f_AOM = S ./ max(S + K_CH4_SO4, 1e-12);
f_SRR(S <= S_cut) = 0;
f_AOM(S <= S_cut) = 0;
R_SRR_col = R_SRR_pot .* f_SRR;
R_AOM_col = R_AOM_pot_col .* f_AOM;
% Top node is a boundary value, not a reactive volume.
R_SRR_col(1) = 0;
R_AOM_col(1) = 0;
Sulfate = S(:).';
R_SRR = R_SRR_col(:).';
R_AOM_actual = R_AOM_col(:).';
idx_front = find(Sulfate <= 1, 1, 'first');
if isempty(idx_front)
    front_depth = z_sed(end);
else
    front_depth = z_sed(idx_front);
end
SO4_diag.iter = it;
SO4_diag.min_SO4 = min(Sulfate);
SO4_diag.front_depth = front_depth;
SO4_diag.I_SRR = trapz(z_sed, R_SRR) .* 1e-3;
SO4_diag.I_AOM = trapz(z_sed, R_AOM_actual) .* 1e-3;
SO4_diag.I_demand_pot = trapz(z_sed, R_SRR_pot(:).') .* 1e-3;
SO4_diag.I_irrig_source = trapz(z_sed, Alpha_Bioirrig .* (SO4init - Sulfate)) .* 1e-3;
SO4_diag.J_top_down = -DSO4 .* poros(1) .* ((Sulfate(2) - Sulfate(1)) ./ dz) .* 1e-3;
end
```

## File: Summarize_PDE_Result.m
```matlab
function Summary = Summarize_PDE_Result(Result)
% SUMMARIZE_PDE_RESULT
% Print compact quantitative diagnostics for the FV-MOL PDE result.
Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;
z = Grid.z(:);
Summary = struct();
Summary.OPD = first_depth_below(z, S.O2, 1, z(end));
Summary.SO4_Depth_10 = first_depth_below(z, S.SO4, 10, z(end));
Summary.SO4_Depth_1  = first_depth_below(z, S.SO4, 1, z(end));
Summary.OM_top_pct = 100 .* (S.OM_lab(1) + S.OM_ref(1));
Summary.OM_bottom_pct = 100 .* (S.OM_lab(end) + S.OM_ref(end));
Summary.O2_bottom = S.O2(end);
Summary.SO4_bottom = S.SO4(end);
Summary.Fe2_max = max(S.Fe2);
Summary.Fe2_bottom = S.Fe2(end);
Summary.HS_max = max(S.HS);
Summary.CH4_max = max(S.CH4);
Summary.CH4_bottom = S.CH4(end);
Summary.FeOOH_top = S.FeOOH(1);
Summary.FeOOH_max = max(S.FeOOH);
Summary.FeOOH_bottom = S.FeOOH(end);
Summary.DIC_bottom = S.DIC(end);
Summary.ALK_bottom = S.ALK(end);
Summary.pH_bottom = D.pH(end);
Summary.sigma_top5 = mean(D.sigma_carb(z <= 5));
Summary.sigma_bottom = D.sigma_carb(end);
Summary.CaCO3_top_pct = 100 .* S.CaCO3(1);
Summary.CaCO3_bottom_pct = 100 .* S.CaCO3(end);
Summary.CaCO3_net_int = D.I_CaCO3_net;
Summary.I_RC = D.I_RC;
Summary.I_respi = D.I_respi;
Summary.I_FeRed_C = D.I_FeRed_C;
Summary.I_SRR = D.I_SRR;
Summary.I_AOM = D.I_AOM;
Summary.I_Meth = D.I_Meth;
Summary.redox_closure = ...
    (D.I_respi + D.I_FeRed_C + 2 .* D.I_SRR + 2 .* D.I_Meth) ./ ...
    max(D.I_RC, 1e-12);
Summary.Fe_supply_scale = D.Fe_supply_scale;
Summary.I_FeRed_pot = D.I_FeRed_pot;
Summary.I_Fe_supply_ext = D.I_Fe_supply_ext;
Summary.last_step_change = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    S_now  = S;
    names = {'O2','Fe2','SO4','HS','CH4','DIC','ALK','FeOOH','CaCO3'};
    changes = zeros(numel(names),1);
    for i = 1:numel(names)
        name = names{i};
        a = S_now.(name)(:);
        b = S_prev.(name)(:);
        changes(i) = max(abs(a - b)) ./ max(max(abs(b)), 1);
    end
    Summary.last_step_change = max(changes);
    Summary.last_step_change_by_species = array2table(changes(:).', ...
        'VariableNames', names);
end
fprintf('\n================ PDE Summary ================\n');
fprintf('OPD <1 uM:                 %.2f cm\n', Summary.OPD);
fprintf('SO4 <10 uM depth:          %.2f cm\n', Summary.SO4_Depth_10);
fprintf('SO4 bottom:                %.2f uM\n', Summary.SO4_bottom);
fprintf('\nOM top / bottom:            %.3f / %.3f %%gDW\n', ...
    Summary.OM_top_pct, Summary.OM_bottom_pct);
fprintf('\nFe2 max / bottom:           %.2f / %.2f uM\n', ...
    Summary.Fe2_max, Summary.Fe2_bottom);
fprintf('FeOOH top / max / bottom:   %.2f / %.2f / %.2f umol/g\n', ...
    Summary.FeOOH_top, Summary.FeOOH_max, Summary.FeOOH_bottom);
fprintf('Fe supply scale:            %.3f\n', Summary.Fe_supply_scale);
fprintf('\nHS max:                     %.2f uM\n', Summary.HS_max);
fprintf('CH4 max / bottom:           %.2f / %.2f uM\n', ...
    Summary.CH4_max, Summary.CH4_bottom);
fprintf('\nDIC bottom:                 %.2f uM\n', Summary.DIC_bottom);
fprintf('ALK bottom:                 %.2f uM\n', Summary.ALK_bottom);
fprintf('pH bottom:                  %.3f\n', Summary.pH_bottom);
fprintf('Mean sigma top 5 cm:        %.4f\n', Summary.sigma_top5);
fprintf('CaCO3 top / bottom:         %.4f / %.4f %%gDW\n', ...
    Summary.CaCO3_top_pct, Summary.CaCO3_bottom_pct);
fprintf('\nIntegrated RC:              %.3f umol C/cm2/yr\n', Summary.I_RC);
fprintf('Integrated O2 resp:         %.3f\n', Summary.I_respi);
fprintf('Integrated Fe reduction C:  %.3f\n', Summary.I_FeRed_C);
fprintf('Integrated SRR:             %.3f umol SO4/cm2/yr\n', Summary.I_SRR);
fprintf('Integrated AOM:             %.3f umol SO4/cm2/yr\n', Summary.I_AOM);
fprintf('Integrated methanogenesis:  %.3f umol CH4/cm2/yr\n', Summary.I_Meth);
fprintf('Redox closure diagnostic:   %.3f\n', Summary.redox_closure);
if ~isnan(Summary.last_step_change)
    fprintf('\nMax relative change last output step: %.3e\n', Summary.last_step_change);
    disp(Summary.last_step_change_by_species);
end
fprintf('=============================================\n\n');
end
function depth = first_depth_below(z, x, threshold, default_depth)
idx = find(x(:) < threshold, 1, 'first');
if isempty(idx)
    depth = default_depth;
else
    depth = z(idx);
end
end
```

## File: Unpack_State.m
```matlab
function State = Unpack_State(Y, Grid)
% UNPACK_STATE
% Convert ode15s vector back to state struct.
n = Grid.n;
Y = Y(:);
expected_len = 11 * n;
if numel(Y) ~= expected_len
    error('State vector length mismatch. Expected %d, got %d.', expected_len, numel(Y));
end
i1 = 1;
State.OM_lab = Y(i1:i1+n-1); i1 = i1 + n;
State.OM_ref = Y(i1:i1+n-1); i1 = i1 + n;
State.FeOOH  = Y(i1:i1+n-1); i1 = i1 + n;
State.CaCO3  = Y(i1:i1+n-1); i1 = i1 + n;
State.O2  = Y(i1:i1+n-1); i1 = i1 + n;
State.Fe2 = Y(i1:i1+n-1); i1 = i1 + n;
State.SO4 = Y(i1:i1+n-1); i1 = i1 + n;
State.HS  = Y(i1:i1+n-1); i1 = i1 + n;
State.CH4 = Y(i1:i1+n-1); i1 = i1 + n;
State.DIC = Y(i1:i1+n-1); i1 = i1 + n;
State.ALK = Y(i1:i1+n-1);
end
```

## File: untitled.m
```matlab

clc; clear; close all;
% 1. Define Parameter Space: {Name, BaseValue, MinVal, MaxVal}
% Grounded in realistic estuarine/riverine bounds
Param_Space = {
    'NPP',        200,   50,    600;   % Primary Production (g/m2/yr)
    'BE',         0.05,  0.01,  0.20;  % Burial Efficiency (fraction)
    'vbottom',    0.5,   0.1,   2.0;   % Sedimentation Rate (cm/yr)
    'F_FeOx',     1.0,   0.1,   5.0;   % Fe(III) Flux (mmol/m2/d)
    'SO4init',    200,   50,    1000;  % Boundary SO4 (uM)
    'Bioturbtop', 10,    1,     30     % Bioturbation (cm2/yr)
};
num_params = size(Param_Space, 1);
range_fraction = 0.10; % Perturb by 10% of the total physical range
% 2. Execute Baseline Run
Base_Config = Config_Baseline();
fprintf('Executing Baseline Run...\n');
try
    Base_Outputs = Run_RTM_1D(Base_Config);
catch
    error('Baseline run failed. Ensure Run_RTM_1D is stable.');
end
% Extract Baseline Targets
Base_Max_CH4   = Base_Outputs.Max_CH4;
Base_Bottom_pH = Base_Outputs.Bottom_pH;
% Preallocate
S_Max_CH4   = zeros(num_params, 1);
S_Bottom_pH = zeros(num_params, 1);
Param_Names = cell(num_params, 1);
% 3. Execute Range-Scaled Perturbation Loop
fprintf('Starting Range-Scaled Scan (%.0f%% of feasible range)...\n', range_fraction * 100);
tic;
for i = 1:num_params
    Param_Names{i} = Param_Space{i, 1};
    base_val = Param_Space{i, 2};
    min_val  = Param_Space{i, 3};
    max_val  = Param_Space{i, 4};
    % Calculate delta based on RANGE, not baseline
    delta_X = range_fraction * (max_val - min_val);
    perturb_val = base_val + delta_X;
    Run_Config = Base_Config;
    Run_Config.(Param_Names{i}) = perturb_val;
    fprintf('Testing %s: %.3f -> %.3f ... ', Param_Names{i}, base_val, perturb_val);
    try
        Outputs = Run_RTM_1D(Run_Config);
        % Calculate Range-Scaled Sensitivity (S)
        delta_CH4 = Outputs.Max_CH4 - Base_Max_CH4;
        % Mathematical safeguard for zero baseline methane (avoids Inf)
        denom_CH4 = max(Base_Max_CH4, 1e-6);
        S_Max_CH4(i) = (delta_CH4 / denom_CH4) / range_fraction;
        delta_pH = Outputs.Bottom_pH - Base_Bottom_pH;
        S_Bottom_pH(i) = (delta_pH / Base_Bottom_pH) / range_fraction;
        fprintf('Done.\n');
    catch ME
        fprintf('FAILED (Stiff ODE). S assigned as NaN.\n');
        S_Max_CH4(i)   = NaN;
        S_Bottom_pH(i) = NaN;
    end
end
exec_time = toc;
fprintf('Scan Complete in %.2f seconds.\n', exec_time);
% 4. Visualization (Tornado Plots)
figure('Name', 'Range-Scaled Sensitivity Analysis', 'Color', 'w', 'Position', [150, 150, 1000, 450]);
% Subplot 1: Sensitivity of Max CH4
subplot(1,2,1);
[sorted_S_CH4, idx_CH4] = sort(S_Max_CH4, 'ascend');
sorted_Names_CH4 = Param_Names(idx_CH4);
barh(sorted_S_CH4, 'FaceColor', [0.85 0.32 0.09], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_CH4, 'TickLabelInterpreter', 'none');
xlabel('Range-Scaled Sensitivity (S_{range})');
title('Sensitivity of Max CH_4');
xline(0, 'k--', 'LineWidth', 1.5);
grid on;
% Subplot 2: Sensitivity of Bottom pH
subplot(1,2,2);
[sorted_S_pH, idx_pH] = sort(S_Bottom_pH, 'ascend');
sorted_Names_pH = Param_Names(idx_pH);
barh(sorted_S_pH, 'FaceColor', [0.0 0.44 0.74], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_pH, 'TickLabelInterpreter', 'none');
xlabel('Range-Scaled Sensitivity (S_{range})');
title('Sensitivity of Bottom pH');
xline(0, 'k--', 'LineWidth', 1.5);
grid on;
sgtitle(sprintf('Range-Scaled Sensitivity (+%.0f%% of Physical Bound)', range_fraction*100), 'FontWeight', 'bold');
```

