# MATLAB Model Source Code

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
SO4 = interp1(z_sed,Sulfate,x);
% Inh = (k_O2./(O2+k_O2));
% Inh1 = (k_SO4./(SO4+k_SO4));
RC1 = interp1(z_sed,RC,x);
R_Meth_current = double(interp1(z_sed, Rate_Meth, x));
% NR = + v_burial_f.* CH4(2) - 0.5.*RC1.*Inh.*Inh1.* 1E9 - (Alpha_Bioirrig_1.*(CH4init-CH4(1)))...
%       + k_AOM.* CH4(1).* (SO4./(SO4+K_CH4_SO4)) + k_aerobic_CH4.* CH4(1).* (O2./(O2+K_CH4_O2)); % umol/l/year
NR = + v_burial_f.* (CH4(2)/(fi*DCH4)) - R_Meth_current - (Alpha_Bioirrig_1.*(CH4init-CH4(1)))...
      + k_AOM.* CH4(1).* (SO4./(SO4+K_CH4_SO4)) + k_aerobic_CH4.* CH4(1).* (O2./(O2+K_CH4_O2)); % umol/l/year
dydx = [ CH4(2) /fi/DCH4
           NR];
end
```

## File: Config_Baseline.m
```matlab
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
    Config.BE  = 0.05;
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
```

## File: Coupled_Carbonate_ODE.m
```matlab
function dYdx = Coupled_Carbonate_ODE(x, Y)
% Y(1) = DIC, Y(2) = dDIC/dx flux
% Y(3) = ALK, Y(4) = dALK/dx flux
% Y(5) = CaCO3, Y(6) = dCaCO3/dx (dummy/solid flux)
global DHCO3 RC R_SRR kFeS C_Fe C_HS Alpha_Bioirrig DICinit HCO3init
global v_burial_Fluid v_burial z_sed poros rho
global k_calcite k_calcite_dis1 k_calcite_dis2 n_power_CaCO31 n_power_CaCO32 n_power_CaCO33
global Calcium Calcium_activity CO3_activity Ksp_ca
global Rate_Meth T_future R_HS_Ox
    v_burial_f = double(interp1(z_sed, v_burial_Fluid, x));
    v_burial_s = double(interp1(z_sed, v_burial, x));
    fi = double(interp1(z_sed, poros, x));
    Alpha_Bioirrig_1 = double(interp1(z_sed, Alpha_Bioirrig, x));
    RC1 = double(interp1(z_sed, RC, x));
    R_SRR1 = double(interp1(z_sed, R_SRR, x));
    R_HS_Ox_1 = double(interp1(z_sed, R_HS_Ox, x));
    C_Fe_1 = double(interp1(z_sed, C_Fe, x));
    C_HS_1 = double(interp1(z_sed, C_HS, x));
    DIC = max(real(Y(1)), 1e-12);
    ALK = max(real(Y(3)), 1e-12);
    CaCO3 = max(real(Y(5)), 0);
    % calculate real time CO3
    [~, CO3_current, ~] = River_Carbonate(ALK, DIC, T_future, 0.1, 1);
    if isempty(CO3_current) || isnan(CO3_current)
        CO3_current = 1e-12;
    end
    % calculate real time sigma_carb (Omega)
    sigma_carb = (Calcium_activity * Calcium * CO3_activity * CO3_current) / Ksp_ca - 1;
    % calculate real time CaCO3
    unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho * (1 - fi));
    R_carb_form = (sigma_carb > 0) * abs(sigma_carb)^n_power_CaCO31 * k_calcite * unit_conversion;
    R_carb_disso = (-0.2 < sigma_carb & sigma_carb < 0) * abs(sigma_carb)^n_power_CaCO32 * k_calcite_dis1 * CaCO3 ...
                 + (sigma_carb <= -0.2) * abs(sigma_carb)^n_power_CaCO33 * k_calcite_dis2 * CaCO3;
    R1_carb_total = R_carb_form - R_carb_disso * (1E3 * 1E6 * 1E-2 * rho * (1 - fi));
    % account for DIC loss in methanogenesis
    R_Meth_current = double(interp1(z_sed, Rate_Meth, x));
Advection_DIC = v_burial_f .* (Y(2) / (fi * DHCO3));
    % Advection_DIC = v_burial_f .* Y(2) ;
    NR_DIC = Advection_DIC - (RC1*1E9-R_Meth_current) + R1_carb_total - (Alpha_Bioirrig_1*(DICinit - DIC));
    Advection_ALK = v_burial_f .* (Y(4) / (fi * DHCO3));
%     NR_ALK = Advection_ALK - 2*(kFeS*C_Fe_1*C_HS_1) + 2*R1_carb_total - (Alpha_Bioirrig_1*(HCO3init - ALK));
    eta_s = 0.5;   % 先试 0.5，再看是否要到 1.0
    Net_S_term = eta_s * ( R_SRR1 ...
                         - 2*(kFeS*C_Fe_1*C_HS_1) ...
                         - R_HS_Ox_1 );
    NR_ALK = Advection_ALK +Net_S_term - 2*(kFeS*C_Fe_1*C_HS_1) + 2*R1_carb_total - (Alpha_Bioirrig_1*(HCO3init - ALK));
    NR_CaCO3 = R_carb_form - R_carb_disso;
    %
    % dYdx = [ Y(2) / (fi * DHCO3);
    %          NR_DIC;
    %          Y(4) / (fi * DHCO3);
    %          NR_ALK;
    %          NR_CaCO3 / v_burial_s;
    %          0 ];
    dYdx = zeros(6,1);
    dYdx(1) = Y(2) / (fi * DHCO3);
    dYdx(2) = NR_DIC;
    dYdx(3) = Y(4) / (fi * DHCO3);
    dYdx(4) = NR_ALK;
    dYdx(5) = NR_CaCO3 / v_burial_s;
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
global z_sed Bioturb RC
global rho Oxygen KFEMonod k_O2 v_burial C_Fe kFeOx poros
RC_1       = interp1(z_sed,RC,x);
poros_1    = interp1(z_sed,poros,x);
O2         = interp1(z_sed,Oxygen,x);
C_Fe_1     = interp1(z_sed,C_Fe,x);
Inh        = (k_O2./(O2+k_O2));
Db         = interp1(z_sed,Bioturb,x);
v_burial_1 = interp1(z_sed,v_burial,x);
sigh = 1 - poros_1;
NR = - 4.*RC_1.*Inh.*1E9.*(Fe3(1)./(Fe3(1)+KFEMonod)).*(poros_1./(1-poros_1)).*1E-3.*(1./rho) ...
     + (kFeOx.*C_Fe_1.*O2).*(poros_1./(1-poros_1)).*1E-3.*(1./rho);   % umol/g/yr
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
global RC Oxygen k_O2 DH2S v_burial_Fluid Alpha_Bioirrig z_sed poros Feinit
global FeooH kFeOx KFEMonod R_FeS kFeS C_HS
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
C_HS_1 = interp1(z_sed,C_HS,x);
FeOx = interp1(z_sed,FeooH,x);
Inh = (k_O2./(O2+k_O2));
% R_FeS_1 = interp1(z_sed,R_FeS,x);
RC1 = interp1(z_sed,RC,x);
% NR = + v_burial_f.* Fe2(2) - 4.*RC1.*Inh.* 1E9.* (FeOx./(FeOx+KFEMonod)) - (Alpha_Bioirrig_1.*(Feinit-Fe2(1)))...
%      + (kFeS.*Fe2(1).*C_HS_1) + (kFeOx.*Fe2(1).*O2); % umol/l/year
NR = + v_burial_f.* (Fe2(2)/(fi*DH2S)) - 4.*RC1.*Inh.* 1E9.* (FeOx./(FeOx+KFEMonod)) - (Alpha_Bioirrig_1.*(Feinit-Fe2(1)))...
     + (kFeS.*Fe2(1).*C_HS_1) + (kFeOx.*Fe2(1).*O2); % umol/l/year
dydx = [ Fe2(2) /fi/DH2S
           NR];
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
function dydx = HS_ODE(x,H2S)
global k_SO4 RC Oxygen k_O2 v_burial_Fluid HSinit
global Kreox Alpha_Bioirrig z_sed poros Sulfate DH2S R_FeS C_Fe kFeS
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
SO4 = interp1(z_sed,Sulfate,x);
C_Fe_1 = interp1(z_sed,C_Fe,x);
Inh = (k_O2./(O2+k_O2));
% R_FeS_1 = interp1(z_sed,R_FeS,x);
RC1 = interp1(z_sed,RC,x);
% NR = + v_burial_f.* H2S(2) - 0.5.*RC1.*Inh.* (SO4/(SO4+k_SO4)).* 1E9 ...
%      + (kFeS.*C_Fe_1.*H2S(1)) + (Kreox.*H2S(1).*O2) - (Alpha_Bioirrig_1.*(HSinit-H2S(1))); % umol/l/year
NR = + v_burial_f.* (H2S(2)/(fi*DH2S)) - 0.5.*RC1.*Inh.* (SO4/(SO4+k_SO4)).* 1E9 ...
     + (kFeS.*C_Fe_1.*H2S(1)) + (Kreox.*H2S(1).*O2) - (Alpha_Bioirrig_1.*(HSinit-H2S(1))); % umol/l/year
dydx = [ H2S(2) /fi/DH2S
           NR];
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

## File: Params_Static.m
```matlab
function Params = Params_Static()
% PARAMS_STATIC
% Constants that should not change from site to site unless explicitly audited.
    Params.rho = 2.73;           % g / cm3
    Params.k_O2 = 2;             % uM
    Params.k_SO4 = 20;           % uM
    Params.KFEMonod = 200;       % umol / g
    Params.DSO4 = 300;           % cm2 / yr
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

## File: River_Carbonate.m
```matlab
function [pH, CO3, H2CO3] = River_Carbonate(ALK, DIC, T, S, P)
% A unified carbonate thermodynamic solver for USGS river modifications
% Input: ALK (umol/kg), DIC (umol/kg), T ( °C), S (salinity ppt), P (pressure at depth atm)
    % If no input, use default value
    if nargin < 3
        T = 20;
        S = 0.1;
        P = 1;
    end
    T_K = T + 273.15;
    B_T = 400 * (S / 35); % total boron in river [umol/kg]
    % Calculate K
    K0 = exp((9345.17/T_K) - 60.2409 + 23.3585*log(T_K/100) + S*(0.023517 - 0.00023656*T_K + 0.00000047036*T_K^2));
    RGAS = 8.314510; R = 83.131;
    % K1
    lnK1 = 2.83655 - 2307.1266/T_K - 1.5529413*log(T_K) - (0.20760841 + 4.0484/T_K)*sqrt(S) + 0.08468345*S - 0.00654208*S^1.5 + log(1 - 0.001005*S);
    K1 = exp(lnK1) * 1e6;
    % K2
    lnK2 = -9.226508 - 3351.6106/T_K - 0.2005743*log(T_K) + (-0.106901773 - 23.9722/T_K)*sqrt(S) + 0.1130822*S - 0.00846934*S^1.5 + log(1 - 0.001005*S);
    K2 = exp(lnK2) * 1e6;
    % Kb
    lnKb = (-8966.90 - 2890.53*sqrt(S) - 77.942*S + 1.728*S^1.5 - 0.0996*S^2)/T_K + 148.0248 + 137.1942*sqrt(S) + 1.62142*S + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T_K) + 0.053105*sqrt(S)*T_K;
    Kb = exp(lnKb) * 1e6;
    % Kw
    Kw = exp(148.96502 - 13847.26/T_K - 23.6521*log(T_K) + (118.67/T_K - 5.977 + 1.0495*log(T_K))*sqrt(S) - 0.01615*S) * 1e12;
    % calculate H+ abundance based on DIC and Alk:
    p5 = -1;
    p4 = -ALK - Kb - K1;
    p3 = DIC*K1 - ALK*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = DIC*(Kb*K1 + 2*K1*K2) - ALK*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*DIC*Kb*K1*K2 - ALK*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    r = roots([p5 p4 p3 p2 p1 p0]);
    H = max(real(r));
    % if isempty(r)
    %     H = 1e-8;
    % else
    %     H = max(real(r));
    % end
    pH = -log10(H) + 6;
    H2CO3 = DIC / (1 + K1/H + K1*K2/H^2);
    CO3 = DIC / (1 + H/K2 + H^2/(K1*K2));
end
```

## File: Run_RTM_1D.m
```matlab

clear all
tic
% ----------------------------- INPUT PARAMETERS ---------------------------
global v_burial Mineral_Mass z_sed KFe_HS Oxygen Sulfate
global k_sed k_O2 DSO4 DH2S DO2 DPO4 k_SO4 Kreox Iron_conc Bioturb Calcium DHCO3 HCO3init
global O2init SO4init HSinit C_organic rho poros RC Alpha_Bioirrig
global R_respi R_SRR Ksp_ca k_calcite DICinit R1_carb CO3_1 BE P_C_ratio Rviv1 R_FeS R_iron R_FeOx Fe_3_init
global v_burial_Fluid CO3_activity Calcium_activity NPP kFeS FeooH Feinit Iron_C R_HS_Ox kapatite P_apaeq R1_carb_disso R1_carb_form
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2 CH4init Pinitial DCH4 kFeOx KFEMonod Sulfide Rapat CaCO3 F_CaCO3 O2_root
global C_HS C_Fe n_power_CaCO31 n_power_CaCO32 k_calcite_dis1 n_power_CaCO33 k_calcite_dis2 CaCO3_init Temp_factor T_future
    Params = Params_Static();
    Config = Config_Baseline();
    Hydro  = Hydro_Preprocessor(Config, Params);
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
    KFeS = Params.KFeS;
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
    BE              = Config.BE;
    NPP             = Config.NPP;
    F_CaCO3         = Config.F_CaCO3;
    T_future        = Config.T_future;
    BE = Config.BE;
    NPP = Config.NPP;
    if Config.use_hydro_npp_multiplier
        NPP = NPP * Hydro.NPP_multiplier;
    end
    F_FeOx = Config.F_FeOx;
    F_CaCO3 = Config.F_CaCO3;
    T_future = Config.T_future;
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
[~, CO3_top, ~] = River_Carbonate(HCO3init, DICinit, T_future, 0.1, 1);
% CO3_top = Carb_CO3(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
CO3_1 = CO3_top*ones(1,n);
C_HS  = zeros(1,n); %intial value for sulfide
R_HS_Ox = zeros(1,n);
CaCO3 = 1E5.*zeros(1,n);
Rapat = zeros(1,n);
Rviv1 = zeros(1,n); %umol/Lsed/yr
R_FeS = ones(1,n); %umol/Lsed/yr
pH    = 7.*ones(1,n);
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
BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic
end
% ---------------------------- OXYGEN -------------------------------------
RC = Temp_factor.*k_sed.*C_organic.*rho.*((1-poros)./(12)); % molCorg/cm3/yr mineralization rate
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
% ---------------------------- OXYGEN PENTRATION DEPTH --------------------
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
    FeOxInit  = (F_FeOx.*36.5)./(v_burial(1) * rho * (1-poros(1)));
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
k_sed(1,1:num_OPD) = 10.^(-b_oxic*log10(age(1,1:num_OPD)) - a_oxic);  % Oxic
k_sed(1,(num_OPD+1):n) = 10.^(-b_anoxic*log10(age(1,(num_OPD+1):n)) - a_anoxic);  % Anoxic
% ------------------- Blue Carbon ------------------------------------
% DOC release
depth_rootzone = 10; % seagrass root length (cm)
z_root = 0:0.1:depth_rootzone;
mu_root = 4; % value for the center of rootzone in normal distribution
sigma_root = 6; % sigma for normal distribution of flux in the rootzone
DOC_root_2 = DOC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/year
DOC_root_22 = interp1(z_root,DOC_root_2,z_sed);
poros_root = interp1(z_sed,poros,z_sed);
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
K_converge = 1;
iteration = 1;
iteration_tolerance = 0.3;
count_loop = 1;
while abs(K_converge) > iteration_tolerance %for count_loop = 1:5
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
BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic
end
% ---------------------------- OXYGEN -------------------------------------
RC = Temp_factor.*k_sed.*C_organic.*rho.*((1-poros)./(12)) + RC_root; % molCorg/cm3/yr mineralization rate
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
R_respi = RC.* (Oxygen./(Oxygen+k_O2)).*1E9; %rate of aerobic respiration umol/l/year
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
Iron_C(iteration,:) = C_Fe;
for i=1:n
  if Iron_C(iteration,i) < 0
      Iron_C(iteration,i) = 0;
  end
end
F_diff_Fe(1,count_loop) = DH2S.*((C_Fe(1,2) - C_Fe(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr
% ------------------------ IRON(III) ---------------------------------------
Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
R_iron(count_loop,:) = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9; %rate of iron reduction umol/l/year
R_FeOx = (kFeOx.*C_Fe.*Oxygen);
R_FeOx_1(count_loop,:) = R_FeOx;
% Fe_3_init  = 365.*1E2.*(F_FeOx)./(v_burial(1));  %umol/l
Fe_3_init = 36.5.*(F_FeOx).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %umol/l
% KFEMonod = 2000;
% % Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[Fe_3_init 0]);
sol = bvp4c(@Fe3_ODE,@Fe3_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('Fe3：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_Fe_3 = y(1,:);
FeooH = C_Fe_3;
% FeooH = max(C_Fe_3, 0.2 .* FeOx);
% ------------------------ SULFATE ---------------------------------------
% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@SO4_ODE,@SO4_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);
if min(y) < 0
    fprintf('SO4 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);
C_SO4 = y(1,:);
Sulfate = C_SO4;
R_SRR = RC.*Inhib.* (Sulfate./(Sulfate+k_SO4)).*1E9; %rate of sulfate reduction umol/l/year
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
C_HS = y(1,:);
Sulfide(iteration,:) = C_HS;
for i=1:n
  if Sulfide(iteration,i) < 0
      Sulfide(iteration,i) = 0;
  end
end
HS_conc = Sulfide(iteration,:)./(1+((10.^(6-pH))./K_HS));
sigma_FeS_1 =  (Iron_C(iteration,:).*HS_conc)./((10.^(6-pH)).*KFeS);
delta_FeS  = (sigma_FeS_1 - 1);
for i=1:n
  if delta_FeS (1,i) > 0
      delta_FeS1(1,i) = 1;
  else
      delta_FeS1(1,i) = 0;
  end
end
R_FeS = kFeS.*C_Fe.*C_HS;
R_FeS_1(count_loop,:) = R_FeS;
R_FeS_store(iteration,:) = R_FeS;
for i=1:n
  if R_FeS(1,i) < 0
      R_FeS(1,i) = 0;
  end
end
R_HS_Ox = Sulfide(iteration,:).* Oxygen.*Kreox; %rate of sulfide reduction umol/l/year
F_diff_HS = DH2S.*((C_HS(1,2) - C_HS(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr
% ------------------------ METHANE ---------------------------------------
%XL calculating inhibition
Inhib_O2_meth  = (k_O2 ./ (Oxygen + k_O2));
Inhib_Fe_meth  = (KFEMonod ./ (FeooH + KFEMonod));
Inhib_SO4_meth = (k_SO4 ./ (Sulfate + k_SO4));
global Rate_Meth
Rate_Meth= 0.5 .* RC .* Inhib_O2_meth .* Inhib_Fe_meth .* Inhib_SO4_meth .* 1E9;
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
% ------------------------------- Coupled Carbonate -----------------------------------
CaCO3_init = 1E-4.*(F_CaCO3).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %gr/grDw
solinit_coupled = bvpinit(linspace(0, Lbottom, nmesh), [DICinit, 0, HCO3init, 0, CaCO3_init, 0]);
sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled);
% options_coupled = bvpset('NMax', 5000, 'RelTol', 1e-2);
% sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled, options_coupled);
y_coupled = deval(sol_coupled, x);
C_DIC  = y_coupled(1,:);
C_alka = y_coupled(3,:);
CaCO3  = y_coupled(5,:);
ALK    = C_alka;
F_diff_DIC = DHCO3 .* ((C_DIC(1,2) - C_DIC(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;
F_diff     = DHCO3 .* ((C_alka(1,2) - C_alka(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;
% -------------- pH and H2CO3 (2 for 6 calculation) -----------------------
for i=1:n
           [pH_1(1,i), CO3_1(1,i), C_H2CO3(1,i)] = River_Carbonate(ALK(1,i), C_DIC(1,i), T_future, 0.1, 1);
end
pH = pH_1;
sigma_carb = (Calcium_activity .* Calcium .* CO3_activity .* CO3_1) ./ Ksp_ca - 1;
unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho .* (1 - poros));
R_carb_form = (sigma_carb > 0) .* abs(sigma_carb).^n_power_CaCO31 .* k_calcite .* unit_conversion;
R_carb_disso = (-0.2 < sigma_carb & sigma_carb < 0) .* abs(sigma_carb).^n_power_CaCO32 .* k_calcite_dis1 .* CaCO3 ...
             + (sigma_carb <= -0.2) .* abs(sigma_carb).^n_power_CaCO33 .* k_calcite_dis2 .* CaCO3;
R1_carb = R_carb_form - R_carb_disso .* (1E3 * 1E6 * 1E-2 * rho .* (1 - poros));
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
alpha_converge(1,iteration) = F_diff;
% alpha_converge(isnan(alpha_converge))=[];
     if size(alpha_converge,2) > 2
       K_converge = (alpha_converge(1,iteration)-alpha_converge(1,iteration-1))/alpha_converge(1,iteration-1);
     end
       iteration = iteration + 1;
count_loop = count_loop + 1;
end  % iteration ends here
% Storing steady-state solutions as initial conditions for PDEs
Diff_fluxes = [F_diff_DIC F_diff F_diff_DIC./F_diff]';
R_ALK_DIC = F_diff./F_diff_DIC;
F_diff_CH4 = DCH4.*((CH4(1,2) - CH4(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr
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
% Mineral saturation indices
% subplot(m_plot,n_plot,16);
% plot(delta_viv1,z_sed,delta_apat2,z_sed,delta_FeS,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Mineral saturation (\Omega - 1)')
% legend('Vivianite','Apatite','FeS');
%
% box on
% grid on
%
% ax.LineWidth = 2;
%
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
% R_HSOX_integ = cumsum(R_HS_Ox.*dz_sed.*1E-3);
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
toc
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

function dydx = SO4_ODE(x,SO4)
global k_SO4 RC Oxygen k_O2 DSO4 v_burial_Fluid SO4init Alpha_Bioirrig z_sed poros
v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
Inh = (k_O2./(O2+k_O2));
RC1 = interp1(z_sed,RC,x);
% NR = + v_burial_f.* (SO4(2)) + 0.5.*RC1.*Inh.* (SO4(1)/(SO4(1)+k_SO4)) * 1E9 - (Alpha_Bioirrig_1.*(SO4init-SO4(1))); % umol/l/year
NR = + v_burial_f.* (SO4(2)/(fi*DSO4)) + 0.5.*RC1.*Inh.* (SO4(1)/(SO4(1)+k_SO4)) * 1E9 - (Alpha_Bioirrig_1.*(SO4init-SO4(1))); % umol/l/year
dydx = [ SO4(2) /fi/DSO4
           NR];
end
```

