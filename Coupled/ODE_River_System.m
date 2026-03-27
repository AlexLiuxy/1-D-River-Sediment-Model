function dYdx = ODE_River_System(x, Y, Config, Params)
    % calc physical properties (porosity and bioirrigation)
    % assuming exponential drop with depth
    phi = Config.porosbottom + (Config.porostop - Config.porosbottom) * exp(-x / Config.porosscale);
    alpha_bio = Config.Bioirrig_top * exp(-x / Config.Bioirrig_scale);
    v_burial = Config.vbottom; 
    
    % unpack variables from Y
    % odd = conc, even = flux
    O2  = Y(1);  dO2dx  = Y(2);
    SO4 = Y(3);  dSO4dx = Y(4);
    CH4 = Y(5);  dCH4dx = Y(6);
    DIC = Y(7);  dDICdx = Y(8);
    ALK = Y(9);  dALKdx = Y(10);
    Ca  = Y(11); dCadx  = Y(12);
    Fe  = Y(13); dFedx  = Y(14);
    HS  = Y(15); dHSdx  = Y(16);

    % force positive values to avoid kinetic blowups
    O2  = max(O2,  1e-9);
    SO4 = max(SO4, 1e-9);
    CH4 = max(CH4, 1e-9);
    DIC = max(DIC, 1e-9);
    Ca  = max(Ca,  1e-9);
    Fe  = max(Fe,  1e-9);
    HS  = max(HS,  1e-9);

    % ---------------------------------------------------------------------
    % CARBONATE SYSTEM & pH (copied from original code)
    % ---------------------------------------------------------------------
    % standard thermo constants for freshwater (can be moved to Params later)
    K1 = 1.18e-6;   % H2CO3 -> HCO3
    K2 = 4.36e-11;  % HCO3 -> CO3
    Kw = 1e-14;     % water
    Kb = 2.3e-9;    % borate
    B_T = 4e-4;     % total boron
    
    % polynomial roots for H+
    p5 = -1;
    p4 = -ALK - Kb - K1;
    p3 = DIC*K1 - ALK*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = DIC*(Kb*K1 + 2*K1*K2) - ALK*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*DIC*Kb*K1*K2 - ALK*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    
    r = roots([p5 p4 p3 p2 p1 p0]);
    H = max(real(r));
    % pH = -log10(H); % keep if needed for plotting later
    
    % carbonate speciation
    CO3 = DIC / (1 + H/K2 + H^2/(K1*K2));
    
    % calcite saturation state (Omega)
    Omega_ca = (Ca * CO3) / Params.Ksp_ca;
    
    % calcite precip/dissolution rate
    if Omega_ca > 1
        R_calcite = Params.k_calcite * (Omega_ca - 1)^2; % precipitation
    else
        R_calcite = -Params.k_calcite_dis1 * (1 - Omega_ca); % dissolution
    end

    % ---------------------------------------------------------------------
    % REACTION RATES
    % ---------------------------------------------------------------------
    % AOM (CH4 + SO4 -> HCO3 + HS + H2O)
    R_AOM = Params.k_AOM * CH4 * (SO4 / (SO4 + Params.K_CH4_SO4));
    
    % Aerobic Methanotrophy (CH4 + 2O2 -> CO2 + 2H2O)
    R_aerobic_CH4 = Params.k_aerobic_CH4 * CH4 * (O2 / (O2 + Params.K_CH4_O2));
    
    % Sulfide oxidation (HS + 2O2 -> SO4 + H+)
    R_HS_ox = Params.Kreox * HS * O2;
    
    % Fe oxidation (Fe2+ + 0.25O2 -> FeOOH)
    R_Fe_ox = Params.kFeOx * Fe * O2;
    
    % FeS precipitation (Fe2+ + HS- -> FeS + H+)
    R_FeS = Params.kFeS * Fe * HS;

    % ---------------------------------------------------------------------
    % ASSEMBLE ODEs
    % governing eq: D * d2C/dx2 - v * dC/dx + alpha*(C0 - C) + Reactions = 0
    % ---------------------------------------------------------------------
    dYdx = zeros(16, 1);
    
    % O2
    dYdx(1) = dO2dx;
    dYdx(2) = (v_burial * dO2dx - alpha_bio * (Config.O2init - O2) ...
               + 2*R_aerobic_CH4 + 2*R_HS_ox + 0.25*R_Fe_ox) / Params.DO2;
               
    % SO4
    dYdx(3) = dSO4dx;
    dYdx(4) = (v_burial * dSO4dx - alpha_bio * (Config.SO4init - SO4) ...
               + R_AOM - R_HS_ox) / Params.DSO4;
               
    % CH4
    dYdx(5) = dCH4dx;
    dYdx(6) = (v_burial * dCH4dx - alpha_bio * (Config.CH4init - CH4) ...
               + R_aerobic_CH4 + R_AOM) / Params.DCH4; % add methanogenesis later if needed
               
    % DIC
    dYdx(7) = dDICdx;
    dYdx(8) = (v_burial * dDICdx - alpha_bio * (Config.DICinit - DIC) ...
               - R_aerobic_CH4 - R_AOM + R_calcite) / Params.DHCO3; 
               
    % ALK 
    dYdx(9) = dALKdx;
    dYdx(10) = (v_burial * dALKdx - alpha_bio * (Config.HCO3init - ALK) ...
                - R_AOM + 2*R_calcite + 2*R_HS_ox + R_FeS) / Params.DHCO3;
                
    % Ca
    dYdx(11) = dCadx;
    dYdx(12) = (v_burial * dCadx - alpha_bio * (Config.Calcium - Ca) ...
                + R_calcite) / Params.DHCO3; % using HCO3 diff coeff for Ca as approx
                
    % Fe
    dYdx(13) = dFedx;
    dYdx(14) = (v_burial * dFedx - alpha_bio * (Config.Feinit - Fe) ...
                + R_Fe_ox + R_FeS) / Params.DHCO3; % add Fe reduction later
                
    % HS
    dYdx(15) = dHSdx;
    dYdx(16) = (v_burial * dHSdx - alpha_bio * (Config.HSinit - HS) ...
                - R_AOM + R_HS_ox + R_FeS) / Params.DH2S;
end
