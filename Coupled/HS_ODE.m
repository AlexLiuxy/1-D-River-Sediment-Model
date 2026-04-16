function dydx = HS_ODE(x,H2S)
global k_SO4 RC_after_Fe Oxygen v_burial_Fluid HSinit
global Kreox Alpha_Bioirrig z_sed poros Sulfate DH2S C_Fe kFeS pH K_HS R_AOM_lag

v_burial_f = interp1(z_sed, v_burial_Fluid, x);
Alpha_Bioirrig_1 = interp1(z_sed, Alpha_Bioirrig, x);
fi = interp1(z_sed, poros, x);
O2 = interp1(z_sed, Oxygen, x);
SO4 = interp1(z_sed, Sulfate, x);
C_Fe_1 = interp1(z_sed, C_Fe, x);
RC1 = interp1(z_sed, RC_after_Fe, x);
R_AOM_1 = interp1(z_sed, R_AOM_lag, x);

pH_x = interp1(z_sed, pH, x);
HS_free = H2S(1) ./ (1 + ((10.^(6 - pH_x)) ./ K_HS));

NR = + v_burial_f .* (H2S(2)/(fi*DH2S)) ...
     - 0.5 .* RC1 .* (SO4/(SO4+k_SO4)) ...
     - R_AOM_1 ...
     + (kFeS .* C_Fe_1 .* HS_free) ...
     + (Kreox .* HS_free .* O2) ...
     - (Alpha_Bioirrig_1 .* (HSinit-H2S(1)));

dydx = [ H2S(2) /fi/DH2S
         NR ];
end