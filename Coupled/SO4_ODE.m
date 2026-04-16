function dydx = SO4_ODE(x,SO4)
global k_SO4 RC_after_Fe DSO4 v_burial_Fluid SO4init Alpha_Bioirrig z_sed poros R_AOM_lag

v_burial_f = interp1(z_sed, v_burial_Fluid, x);
Alpha_Bioirrig_1 = interp1(z_sed, Alpha_Bioirrig, x);
fi = interp1(z_sed, poros, x);
RC1 = interp1(z_sed, RC_after_Fe, x);
R_AOM_1 = interp1(z_sed, R_AOM_lag, x);

NR = + v_burial_f .* (SO4(2)/(fi*DSO4)) ...
     + 0.5 .* RC1 .* (SO4(1)/(SO4(1)+k_SO4)) ...
     + R_AOM_1 ...
     - (Alpha_Bioirrig_1 .* (SO4init-SO4(1)));

dydx = [ SO4(2) /fi/DSO4
         NR ];
end