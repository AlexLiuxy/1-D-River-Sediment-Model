function dydx = HS_ODE(x,H2S)
% global k_SO4 RC Oxygen k_O2 v_burial_Fluid HSinit 
% global Kreox Alpha_Bioirrig z_sed poros Sulfate DH2S R_FeS C_Fe kFeS
global k_SO4 Oxygen v_burial_Fluid HSinit
global Kreox Alpha_Bioirrig z_sed poros Sulfate DH2S R_FeS C_Fe kFeS RC_after_Fe

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
SO4 = interp1(z_sed,Sulfate,x);
C_Fe_1 = interp1(z_sed,C_Fe,x);
% Inh = (k_O2./(O2+k_O2));
% R_FeS_1 = interp1(z_sed,R_FeS,x);
% RC1 = interp1(z_sed,RC,x);

RC_Fe_resid = interp1(z_sed, RC_after_Fe, x);
NR = + v_burial_f.* (H2S(2)/(fi*DH2S)) - 0.5 .* RC_Fe_resid .* (SO4/(SO4+k_SO4)) ...
     + (kFeS.*C_Fe_1.*H2S(1)) + (Kreox.*H2S(1).*O2) - (Alpha_Bioirrig_1.*(HSinit-H2S(1)));

% NR = + v_burial_f.* H2S(2) - 0.5.*RC1.*Inh.* (SO4/(SO4+k_SO4)).* 1E9 ...
%      + (kFeS.*C_Fe_1.*H2S(1)) + (Kreox.*H2S(1).*O2) - (Alpha_Bioirrig_1.*(HSinit-H2S(1))); % umol/l/year
% NR = + v_burial_f.* (H2S(2)/(fi*DH2S)) - 0.5.*RC1.*Inh.* (SO4/(SO4+k_SO4)).* 1E9 ...
%      + (kFeS.*C_Fe_1.*H2S(1)) + (Kreox.*H2S(1).*O2) - (Alpha_Bioirrig_1.*(HSinit-H2S(1))); % umol/l/year
dydx = [ H2S(2) /fi/DH2S
           NR];

end