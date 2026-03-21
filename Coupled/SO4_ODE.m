
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