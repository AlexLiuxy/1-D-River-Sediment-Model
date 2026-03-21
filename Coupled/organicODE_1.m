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