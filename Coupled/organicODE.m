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

