function dydx = Fe3_ODE(x,Fe3)
global z_sed Bioturb RC R_FeOx 
global rho Oxygen KFEMonod k_O2 v_burial C_Fe kFeOx poros

RC_1 = interp1(z_sed,RC,x);
R_FeOx_1 = interp1(z_sed,R_FeOx,x);
poros_1 = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
C_Fe_1 = interp1(z_sed,C_Fe,x);
Inh = (k_O2./(O2+k_O2));
Db = interp1(z_sed,Bioturb,x);
v_burial_1 = interp1(z_sed,v_burial,x);


NR = - 4.*RC_1.*Inh.*1E9.*(Fe3(1)./(Fe3(1)+KFEMonod)).*(poros_1./(1-poros_1)).*1E-3.*(1./rho) + ...
       (kFeOx.*C_Fe_1.*O2).*(poros_1./(1-poros_1)).*1E-3.*(1./rho); % umol/g/year
 
dydx = [ NR / v_burial_1
         0];

end




