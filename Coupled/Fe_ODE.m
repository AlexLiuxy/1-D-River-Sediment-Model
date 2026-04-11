function dydx = Fe_ODE(x,Fe2)
% global RC Oxygen k_O2 DH2S v_burial_Fluid Alpha_Bioirrig z_sed poros Feinit
% global FeooH kFeOx KFEMonod R_FeS kFeS C_HS
global Oxygen DH2S v_burial_Fluid Alpha_Bioirrig z_sed poros Feinit
global kFeOx R_FeRed kFeS C_HS

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);
fi = interp1(z_sed,poros,x);
O2 = interp1(z_sed,Oxygen,x);
C_HS_1 = interp1(z_sed,C_HS,x);
% FeOx = interp1(z_sed,FeooH,x);

% Inh = (k_O2./(O2+k_O2));
% % R_FeS_1 = interp1(z_sed,R_FeS,x);
% 
% RC1 = interp1(z_sed,RC,x);
% 
% 
% % NR = + v_burial_f.* Fe2(2) - 4.*RC1.*Inh.* 1E9.* (FeOx./(FeOx+KFEMonod)) - (Alpha_Bioirrig_1.*(Feinit-Fe2(1)))...
% %      + (kFeS.*Fe2(1).*C_HS_1) + (kFeOx.*Fe2(1).*O2); % umol/l/year 
% 
% NR = + v_burial_f.* (Fe2(2)/(fi*DH2S)) - 4.*RC1.*Inh.* 1E9.* (FeOx./(FeOx+KFEMonod)) - (Alpha_Bioirrig_1.*(Feinit-Fe2(1)))...
%      + (kFeS.*Fe2(1).*C_HS_1) + (kFeOx.*Fe2(1).*O2); % umol/l/year 

R_FeRed_1 = interp1(z_sed, R_FeRed, x);
NR = + v_burial_f.* (Fe2(2)/(fi*DH2S)) - R_FeRed_1 ...
     - (Alpha_Bioirrig_1.*(Feinit-Fe2(1))) + (kFeS.*Fe2(1).*C_HS_1) + (kFeOx.*Fe2(1).*O2);

dydx = [ Fe2(2) /fi/DH2S
           NR];

end