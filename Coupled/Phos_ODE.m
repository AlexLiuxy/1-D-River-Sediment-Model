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