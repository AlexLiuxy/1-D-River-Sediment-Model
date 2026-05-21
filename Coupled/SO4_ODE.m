% function dydx = SO4_ODE(x,SO4)
% global k_SO4 RC_after_Fe DSO4 v_burial_Fluid SO4init Alpha_Bioirrig z_sed poros R_AOM_lag
% 
% v_burial_f = interp1(z_sed, v_burial_Fluid, x);
% Alpha_Bioirrig_1 = interp1(z_sed, Alpha_Bioirrig, x);
% fi = interp1(z_sed, poros, x);
% RC1 = interp1(z_sed, RC_after_Fe, x);
% R_AOM_1 = interp1(z_sed, R_AOM_lag, x);
% 
% NR = + v_burial_f .* (SO4(2)/(fi*DSO4)) ...
%      + 0.5 .* RC1 .* (SO4(1)/(SO4(1)+k_SO4)) ...
%      + R_AOM_1 ...
%      - (Alpha_Bioirrig_1 .* (SO4init-SO4(1)));
% 
% dydx = [ SO4(2) /fi/DSO4
%          NR ];
% end

function dydx = SO4_ODE(x, SO4)
global k_SO4 K_CH4_SO4 RC_after_Fe DSO4 v_burial_Fluid SO4init Alpha_Bioirrig z_sed poros
global R_AOM_pot

v_burial_f = interp1(z_sed, v_burial_Fluid, x, 'linear', 'extrap');
Alpha      = interp1(z_sed, Alpha_Bioirrig, x, 'linear', 'extrap');
fi         = interp1(z_sed, poros, x, 'linear', 'extrap');
RC1        = interp1(z_sed, RC_after_Fe, x, 'linear', 'extrap');

% potential AOM from previous CH4 profile
R_AOM_pot_1 = interp1(z_sed, R_AOM_pot, x, 'linear', 'extrap');

% True nonnegative sulfate for reaction limitation.
% Important: no softplus. If SO4 <= 0, reactions using SO4 shut down.
SO4_pos = max(real(SO4(1)), 0);

f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);

R_SRR_local = 0.5 .* RC1 .* f_SRR;
R_AOM_local = R_AOM_pot_1 .* f_AOM;

NR = + v_burial_f .* (SO4(2) ./ (fi .* DSO4)) ...
     + R_SRR_local ...
     + R_AOM_local ...
     - Alpha .* (SO4init - SO4(1));

dydx = [ SO4(2) ./ (fi .* DSO4)
         NR ];
end