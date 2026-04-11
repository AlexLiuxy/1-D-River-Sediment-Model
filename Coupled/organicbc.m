% function res = organicbc(C_orga,C_orgb)
% global NPP v_burial poros rho Bioturb BE
% 
% NPP1 = BE * NPP * 1E-4; %gram/cm2/year
% v_burial1 = v_burial(1,1);  %cm/year
% poros1 = poros(1,1); 
% A1 = rho * (1-poros1);
% Bioturb1 = Bioturb(1,1);
% 
% BC_1 = - Bioturb1 * A1 * C_orga(2) + A1 * v_burial1 * C_orga(1) - NPP1;
% 
% res = [ BC_1 
%         C_orgb(2)];
% 
% end

function res = organicbc(C_orga,C_orgb)
global Corg_top
res = [ C_orga(1) - Corg_top     % top: fixed solid-phase OM concentration
        C_orgb(2) ];             % bottom: zero gradient / zero diffusive flux
end