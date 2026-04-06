function res = Coupled_Carbonate_bc(Ya, Yb)
global DICinit HCO3init CaCO3_init

  
  res = [ Ya(1) - DICinit;       % 1: DIC top
          Yb(2);                 % 2: DIC bottom flux = 0
          Ya(3) - HCO3init;      % 3: ALK top
          Yb(4);                 % 4: ALK bottom flux = 0
          Ya(5) - CaCO3_init;    % 5: CaCO3 top
          Yb(6) ];               % 6: CaCO3 bottom flux = 0
end

% function res = Coupled_Carbonate_bc(Ya, Yb)
% global DICinit HCO3init Calcium CaCO3_init
% 
% res = [ Ya(1) - DICinit;    % DIC top
%         Yb(2);              % DIC bottom dissolved flux = 0
%         Ya(3) - HCO3init;   % ALK top
%         Yb(4);              % ALK bottom dissolved flux = 0
%         Ya(5) - Calcium;    % Ca top
%         Yb(6);              % Ca bottom dissolved flux = 0
%         Ya(7) - CaCO3_init; % CaCO3 top solid concentration
%         Yb(8) ];            % CaCO3 bottom mixing flux = 0
% end