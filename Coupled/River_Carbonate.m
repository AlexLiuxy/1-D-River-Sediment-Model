function [pH, CO3, H2CO3] = River_Carbonate(ALK, DIC, T, S, P)
% A unified carbonate thermodynamic solver for USGS river modifications
% Input: ALK (umol/kg), DIC (umol/kg), T ( °C), S (salinity ppt), P (pressure at depth atm)


    % If no input, use default value
    if nargin < 3
        T = 20;   
        S = 0.1; 
        P = 1;    
    end

    T_K = T + 273.15; 
    B_T = 400 * (S / 35); % total boron in river [umol/kg]
    
    % Calculate K
    K0 = exp((9345.17/T_K) - 60.2409 + 23.3585*log(T_K/100) + S*(0.023517 - 0.00023656*T_K + 0.00000047036*T_K^2));
    
    
    RGAS = 8.314510; R = 83.131;
    
    % K1
    lnK1 = 2.83655 - 2307.1266/T_K - 1.5529413*log(T_K) - (0.20760841 + 4.0484/T_K)*sqrt(S) + 0.08468345*S - 0.00654208*S^1.5 + log(1 - 0.001005*S);
    K1 = exp(lnK1) * 1e6; 
    
    % K2
    lnK2 = -9.226508 - 3351.6106/T_K - 0.2005743*log(T_K) + (-0.106901773 - 23.9722/T_K)*sqrt(S) + 0.1130822*S - 0.00846934*S^1.5 + log(1 - 0.001005*S);
    K2 = exp(lnK2) * 1e6;
    
    % Kb 
    lnKb = (-8966.90 - 2890.53*sqrt(S) - 77.942*S + 1.728*S^1.5 - 0.0996*S^2)/T_K + 148.0248 + 137.1942*sqrt(S) + 1.62142*S + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T_K) + 0.053105*sqrt(S)*T_K;
    Kb = exp(lnKb) * 1e6;
    
    % Kw 
    Kw = exp(148.96502 - 13847.26/T_K - 23.6521*log(T_K) + (118.67/T_K - 5.977 + 1.0495*log(T_K))*sqrt(S) - 0.01615*S) * 1e12;

    % calculate H+ abundance based on DIC and Alk:
    p5 = -1;
    p4 = -ALK - Kb - K1;
    p3 = DIC*K1 - ALK*(Kb+K1) + Kb*B_T + Kw - Kb*K1 - K1*K2;
    p2 = DIC*(Kb*K1 + 2*K1*K2) - ALK*(Kb*K1 + K1*K2) + Kb*B_T*K1 + Kw*Kb + Kw*K1 - Kb*K1*K2;
    p1 = 2*DIC*Kb*K1*K2 - ALK*Kb*K1*K2 + Kb*B_T*K1*K2 + Kw*Kb*K1 + Kw*K1*K2;
    p0 = Kw*Kb*K1*K2;
    
    r = roots([p5 p4 p3 p2 p1 p0]);
   
    H = max(real(r));

    % if isempty(r)
    %     H = 1e-8; 
    % else
    %     H = max(real(r)); 
    % end
    
    
    pH = -log10(H) + 6;
    H2CO3 = DIC / (1 + K1/H + K1*K2/H^2);
    CO3 = DIC / (1 + H/K2 + H^2/(K1*K2));

end