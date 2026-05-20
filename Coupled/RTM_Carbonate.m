function Carb = RTM_Carbonate(ALK, DIC, T, S)
% RTM_CARBONATE
% Vectorized carbonate speciation solver for the transient PDE RHS.
%
% Uses fixed-iteration vectorized bisection in log10([H+] in uM).
% This avoids fzero and roots inside ode15s RHS calls.
%
% Inputs:
%   ALK, DIC : uM, vectors or scalars
%   T        : deg C
%   S        : salinity
%
% Outputs in Carb:
%   pH, CO3, H2CO3, HCO3, BOH4, OH

ALK = max(real(ALK(:)), 1e-12);
DIC = max(real(DIC(:)), 1e-12);

if nargin < 3 || isempty(T)
    T = 25;
end
if nargin < 4 || isempty(S)
    S = 0.1;
end

T_K = T + 273.15;

% River/low-salinity boron approximation inherited from current code.
B_T = 400 .* (S ./ 35);

lnK1 = 2.83655 - 2307.1266 ./ T_K - 1.5529413 .* log(T_K) ...
     - (0.20760841 + 4.0484 ./ T_K) .* sqrt(S) ...
     + 0.08468345 .* S - 0.00654208 .* S.^1.5 ...
     + log(1 - 0.001005 .* S);
K1 = exp(lnK1) .* 1e6;

lnK2 = -9.226508 - 3351.6106 ./ T_K - 0.2005743 .* log(T_K) ...
     + (-0.106901773 - 23.9722 ./ T_K) .* sqrt(S) ...
     + 0.1130822 .* S - 0.00846934 .* S.^1.5 ...
     + log(1 - 0.001005 .* S);
K2 = exp(lnK2) .* 1e6;

lnKb = (-8966.90 - 2890.53 .* sqrt(S) - 77.942 .* S ...
      + 1.728 .* S.^1.5 - 0.0996 .* S.^2) ./ T_K ...
     + 148.0248 + 137.1942 .* sqrt(S) + 1.62142 .* S ...
     + (-24.4344 - 25.085 .* sqrt(S) - 0.2474 .* S) .* log(T_K) ...
     + 0.053105 .* sqrt(S) .* T_K;
Kb = exp(lnKb) .* 1e6;

Kw = exp(148.96502 - 13847.26 ./ T_K - 23.6521 .* log(T_K) ...
   + (118.67 ./ T_K - 5.977 + 1.0495 .* log(T_K)) .* sqrt(S) ...
   - 0.01615 .* S) .* 1e12;

% Broadcast constants to vector size.
K1 = K1 .* ones(size(ALK));
K2 = K2 .* ones(size(ALK));
Kb = Kb .* ones(size(ALK));
Kw = Kw .* ones(size(ALK));
B_T = B_T .* ones(size(ALK));

% Bracket log10(H_uM).
% Current River_Carbonate used -7 to 2; keep same safe range.
lo = -7 .* ones(size(ALK));
hi =  2 .* ones(size(ALK));

flo = alk_residual(10.^lo, ALK, DIC, B_T, K1, K2, Kb, Kw);
fhi = alk_residual(10.^hi, ALK, DIC, B_T, K1, K2, Kb, Kw);

% If a rare point is not bracketed, keep the nearest endpoint.
bad = flo .* fhi > 0;

% Fixed iteration count is predictable and ode15s-friendly.
for iter = 1:50
    mid = 0.5 .* (lo + hi);
    fmid = alk_residual(10.^mid, ALK, DIC, B_T, K1, K2, Kb, Kw);

    % Residual generally decreases with increasing H.
    go_right = fmid > 0;
    lo(go_right) = mid(go_right);
    hi(~go_right) = mid(~go_right);
end

logH = 0.5 .* (lo + hi);

% Fallback for non-bracketed cells.
if any(bad)
    use_lo = abs(flo) < abs(fhi);
    logH(bad & use_lo) = lo(bad & use_lo);
    logH(bad & ~use_lo) = hi(bad & ~use_lo);
end

H = max(10.^logH, 1e-12);

denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
denom = max(denom, 1e-30);

H2CO3 = DIC ./ denom;
HCO3  = DIC .* (K1 ./ H) ./ denom;
CO3   = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
BOH4  = B_T .* (Kb ./ H) ./ (1 + Kb ./ H);
OH    = Kw ./ H;

Carb = struct();
Carb.pH    = 6 - log10(H);
Carb.H2CO3 = max(real(H2CO3), 1e-12);
Carb.HCO3  = max(real(HCO3),  1e-12);
Carb.CO3   = max(real(CO3),   1e-12);
Carb.BOH4  = max(real(BOH4),  0);
Carb.OH    = max(real(OH),    0);

end

function res = alk_residual(H, ALK, DIC, B_T, K1, K2, Kb, Kw)
    H = max(H, 1e-12);

    denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
    denom = max(denom, 1e-30);

    HCO3 = DIC .* (K1 ./ H) ./ denom;
    CO3  = DIC .* (K1 .* K2 ./ H.^2) ./ denom;
    BOH4 = B_T .* (Kb ./ H) ./ (1 + Kb ./ H);
    OH   = Kw ./ H;

    TA_calc = HCO3 + 2 .* CO3 + BOH4 + OH - H;
    res = TA_calc - ALK;
end