function [pH, CO3, H2CO3] = River_Carbonate(ALK, DIC, T, S, P)

    if nargin < 3 || isempty(T), T = 20; end
    if nargin < 4 || isempty(S), S = 0.1; end
    if nargin < 5 || isempty(P), P = 1; end %#ok<NASGU>

    ALK = max(real(ALK), 1e-12);
    DIC = max(real(DIC), 1e-12);

    T_K = T + 273.15;
    B_T = 400 * (S / 35);   % umol/kg

    lnK1 = 2.83655 - 2307.1266 / T_K - 1.5529413 * log(T_K) ...
         - (0.20760841 + 4.0484 / T_K) * sqrt(S) ...
         + 0.08468345 * S - 0.00654208 * S^(1.5) ...
         + log(1 - 0.001005 * S);
    K1 = exp(lnK1) * 1e6;

    lnK2 = -9.226508 - 3351.6106 / T_K - 0.2005743 * log(T_K) ...
         + (-0.106901773 - 23.9722 / T_K) * sqrt(S) ...
         + 0.1130822 * S - 0.00846934 * S^(1.5) ...
         + log(1 - 0.001005 * S);
    K2 = exp(lnK2) * 1e6;

    lnKb = (-8966.90 - 2890.53 * sqrt(S) - 77.942 * S + 1.728 * S^(1.5) - 0.0996 * S^2) / T_K ...
         + 148.0248 + 137.1942 * sqrt(S) + 1.62142 * S ...
         + (-24.4344 - 25.085 * sqrt(S) - 0.2474 * S) * log(T_K) ...
         + 0.053105 * sqrt(S) * T_K;
    Kb = exp(lnKb) * 1e6;

    Kw = exp(148.96502 - 13847.26 / T_K - 23.6521 * log(T_K) ...
       + (118.67 / T_K - 5.977 + 1.0495 * log(T_K)) * sqrt(S) - 0.01615 * S) * 1e12;

    f = @(logH) alk_balance_residual(10.^logH, ALK, DIC, B_T, K1, K2, Kb, Kw);

    % lighter bracket search
    logH_grid = linspace(-7, 2, 81);
    f_grid = zeros(size(logH_grid));
    for i = 1:numel(logH_grid)
        f_grid(i) = f(logH_grid(i));
    end

    idx = find(f_grid(1:end-1) .* f_grid(2:end) <= 0, 1, 'first');

    if ~isempty(idx)
        logH = fzero(f, [logH_grid(idx), logH_grid(idx+1)]);
    else
        [~, imin] = min(abs(f_grid));
        logH = logH_grid(imin);
    end

    H = max(10.^logH, 1e-12);
    pH = 6 - log10(H);

    denom = 1 + K1 ./ H + K1 .* K2 ./ H.^2;
    denom = max(denom, 1e-30);

    H2CO3 = DIC ./ denom;
    CO3   = DIC .* (K1 .* K2 ./ H.^2) ./ denom;

    H2CO3 = max(real(H2CO3), 1e-12);
    CO3   = max(real(CO3),   1e-12);
end

function res = alk_balance_residual(H, ALK, DIC, B_T, K1, K2, Kb, Kw)
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