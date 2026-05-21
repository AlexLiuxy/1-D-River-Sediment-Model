function [Sulfate, R_SRR, R_AOM_actual, SO4_diag] = Solve_SO4_FV_Front()
% Positivity-preserving finite-volume SO4 solver with depletion front behavior.
%
% This replaces the ordinary SO4 bvp4c solve.
% It keeps SO4 >= 0 and makes SRR/AOM vanish in sulfate-depleted cells.

global z_sed poros DSO4 SO4init Alpha_Bioirrig
global k_SO4 K_CH4_SO4 RC_after_Fe R_AOM_pot

z = z_sed(:);
n = numel(z);

if n < 3
    error('Solve_SO4_FV_Front requires at least 3 grid cells.');
end

dz_all = diff(z);
if max(abs(dz_all - dz_all(1))) > 1e-10
    error('Solve_SO4_FV_Front currently assumes uniform z_sed.');
end
dz = dz_all(1);

phi = poros(:);
Alpha = Alpha_Bioirrig(:);

Knode = max(phi .* DSO4, 1e-12);
Kface = 0.5 .* (Knode(1:end-1) + Knode(2:end));

R_SRR_pot = max(0, 0.5 .* RC_after_Fe(:));  % umol SO4/L/yr
R_AOM_pot_col = max(0, R_AOM_pot(:));       % umol SO4/L/yr potential

S_floor = 0;
S_cut = 1e-6;
tol = 1e-8;
max_iter = 200;
relax = 0.7;

% Initial guess: positive, decreasing sulfate.
S = SO4init .* exp(-z ./ max(0.25*z(end), dz));
S(1) = SO4init;

for it = 1:max_iter
    S_old = max(S, S_cut);

    % Linearized Monod sinks:
    % R ~= lambda * S
    lambda_srr = R_SRR_pot ./ max(S_old + k_SO4, 1e-12);
    lambda_aom = R_AOM_pot_col ./ max(S_old + K_CH4_SO4, 1e-12);
    lambda = lambda_srr + lambda_aom;

    % Depleted cells remain non-reactive.
    depleted_old = S <= S_cut;
    lambda(depleted_old) = 0;

    A = spalloc(n, n, 3*n);
    b = zeros(n,1);

    % Top Dirichlet: bottom-water sulfate.
    A(1,1) = 1;
    b(1) = SO4init;

    for i = 2:n
        Km = Kface(i-1);

        if i < n
            Kp = Kface(i);
        else
            Kp = 0;  % bottom no-flux
        end

        % Equation:
        % -div(K grad S) + (Alpha + lambda)*S = Alpha*SO4init
        A(i,i-1) = -Km / dz^2;
        A(i,i)   = (Km + Kp) / dz^2 + Alpha(i) + lambda(i);

        if i < n
            A(i,i+1) = -Kp / dz^2;
        end

        b(i) = Alpha(i) .* SO4init;
    end

    S_lin = A \ b;
    S_new = max(S_floor, relax .* S_lin + (1 - relax) .* S);
    S_new(1) = SO4init;

    err = norm(S_new - S, inf) ./ max(SO4init, 1);
    S = S_new;

    if err < tol
        break
    end
end

% Active-set depletion projection.
S(S <= S_cut) = 0;
S(1) = SO4init;

% Actual rates from final nonnegative sulfate.
f_SRR = S ./ max(S + k_SO4, 1e-12);
f_AOM = S ./ max(S + K_CH4_SO4, 1e-12);

f_SRR(S <= S_cut) = 0;
f_AOM(S <= S_cut) = 0;

R_SRR_col = R_SRR_pot .* f_SRR;
R_AOM_col = R_AOM_pot_col .* f_AOM;

% Top node is a boundary value, not a reactive volume.
R_SRR_col(1) = 0;
R_AOM_col(1) = 0;

Sulfate = S(:).';
R_SRR = R_SRR_col(:).';
R_AOM_actual = R_AOM_col(:).';

idx_front = find(Sulfate <= 1, 1, 'first');
if isempty(idx_front)
    front_depth = z_sed(end);
else
    front_depth = z_sed(idx_front);
end

SO4_diag.iter = it;
SO4_diag.min_SO4 = min(Sulfate);
SO4_diag.front_depth = front_depth;
SO4_diag.I_SRR = trapz(z_sed, R_SRR) .* 1e-3;
SO4_diag.I_AOM = trapz(z_sed, R_AOM_actual) .* 1e-3;
SO4_diag.I_demand_pot = trapz(z_sed, R_SRR_pot(:).') .* 1e-3;
SO4_diag.I_irrig_source = trapz(z_sed, Alpha_Bioirrig .* (SO4init - Sulfate)) .* 1e-3;
SO4_diag.J_top_down = -DSO4 .* poros(1) .* ((Sulfate(2) - Sulfate(1)) ./ dz) .* 1e-3;
end