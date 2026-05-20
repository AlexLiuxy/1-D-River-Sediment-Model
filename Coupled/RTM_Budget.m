function Budget = RTM_Budget(Result)
% RTM_BUDGET
% PDE-consistent Fe, CH4, and FeOOH budget diagnostics.
%
% Solute budgets use phi-weighted reaction integrals because PDE storage is:
%   d(phi*C)/dt = transport + phi*reaction + phi*exchange
%
% Solid budgets use A-weighted reaction integrals:
%   A = rho*(1-phi)

Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;
P = Result.Params;
F = Result.Forcing;

z = Grid.z(:);
dz = Grid.dz;
phi = Grid.poros(:);
A = Grid.rho .* max(1 - phi, 1e-6);

Budget = struct();

% =========================
% CH4 budget
% =========================

CH4 = S.CH4(:);

J_CH4_top_down = top_solute_flux(CH4, F.CH4_top, phi, P.DCH4, Grid.v_fluid, dz);
J_CH4_top_up = max(0, -J_CH4_top_down);

I_CH4_irrig_sink = trapz(z, phi .* Grid.Alpha_exchange(:) .* max(CH4 - F.CH4_top, 0)) .* 1e-3;

Budget.CH4.I_Meth = trapz(z, D.R_Meth(:)) .* 1e-3;
Budget.CH4.I_AOM  = trapz(z, D.R_AOM(:)) .* 1e-3;
Budget.CH4.I_Ox   = trapz(z, D.R_CH4Ox(:)) .* 1e-3;

Budget.CH4.I_IrrigSink = I_CH4_irrig_sink;
Budget.CH4.J_TopUp = J_CH4_top_up;

Budget.CH4.Source = Budget.CH4.I_Meth;
Budget.CH4.Sink = Budget.CH4.I_AOM + Budget.CH4.I_Ox + ...
                  Budget.CH4.I_IrrigSink + Budget.CH4.J_TopUp;

Budget.CH4.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dCH4dt = (S.CH4(:) - S_prev.CH4(:)) ./ max(dt, 1e-12);
    Budget.CH4.Storage = trapz(z, phi .* dCH4dt) .* 1e-3;
end

Budget.CH4.Residual = Budget.CH4.Source - Budget.CH4.Sink - nan0(Budget.CH4.Storage);
Budget.CH4.SinkSourceRatio = Budget.CH4.Sink ./ max(Budget.CH4.Source, 1e-12);

% =========================
% Fe2 budget
% =========================

Fe2 = S.Fe2(:);

J_Fe2_top_down = top_solute_flux(Fe2, F.Fe_top, phi, P.DH2S, Grid.v_fluid, dz);
J_Fe2_top_up = max(0, -J_Fe2_top_down);

I_Fe2_irrig_sink = trapz(z, phi .* Grid.Alpha_exchange(:) .* max(Fe2 - F.Fe_top, 0)) .* 1e-3;

Budget.Fe2.I_FeRed = trapz(z, D.R_FeRed(:)) .* 1e-3;
Budget.Fe2.I_FeOx  = trapz(z, D.R_FeOx(:)) .* 1e-3;
Budget.Fe2.I_FeS   = trapz(z, D.R_FeS(:)) .* 1e-3;

Budget.Fe2.I_IrrigSink = I_Fe2_irrig_sink;
Budget.Fe2.J_TopUp = J_Fe2_top_up;

Budget.Fe2.Source = Budget.Fe2.I_FeRed;
Budget.Fe2.Sink = Budget.Fe2.I_FeOx + Budget.Fe2.I_FeS + ...
                  Budget.Fe2.I_IrrigSink + Budget.Fe2.J_TopUp;

Budget.Fe2.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dFe2dt = (S.Fe2(:) - S_prev.Fe2(:)) ./ max(dt, 1e-12);
    Budget.Fe2.Storage = trapz(z, phi .* dFe2dt) .* 1e-3;
end

Budget.Fe2.Residual = Budget.Fe2.Source - Budget.Fe2.Sink - nan0(Budget.Fe2.Storage);
Budget.Fe2.SinkSourceRatio = Budget.Fe2.Sink ./ max(Budget.Fe2.Source, 1e-12);

% =========================
% FeOOH solid budget
% =========================

FeOOH = S.FeOOH(:);

F_FeOOH_top = F.F_FeOx .* 36.5;
if isfield(P, 'Fe_inventory_factor')
    F_FeOOH_top = P.Fe_inventory_factor .* F_FeOOH_top;
end

F_FeOOH_bottom = A(end) .* Grid.v_solid(end) .* FeOOH(end);

% FeOOH reaction terms are equivalent to phi-weighted Fe solute rates.
Budget.FeOOH.TopInput = F_FeOOH_top;
Budget.FeOOH.BottomBurial = F_FeOOH_bottom;
Budget.FeOOH.I_FeRed = trapz(z, D.R_FeRed(:)) .* 1e-3;
Budget.FeOOH.I_FeOx  = trapz(z, D.R_FeOx(:)) .* 1e-3;

Budget.FeOOH.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dFeOOHdt = (S.FeOOH(:) - S_prev.FeOOH(:)) ./ max(dt, 1e-12);
    Budget.FeOOH.Storage = trapz(z, A .* dFeOOHdt);
end

Budget.FeOOH.Net = Budget.FeOOH.TopInput + Budget.FeOOH.I_FeOx ...
                   - Budget.FeOOH.I_FeRed - Budget.FeOOH.BottomBurial ...
                   - nan0(Budget.FeOOH.Storage);

% =========================
% Print
% =========================

fprintf('\n================ RTM Budget ================\n');

fprintf('\n--- CH4 budget, umol/cm2/yr ---\n');
fprintf('Meth source:          %.3f\n', Budget.CH4.I_Meth);
fprintf('AOM sink:             %.3f\n', Budget.CH4.I_AOM);
fprintf('O2 oxidation sink:    %.3f\n', Budget.CH4.I_Ox);
fprintf('Irrigation sink:      %.3f\n', Budget.CH4.I_IrrigSink);
fprintf('Top diffusive efflux: %.3f\n', Budget.CH4.J_TopUp);
fprintf('Storage:              %.3f\n', Budget.CH4.Storage);
fprintf('Sink/source ratio:    %.3f\n', Budget.CH4.SinkSourceRatio);
fprintf('Residual:             %.3f\n', Budget.CH4.Residual);

fprintf('\n--- Fe2 budget, umol/cm2/yr ---\n');
fprintf('Fe reduction source:  %.3f\n', Budget.Fe2.I_FeRed);
fprintf('Fe oxidation sink:    %.3f\n', Budget.Fe2.I_FeOx);
fprintf('FeS sink:             %.3f\n', Budget.Fe2.I_FeS);
fprintf('Irrigation sink:      %.3f\n', Budget.Fe2.I_IrrigSink);
fprintf('Top diffusive efflux: %.3f\n', Budget.Fe2.J_TopUp);
fprintf('Storage:              %.3f\n', Budget.Fe2.Storage);
fprintf('Sink/source ratio:    %.3f\n', Budget.Fe2.SinkSourceRatio);
fprintf('Residual:             %.3f\n', Budget.Fe2.Residual);

fprintf('\n--- FeOOH budget, umol/cm2/yr ---\n');
fprintf('Top FeOOH input:      %.3f\n', Budget.FeOOH.TopInput);
fprintf('FeOOH regeneration:   %.3f\n', Budget.FeOOH.I_FeOx);
fprintf('FeOOH reduction sink: %.3f\n', Budget.FeOOH.I_FeRed);
fprintf('Bottom burial loss:   %.3f\n', Budget.FeOOH.BottomBurial);
fprintf('Storage:              %.3f\n', Budget.FeOOH.Storage);
fprintf('Residual:             %.3f\n', Budget.FeOOH.Net);

fprintf('============================================\n\n');

end

function J_top_down = top_solute_flux(C, C_top, phi, D, v, dz)
% Positive downward. Convert to umol/cm2/yr.

C = C(:);
phi = phi(:);
v = v(:);

if v(1) >= 0
    C_adv = C_top;
else
    C_adv = C(1);
end

F_top = -phi(1) .* D .* (C(1) - C_top) ./ (0.5 .* dz) ...
        + phi(1) .* v(1) .* C_adv;

J_top_down = F_top .* 1e-3;
end

function x = nan0(x)
if isnan(x)
    x = 0;
end
end