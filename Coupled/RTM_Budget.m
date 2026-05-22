function Budget = RTM_Budget(Result)
% RTM_BUDGET
% PDE-consistent Fe, CH4, and FeOOH budget diagnostics.
%
% Solute budgets use phi-weighted reaction integrals because PDE storage is:
%   d(phi*C)/dt = transport + reaction + phi*exchange
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

Budget.CH4.I_Meth = sum(D.R_Meth(:)) .* dz .* 1e-3;
Budget.CH4.I_AOM  = sum(D.R_AOM(:)) .* dz .* 1e-3;
Budget.CH4.I_Ox   = sum(D.R_CH4Ox(:)) .* dz .* 1e-3;


Budget.CH4.I_Bubble = sum(D.R_Bubble(:)) .* dz .* 1e-3;


Budget.CH4.I_IrrigSink = I_CH4_irrig_sink;
Budget.CH4.J_TopUp = J_CH4_top_up;

Budget.CH4.Source = Budget.CH4.I_Meth;
Budget.CH4.Sink = Budget.CH4.I_AOM + Budget.CH4.I_Ox + ...
                  Budget.CH4.I_Bubble + ...
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
% OM solid budget
% =========================

OM_lab = S.OM_lab(:);
OM_ref = S.OM_ref(:);

F_OM_lab_top = F.F_lab_OM;
F_OM_ref_top = F.F_ref_OM;

F_OM_lab_bottom = A(end) .* Grid.v_solid(end) .* OM_lab(end);
F_OM_ref_bottom = A(end) .* Grid.v_solid(end) .* OM_ref(end);

% OM reaction rates are in g/gDW/yr.
I_OM_lab_reaction = trapz(z, A .* max(-Result.Rates_final.OM_lab(:), 0));
I_OM_ref_reaction = trapz(z, A .* max(-Result.Rates_final.OM_ref(:), 0));

Budget.OM.TopLab = F_OM_lab_top;
Budget.OM.TopRef = F_OM_ref_top;
Budget.OM.TopTotal = F_OM_lab_top + F_OM_ref_top;

Budget.OM.BottomLab = F_OM_lab_bottom;
Budget.OM.BottomRef = F_OM_ref_bottom;
Budget.OM.BottomTotal = F_OM_lab_bottom + F_OM_ref_bottom;

Budget.OM.ReactLab = I_OM_lab_reaction;
Budget.OM.ReactRef = I_OM_ref_reaction;
Budget.OM.ReactTotal = I_OM_lab_reaction + I_OM_ref_reaction;

Budget.OM.StorageLab = NaN;
Budget.OM.StorageRef = NaN;

if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);

    dOMlabdt = (S.OM_lab(:) - S_prev.OM_lab(:)) ./ max(dt, 1e-12);
    dOMrefdt = (S.OM_ref(:) - S_prev.OM_ref(:)) ./ max(dt, 1e-12);

    Budget.OM.StorageLab = trapz(z, A .* dOMlabdt);
    Budget.OM.StorageRef = trapz(z, A .* dOMrefdt);
end

Budget.OM.ResidualLab = Budget.OM.TopLab ...
                      - Budget.OM.BottomLab ...
                      - Budget.OM.ReactLab ...
                      - nan0(Budget.OM.StorageLab);

Budget.OM.ResidualRef = Budget.OM.TopRef ...
                      - Budget.OM.BottomRef ...
                      - Budget.OM.ReactRef ...
                      - nan0(Budget.OM.StorageRef);

% =========================
% Ca budget
% =========================

Ca = S.Ca(:);

Budget.Ca.Top = Ca(1);
Budget.Ca.Bottom = Ca(end);
Budget.Ca.Min = min(Ca);
Budget.Ca.Max = max(Ca);

% Boundary fluxes. Positive TopFlux is downward into sediment.
% Positive BottomFlux is downward out of the model domain.
Budget.Ca.TopFlux = top_solute_flux(Ca, F.Ca_top, phi, P.DCa, Grid.v_fluid, dz);
Budget.Ca.BottomFlux = bottom_solute_flux(Ca, phi, Grid.v_fluid);

% Bioirrigation/exchange term. Positive means Ca enters porewater.
Budget.Ca.I_Irrig = sum(phi .* Grid.Alpha_exchange(:) .* (F.Ca_top - Ca)) ...
                    .* dz .* 1e-3;

% R_carb_net_uM > 0 means CaCO3 precipitation, consuming Ca.
% Rates.Ca = -R_carb_net_uM.
Budget.Ca.I_CaCO3Net = sum(D.R_carb_net_uM(:)) .* dz .* 1e-3;
Budget.Ca.I_Reaction = -Budget.Ca.I_CaCO3Net;

Budget.Ca.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dCadt = (S.Ca(:) - S_prev.Ca(:)) ./ max(dt, 1e-12);
    Budget.Ca.Storage = sum(phi .* dCadt) .* dz .* 1e-3;
end

Budget.Ca.Residual = Budget.Ca.TopFlux ...
                   - Budget.Ca.BottomFlux ...
                   + Budget.Ca.I_Irrig ...
                   + Budget.Ca.I_Reaction ...
                   - nan0(Budget.Ca.Storage);

% =========================
% Fe cap diagnostics
% =========================


Budget.FeCap.Scale = D.Fe_supply_scale;

Budget.FeCap.I_FeRed_pot = D.I_FeRed_pot;

Budget.FeCap.I_Fe_supply_ext = D.I_Fe_supply_ext;


% =========================
% Redox partition diagnostics
% FV-consistent cell-sum integration.
% =========================

I_RC = sum(D.RC_uM(:)) .* dz .* 1e-3;
I_O2_C = sum(D.R_respi(:)) .* dz .* 1e-3;
I_Fe_C = sum(D.R_FeRed(:) ./ 4) .* dz .* 1e-3;
I_SO4_C = sum(2 .* D.R_SRR(:)) .* dz .* 1e-3;
I_Meth_C = sum(2 .* D.R_Meth(:)) .* dz .* 1e-3;

Budget.Redox.I_RC = I_RC;
Budget.Redox.I_O2_C = I_O2_C;
Budget.Redox.I_Fe_C = I_Fe_C;
Budget.Redox.I_SO4_C = I_SO4_C;
Budget.Redox.I_Meth_C = I_Meth_C;

Budget.Redox.frac_O2 = I_O2_C ./ max(I_RC, 1e-12);
Budget.Redox.frac_Fe = I_Fe_C ./ max(I_RC, 1e-12);
Budget.Redox.frac_SO4 = I_SO4_C ./ max(I_RC, 1e-12);
Budget.Redox.frac_Meth = I_Meth_C ./ max(I_RC, 1e-12);

% =========================
% FeOOH availability diagnostics
% =========================

Fe_gate = S.FeOOH(:) ./ max(S.FeOOH(:) + P.KFEMonod, 1e-12);
RC_after_O2 = max(D.RC_uM(:) - D.R_respi(:), 0);

Budget.FeGate.mean_all = mean(Fe_gate);
Budget.FeGate.mean_reactive = sum(Fe_gate .* RC_after_O2) ./ max(sum(RC_after_O2), 1e-12);
Budget.FeGate.max = max(Fe_gate);
Budget.FeGate.RC_after_O2_int = sum(RC_after_O2) .* dz .* 1e-3;

% =========================
% O2 budget diagnostics
% FV-consistent cell-sum integration.
% =========================

O2 = S.O2(:);

J_O2_top_down = top_solute_flux(O2, F.O2_top, phi, P.DO2, Grid.v_fluid, dz);

I_O2_irrig_source = sum(phi .* Grid.Alpha_exchange(:) .* max(F.O2_top - O2, 0)) ...
                    .* dz .* 1e-3;

I_O2_resp_sink = sum(D.R_respi(:)) .* dz .* 1e-3;


I_O2_secondary_sink = sum(D.O2_secondary_sink(:)) .* dz .* 1e-3;


I_O2_total_sink = I_O2_resp_sink + I_O2_secondary_sink;


Budget.O2.J_TopDown = max(0, J_O2_top_down);
Budget.O2.I_IrrigSource = I_O2_irrig_source;
Budget.O2.I_RespSink = I_O2_resp_sink;
Budget.O2.I_SecondarySink = I_O2_secondary_sink;
Budget.O2.I_TotalSink = I_O2_total_sink;

Budget.O2.J_TopDown = max(0, J_O2_top_down);
Budget.O2.I_IrrigSource = I_O2_irrig_source;
Budget.O2.I_RespSink = I_O2_resp_sink;

Budget.O2.Storage = NaN;
if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    dt = Result.t(end) - Result.t(end-1);
    dO2dt = (S.O2(:) - S_prev.O2(:)) ./ max(dt, 1e-12);
    Budget.O2.Storage = sum(phi .* dO2dt) .* dz .* 1e-3;
end

Budget.O2.Residual = Budget.O2.J_TopDown ...
                   + Budget.O2.I_IrrigSource ...
                   - Budget.O2.I_TotalSink ...
                   - nan0(Budget.O2.Storage);

% =========================
% Print
% =========================

fprintf('\n================ RTM Budget ================\n');

fprintf('\n--- CH4 budget, umol/cm2/yr ---\n');
fprintf('Meth source:          %.3f\n', Budget.CH4.I_Meth);
fprintf('AOM sink:             %.3f\n', Budget.CH4.I_AOM);
fprintf('O2 oxidation sink:    %.3f\n', Budget.CH4.I_Ox);
fprintf('Bubbling sink:        %.3f\n', Budget.CH4.I_Bubble);
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

fprintf('\n--- OM budget, g/cm2/yr ---\n');
fprintf('Top lab OM input:      %.5f\n', Budget.OM.TopLab);
fprintf('Top ref OM input:      %.5f\n', Budget.OM.TopRef);
fprintf('Bottom lab burial:     %.5f\n', Budget.OM.BottomLab);
fprintf('Bottom ref burial:     %.5f\n', Budget.OM.BottomRef);
fprintf('Lab OM reaction:       %.5f\n', Budget.OM.ReactLab);
fprintf('Ref OM reaction:       %.5f\n', Budget.OM.ReactRef);
fprintf('Lab OM storage:        %.5f\n', Budget.OM.StorageLab);
fprintf('Ref OM storage:        %.5f\n', Budget.OM.StorageRef);
fprintf('Lab OM residual:       %.5f\n', Budget.OM.ResidualLab);
fprintf('Ref OM residual:       %.5f\n', Budget.OM.ResidualRef);

fprintf('\n--- Fe cap diagnostics ---\n');
fprintf('Fe supply scale:       %.3f\n', Budget.FeCap.Scale);
fprintf('Potential FeRed:       %.3f umol Fe/cm2/yr\n', Budget.FeCap.I_FeRed_pot);
fprintf('External Fe cap:       %.3f umol Fe/cm2/yr\n', Budget.FeCap.I_Fe_supply_ext);

fprintf('\n--- Redox partition, carbon basis ---\n');
fprintf('I_RC:                  %.3f umol C/cm2/yr\n', Budget.Redox.I_RC);
fprintf('O2 respiration:        %.3f (%.1f%%)\n', Budget.Redox.I_O2_C,   100*Budget.Redox.frac_O2);
fprintf('Fe reduction C:        %.3f (%.1f%%)\n', Budget.Redox.I_Fe_C,   100*Budget.Redox.frac_Fe);
fprintf('SO4 reduction C:       %.3f (%.1f%%)\n', Budget.Redox.I_SO4_C,  100*Budget.Redox.frac_SO4);
fprintf('Methanogenesis C:      %.3f (%.1f%%)\n', Budget.Redox.I_Meth_C, 100*Budget.Redox.frac_Meth);

fprintf('\n--- FeOOH availability ---\n');
fprintf('Mean Fe gate:          %.4f\n', Budget.FeGate.mean_all);
fprintf('RC-weighted Fe gate:   %.4f\n', Budget.FeGate.mean_reactive);
fprintf('Max Fe gate:           %.4f\n', Budget.FeGate.max);
fprintf('RC after O2 integral:  %.3f umol C/cm2/yr\n', Budget.FeGate.RC_after_O2_int);

fprintf('\n--- O2 budget, umol/cm2/yr ---\n');
fprintf('Top O2 influx:         %.3f\n', Budget.O2.J_TopDown);
fprintf('Irrigation O2 source:  %.3f\n', Budget.O2.I_IrrigSource);
fprintf('OM respiration sink:   %.3f\n', Budget.O2.I_RespSink);
fprintf('Secondary O2 sink:     %.3f\n', Budget.O2.I_SecondarySink);
fprintf('Total O2 sink:         %.3f\n', Budget.O2.I_TotalSink);
fprintf('Storage:               %.3f\n', Budget.O2.Storage);
fprintf('Residual:              %.3f\n', Budget.O2.Residual);

fprintf('\n--- Ca budget, umol/cm2/yr ---\n');
fprintf('Ca top / bottom:      %.2f / %.2f uM\n', Budget.Ca.Top, Budget.Ca.Bottom);
fprintf('Ca min / max:         %.2f / %.2f uM\n', Budget.Ca.Min, Budget.Ca.Max);
fprintf('Top flux (+down):     %.3f\n', Budget.Ca.TopFlux);
fprintf('Bottom flux (+down):  %.3f\n', Budget.Ca.BottomFlux);
fprintf('Irrigation exchange:  %.3f\n', Budget.Ca.I_Irrig);
fprintf('CaCO3 net rxn:        %.3f\n', Budget.Ca.I_CaCO3Net);
fprintf('Ca reaction source:   %.3f\n', Budget.Ca.I_Reaction);
fprintf('Storage:              %.3f\n', Budget.Ca.Storage);
fprintf('Residual:             %.3f\n', Budget.Ca.Residual);

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

function J_bottom_down = bottom_solute_flux(C, phi, v)
% Positive downward. Convert to umol/cm2/yr.
C = C(:);
phi = phi(:);
v = v(:);

F_bottom = phi(end) .* v(end) .* C(end);
J_bottom_down = F_bottom .* 1e-3;
end