function Compare = RTM_Compare(Ode, Pde)
% RTM_COMPARE
% Quantitative ODE-vs-PDE benchmark table.
%
% Usage:
%   Run ODE script first and export Ode struct.
%   Run PDE and get Result.
%   Compare = RTM_Compare(Ode, Result);

P = Pde.State_final;
D = Pde.Diag_final;
G = Pde.Grid;

z_pde = G.z(:);

Compare = table();

Compare = add_row(Compare, 'OPD_cm', ...
    first_depth_below(Ode.z, Ode.O2, 1, Ode.z(end)), ...
    first_depth_below(z_pde, P.O2, 1, z_pde(end)));

Compare = add_row(Compare, 'SO4_depth_10_cm', ...
    first_depth_below(Ode.z, Ode.SO4, 10, Ode.z(end)), ...
    first_depth_below(z_pde, P.SO4, 10, z_pde(end)));

Compare = add_row(Compare, 'O2_bottom_uM', ...
    Ode.O2(end), P.O2(end));

Compare = add_row(Compare, 'SO4_bottom_uM', ...
    Ode.SO4(end), P.SO4(end));

Compare = add_row(Compare, 'Fe2_max_uM', ...
    max(Ode.Fe2), max(P.Fe2));

Compare = add_row(Compare, 'Fe2_bottom_uM', ...
    Ode.Fe2(end), P.Fe2(end));

Compare = add_row(Compare, 'HS_max_uM', ...
    max(Ode.HS), max(P.HS));

Compare = add_row(Compare, 'CH4_max_uM', ...
    max(Ode.CH4), max(P.CH4));

Compare = add_row(Compare, 'CH4_bottom_uM', ...
    Ode.CH4(end), P.CH4(end));

Compare = add_row(Compare, 'FeOOH_top_umol_g', ...
    Ode.FeOOH(1), P.FeOOH(1));

Compare = add_row(Compare, 'FeOOH_max_umol_g', ...
    max(Ode.FeOOH), max(P.FeOOH));

Compare = add_row(Compare, 'FeOOH_bottom_umol_g', ...
    Ode.FeOOH(end), P.FeOOH(end));

Compare = add_row(Compare, 'DIC_bottom_uM', ...
    Ode.DIC(end), P.DIC(end));

Compare = add_row(Compare, 'ALK_bottom_uM', ...
    Ode.ALK(end), P.ALK(end));

Compare = add_row(Compare, 'pH_bottom', ...
    Ode.pH(end), D.pH(end));

Compare = add_row(Compare, 'sigma_top5_mean', ...
    mean(Ode.sigma(Ode.z <= 5)), ...
    mean(D.sigma_carb(z_pde <= 5)));

Compare = add_row(Compare, 'sigma_bottom', ...
    Ode.sigma(end), D.sigma_carb(end));

Compare = add_row(Compare, 'CaCO3_top_pct', ...
    100 .* Ode.CaCO3(1), ...
    100 .* P.CaCO3(1));

Compare = add_row(Compare, 'CaCO3_bottom_pct', ...
    100 .* Ode.CaCO3(end), ...
    100 .* P.CaCO3(end));

Compare = add_row(Compare, 'OM_top_pct', ...
    100 .* max(getfield_or_nan(Ode, 'OM_total_top'), NaN), ...
    100 .* (P.OM_lab(1) + P.OM_ref(1)));

Compare = add_row(Compare, 'OM_bottom_pct', ...
    100 .* max(getfield_or_nan(Ode, 'OM_total_bottom'), NaN), ...
    100 .* (P.OM_lab(end) + P.OM_ref(end)));

Compare = add_row(Compare, 'OM_lab_top_pct_PDE_only', ...
    NaN, 100 .* P.OM_lab(1));

Compare = add_row(Compare, 'OM_ref_top_pct_PDE_only', ...
    NaN, 100 .* P.OM_ref(1));

% Optional integrated diagnostics if available.
if isfield(Ode, 'RC') && isfield(D, 'RC_uM')
    Compare = add_row(Compare, 'I_RC_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.RC .* 1e9) .* 1e-3, ...
        trapz(z_pde,  D.RC_uM(:)) .* 1e-3);
end

if isfield(Ode, 'R_SRR') && isfield(D, 'R_SRR')
    Compare = add_row(Compare, 'I_SRR_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_SRR) .* 1e-3, ...
        trapz(z_pde,  D.R_SRR(:)) .* 1e-3);
end

if isfield(Ode, 'Rate_Meth') && isfield(D, 'R_Meth')
    Compare = add_row(Compare, 'I_Meth_umolCH4_cm2_yr', ...
        trapz(Ode.z, Ode.Rate_Meth) .* 1e-3, ...
        trapz(z_pde,  D.R_Meth(:)) .* 1e-3);
end
% Redox partition diagnostics.
if isfield(Ode, 'R_respi') && isfield(D, 'R_respi')
    Compare = add_row(Compare, 'I_O2_resp_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.R_respi) .* 1e-3, ...
        trapz(z_pde,  D.R_respi(:)) .* 1e-3);
end

if isfield(Ode, 'R_FeRed') && isfield(D, 'R_FeRed')
    Compare = add_row(Compare, 'I_FeRed_umolFe_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeRed) .* 1e-3, ...
        trapz(z_pde,  D.R_FeRed(:)) .* 1e-3);

    Compare = add_row(Compare, 'I_FeRed_C_umolC_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeRed ./ 4) .* 1e-3, ...
        trapz(z_pde,  D.R_FeRed(:) ./ 4) .* 1e-3);
end

if isfield(Ode, 'R_AOM_actual') && isfield(D, 'R_AOM')
    Compare = add_row(Compare, 'I_AOM_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_AOM_actual) .* 1e-3, ...
        trapz(z_pde,  D.R_AOM(:)) .* 1e-3);
elseif isfield(Ode, 'R_AOM') && isfield(D, 'R_AOM')
    Compare = add_row(Compare, 'I_AOM_umolSO4_cm2_yr', ...
        trapz(Ode.z, Ode.R_AOM) .* 1e-3, ...
        trapz(z_pde,  D.R_AOM(:)) .* 1e-3);
end

if isfield(Ode, 'R_CH4Ox') && isfield(D, 'R_CH4Ox')
    Compare = add_row(Compare, 'I_CH4Ox_umolCH4_cm2_yr', ...
        trapz(Ode.z, Ode.R_CH4Ox) .* 1e-3, ...
        trapz(z_pde,  D.R_CH4Ox(:)) .* 1e-3);
end

if isfield(Ode, 'R_FeS') && isfield(D, 'R_FeS')
    Compare = add_row(Compare, 'I_FeS_umolFe_cm2_yr', ...
        trapz(Ode.z, Ode.R_FeS) .* 1e-3, ...
        trapz(z_pde,  D.R_FeS(:)) .* 1e-3);
end

if isfield(Ode, 'R_HS_Ox') && isfield(D, 'R_HSOx')
    Compare = add_row(Compare, 'I_HSOx_umolS_cm2_yr', ...
        trapz(Ode.z, Ode.R_HS_Ox) .* 1e-3, ...
        trapz(z_pde,  D.R_HSOx(:)) .* 1e-3);
end

% Derived comparison columns.
Compare.AbsDiff = Compare.PDE - Compare.ODE;
Compare.RelDiff_pct = 100 .* Compare.AbsDiff ./ max(abs(Compare.ODE), 1e-12);

disp(Compare);

end

function T = add_row(T, name, ode_val, pde_val)
newrow = table(string(name), ode_val, pde_val, ...
    'VariableNames', {'Metric','ODE','PDE'});
T = [T; newrow];
end

function depth = first_depth_below(z, x, threshold, default_depth)
z = z(:);
x = x(:);
idx = find(x < threshold, 1, 'first');
if isempty(idx)
    depth = default_depth;
else
    depth = z(idx);
end


end

function v = getfield_or_nan(S, name)
if isfield(S, name)
    v = S.(name);
else
    v = NaN;
end

end