function Summary = Summarize_PDE_Result(Result)
% SUMMARIZE_PDE_RESULT
% Print compact quantitative diagnostics for the FV-MOL PDE result.

Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;

z = Grid.z(:);

Summary = struct();

Summary.OPD = first_depth_below(z, S.O2, 1, z(end));
Summary.SO4_Depth_10 = first_depth_below(z, S.SO4, 10, z(end));
Summary.SO4_Depth_1  = first_depth_below(z, S.SO4, 1, z(end));

Summary.OM_top_pct = 100 .* (S.OM_lab(1) + S.OM_ref(1));
Summary.OM_bottom_pct = 100 .* (S.OM_lab(end) + S.OM_ref(end));

Summary.O2_bottom = S.O2(end);
Summary.SO4_bottom = S.SO4(end);
Summary.Fe2_max = max(S.Fe2);
Summary.Fe2_bottom = S.Fe2(end);
Summary.HS_max = max(S.HS);
Summary.CH4_max = max(S.CH4);
Summary.CH4_bottom = S.CH4(end);
if isfield(D, 'CH4_bubble_threshold_top')
    Summary.CH4_bubble_threshold_top = D.CH4_bubble_threshold_top;
    Summary.CH4_bubble_threshold_bottom = D.CH4_bubble_threshold_bottom;
else
    Summary.CH4_bubble_threshold_top = NaN;
    Summary.CH4_bubble_threshold_bottom = NaN;
end

Summary.FeOOH_top = S.FeOOH(1);
Summary.FeOOH_max = max(S.FeOOH);
Summary.FeOOH_bottom = S.FeOOH(end);

Summary.DIC_bottom = S.DIC(end);
Summary.ALK_bottom = S.ALK(end);
Summary.Ca_top = S.Ca(1);
Summary.Ca_bottom = S.Ca(end);
Summary.pH_bottom = D.pH(end);
Summary.sigma_top5 = mean(D.sigma_carb(z <= 5));
Summary.sigma_bottom = D.sigma_carb(end);

Summary.CaCO3_top_pct = 100 .* S.CaCO3(1);
Summary.CaCO3_bottom_pct = 100 .* S.CaCO3(end);
Summary.CaCO3_net_int = D.I_CaCO3_net;

Summary.I_RC = D.I_RC;
Summary.I_respi = D.I_respi;
Summary.I_FeRed_C = D.I_FeRed_C;
Summary.I_SRR = D.I_SRR;
Summary.I_AOM = D.I_AOM;
Summary.I_Meth = D.I_Meth;
if isfield(D, 'I_Bubble')
    Summary.I_Bubble = D.I_Bubble;
else
    Summary.I_Bubble = 0;
end

Summary.redox_closure = ...
    (D.I_respi + D.I_FeRed_C + 2 .* D.I_SRR + 2 .* D.I_Meth) ./ ...
    max(D.I_RC, 1e-12);

Summary.Fe_supply_scale = D.Fe_supply_scale;
Summary.I_FeRed_pot = D.I_FeRed_pot;
Summary.I_Fe_supply_ext = D.I_Fe_supply_ext;

Summary.last_step_change = NaN;

if isfield(Result, 'Y') && size(Result.Y,1) >= 2
    S_prev = Unpack_State(Result.Y(end-1,:).', Grid);
    S_now  = S;

    names = {'O2','Fe2','SO4','HS','CH4','DIC','ALK','Ca','FeOOH','CaCO3'};
    changes = zeros(numel(names),1);

    for i = 1:numel(names)
        name = names{i};
        a = S_now.(name)(:);
        b = S_prev.(name)(:);
        changes(i) = max(abs(a - b)) ./ max(max(abs(b)), 1);
    end

    Summary.last_step_change = max(changes);
    Summary.last_step_change_by_species = array2table(changes(:).', ...
        'VariableNames', names);
end

fprintf('\n================ PDE Summary ================\n');
fprintf('OPD <1 uM:                 %.2f cm\n', Summary.OPD);
fprintf('SO4 <10 uM depth:          %.2f cm\n', Summary.SO4_Depth_10);
fprintf('SO4 bottom:                %.2f uM\n', Summary.SO4_bottom);

fprintf('\nOM top / bottom:            %.3f / %.3f %%gDW\n', ...
    Summary.OM_top_pct, Summary.OM_bottom_pct);

fprintf('\nFe2 max / bottom:           %.2f / %.2f uM\n', ...
    Summary.Fe2_max, Summary.Fe2_bottom);
fprintf('FeOOH top / max / bottom:   %.2f / %.2f / %.2f umol/g\n', ...
    Summary.FeOOH_top, Summary.FeOOH_max, Summary.FeOOH_bottom);
fprintf('Fe supply scale:            %.3f\n', Summary.Fe_supply_scale);

fprintf('\nHS max:                     %.2f uM\n', Summary.HS_max);
fprintf('CH4 max / bottom:           %.2f / %.2f uM\n', ...
    Summary.CH4_max, Summary.CH4_bottom);
fprintf('CH4 bubble threshold top/bottom: %.2f / %.2f uM\n', ...
    Summary.CH4_bubble_threshold_top, Summary.CH4_bubble_threshold_bottom);

fprintf('\nDIC bottom:                 %.2f uM\n', Summary.DIC_bottom);
fprintf('ALK bottom:                 %.2f uM\n', Summary.ALK_bottom);
fprintf('Ca top / bottom:            %.2f / %.2f uM\n', ...
    Summary.Ca_top, Summary.Ca_bottom);
fprintf('pH bottom:                  %.3f\n', Summary.pH_bottom);
fprintf('Mean sigma top 5 cm:        %.4f\n', Summary.sigma_top5);
fprintf('CaCO3 top / bottom:         %.4f / %.4f %%gDW\n', ...
    Summary.CaCO3_top_pct, Summary.CaCO3_bottom_pct);

fprintf('\nIntegrated RC:              %.3f umol C/cm2/yr\n', Summary.I_RC);
fprintf('Integrated O2 resp:         %.3f\n', Summary.I_respi);
fprintf('Integrated Fe reduction C:  %.3f\n', Summary.I_FeRed_C);
fprintf('Integrated SRR:             %.3f umol SO4/cm2/yr\n', Summary.I_SRR);
fprintf('Integrated AOM:             %.3f umol SO4/cm2/yr\n', Summary.I_AOM);
fprintf('Integrated methanogenesis:  %.3f umol CH4/cm2/yr\n', Summary.I_Meth);
fprintf('Integrated bubbling loss:   %.3f umol CH4/cm2/yr\n', Summary.I_Bubble);
fprintf('Redox closure diagnostic:   %.3f\n', Summary.redox_closure);

if ~isnan(Summary.last_step_change)
    fprintf('\nMax relative change last output step: %.3e\n', Summary.last_step_change);
    disp(Summary.last_step_change_by_species);
end

fprintf('=============================================\n\n');

end

function depth = first_depth_below(z, x, threshold, default_depth)
idx = find(x(:) < threshold, 1, 'first');
if isempty(idx)
    depth = default_depth;
else
    depth = z(idx);
end
end