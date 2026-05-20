function [Rates, Diag] = RTM_Reaction_Rates(State, Grid, Forcing, Params, Config)
% RTM_REACTION_RATES
% Coupled reaction-rate calculator for the transient FV-MOL RTM.
%
% This file intentionally extracts mechanism logic from the stable ODE model,
% but does NOT copy bvp4c-specific flux variables or Picard iteration logic.
%
% Required State fields, all column vectors with length Grid.n:
%   Solids:
%     OM_lab      g/gDW
%     OM_ref      g/gDW
%     FeOOH       umol/gDW
%     CaCO3       g/gDW
%
%   Solutes:
%     O2          uM
%     Fe2         uM
%     SO4         uM
%     HS          uM total sulfide
%     CH4         uM
%     DIC         uM
%     ALK         uM
%
% Returned Rates fields have units:
%   Solids: same concentration unit per yr
%   Solutes: uM/yr
%
% Transport, burial, diffusion, bioirrigation, and boundary conditions
% should be handled outside this function in RTM_PDE_RHS.m.

% ---------- vector safety ----------
z   = Grid.z(:);
phi = Grid.poros(:);
ks  = Grid.k_sed(:);
n   = numel(z);

OM_lab = col(State.OM_lab, n, 'OM_lab');
OM_ref = col(State.OM_ref, n, 'OM_ref');
O2     = max(col(State.O2,    n, 'O2'),    0);
FeOOH  = max(col(State.FeOOH, n, 'FeOOH'), 0);
Fe2    = max(col(State.Fe2,   n, 'Fe2'),   0);
SO4    = max(col(State.SO4,   n, 'SO4'),   0);
HS     = max(col(State.HS,    n, 'HS'),    0);
CH4    = max(col(State.CH4,   n, 'CH4'),   0);
DIC    = max(col(State.DIC,   n, 'DIC'),   1e-12);
ALK    = max(col(State.ALK,   n, 'ALK'),   1e-12);
CaCO3  = max(col(State.CaCO3, n, 'CaCO3'), 0);

% ---------- scalar parameters ----------
rho = pick_scalar(Params, Config, Forcing, {'rho'}, 2.73);

k_O2       = pick_scalar(Params, Config, Forcing, {'k_O2'}, 2);
k_SO4      = pick_scalar(Params, Config, Forcing, {'k_SO4'}, 20);
KFEMonod   = pick_scalar(Params, Config, Forcing, {'KFEMonod'}, 1000);
K_HS       = pick_scalar(Params, Config, Forcing, {'K_HS'}, 7);

kFeOx      = pick_scalar(Params, Config, Forcing, {'kFeOx'}, 10);
kFeS       = pick_scalar(Params, Config, Forcing, {'kFeS'}, 1);
Kreox      = pick_scalar(Params, Config, Forcing, {'Kreox'}, 500);

K_CH4_SO4  = pick_scalar(Params, Config, Forcing, {'K_CH4_SO4'}, 100);
K_CH4_O2   = pick_scalar(Params, Config, Forcing, {'K_CH4_O2'}, 1);
k_AOM      = pick_scalar(Params, Config, Forcing, {'k_AOM'}, 0.2);
k_CH4_O2   = pick_scalar(Params, Config, Forcing, {'k_aerobic_CH4'}, 6);

Q10        = pick_scalar(Params, Config, Forcing, {'Q10'}, 2);
T_ref      = pick_scalar(Params, Config, Forcing, {'T_ref'}, 25);
T          = pick_scalar(Forcing, Config, Params, {'T', 'T_future'}, 25);
Salinity   = pick_scalar(Forcing, Config, Params, {'Salinity'}, 0.1);

F_FeOx     = pick_scalar(Forcing, Config, Params, {'F_FeOx'}, 0);
Ca         = pick_scalar(Forcing, Config, Params, {'Ca_top', 'Calcium'}, 1000);

k_ref_factor = pick_scalar(Config, Params, Forcing, {'k_ref_factor'}, 0);

% Carbonate kinetic parameters.
Ksp_ca     = pick_scalar(Params, Config, Forcing, {'Ksp_ca'}, 4.5e5);
k_calcite  = pick_scalar(Params, Config, Forcing, {'k_calcite'}, 1);
k_dis1     = pick_scalar(Params, Config, Forcing, {'k_calcite_dis1'}, 0.005);
k_dis2     = pick_scalar(Params, Config, Forcing, {'k_calcite_dis2'}, 10);
n_form     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO31'}, 1.76);
n_dis1     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO32'}, 0.11);
n_dis2     = pick_scalar(Params, Config, Forcing, {'n_power_CaCO33'}, 4);
gamma_Ca   = pick_scalar(Params, Config, Forcing, {'Calcium_activity'}, 0.6);
gamma_CO3  = pick_scalar(Params, Config, Forcing, {'CO3_activity'}, 0.6);

% Numerical limiter parameters.
SO4_cut = pick_scalar(Config, Params, Forcing, {'SO4_cut'}, 1e-6);

Temp_factor = Q10.^((T - T_ref) ./ 10);

% ========================================================================
% 1. Carbonate speciation first, because HS free fraction needs pH.
% ========================================================================

Carb = RTM_Carbonate(ALK, DIC, T, Salinity);

pH  = Carb.pH;
CO3 = max(Carb.CO3, 1e-12);

sigma_carb = (gamma_Ca .* Ca .* gamma_CO3 .* CO3) ./ Ksp_ca - 1;

% ========================================================================
% 2. OM degradation and residual-carbon cascade.
% ========================================================================

% OM_lab is the redox-active pool.
% Unit follows current ODE convention:
%   RC_mol = mol C / cm3 bulk sediment / yr
%   RC_uM  = umol C / L bulk sediment / yr
RC_mol = Temp_factor .* ks .* OM_lab .* rho .* max(1 - phi, 1e-6) ./ 12;
RC_uM  = RC_mol .* 1e9;

% Solid OM reaction terms.
R_OM_lab = -Temp_factor .* ks .* OM_lab;
R_OM_ref = -Temp_factor .* k_ref_factor .* ks .* OM_ref;

% O2 branch.
R_respi = RC_uM .* O2 ./ max(O2 + k_O2, 1e-12);
R_respi = min(max(R_respi, 0), RC_uM);

RC_after_O2 = max(RC_uM - R_respi, 0);

% Fe branch.
Fe_gate = FeOOH ./ max(FeOOH + KFEMonod, 1e-12);
C_to_Fe_pot = RC_after_O2 .* Fe_gate;
R_FeRed_pot = 4 .* C_to_Fe_pot;  % umol Fe / L / yr

I_FeRed_pot = trapz(z, R_FeRed_pot) .* 1e-3;  % umol Fe / cm2 / yr
I_Fe_supply_ext = F_FeOx .* 36.5;             % mmol/m2/d -> umol/cm2/yr

Fe_supply_scale = min(1, I_Fe_supply_ext ./ max(I_FeRed_pot, 1e-12));
R_FeRed = R_FeRed_pot .* Fe_supply_scale;

C_to_Fe = R_FeRed ./ 4;
RC_after_Fe = max(RC_after_O2 - C_to_Fe, 0);

% SO4 branch with active depletion limiter.
SO4_pos = max(SO4, 0);
f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
f_SRR(SO4_pos <= SO4_cut) = 0;

R_SRR_pot = 0.5 .* RC_after_Fe;
R_SRR = R_SRR_pot .* f_SRR;

% AOM uses current CH4 and current SO4 in the coupled PDE RHS.
f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
f_AOM(SO4_pos <= SO4_cut) = 0;

R_AOM = k_AOM .* CH4 .* f_AOM;

% Methanogenesis branch.
RC_after_SO4 = max(RC_after_Fe - 2 .* R_SRR, 0);
R_Meth = 0.5 .* RC_after_SO4;

% Aerobic methane oxidation.
R_CH4Ox = k_CH4_O2 .* CH4 .* O2 ./ max(O2 + K_CH4_O2, 1e-12);

% ========================================================================
% 3. Fe/S secondary reactions.
% ========================================================================

HS_free = HS ./ (1 + (10.^(6 - pH)) ./ max(K_HS, 1e-12));
HS_free = max(HS_free, 0);

R_FeOx = kFeOx .* Fe2 .* O2;
R_FeS  = kFeS  .* Fe2 .* HS_free;
R_HSOx = Kreox .* HS_free .* O2;

% Convert FeOOH bulk-volume reaction rates to solid concentration rates:
% R_FeRed / R_FeOx: umol/L_bulk/yr
% FeOOH:            umol/gDW
% A = rho*(1-phi):  gDW/cm3_bulk
%
% umol/L_bulk/yr * 1e-3 = umol/cm3_bulk/yr
% solid rate = umol/cm3_bulk/yr / A
bulk_to_solid_Fe = 1e-3 ./ (rho .* max(1 - phi, 1e-6));

R_FeRed_solid = R_FeRed .* bulk_to_solid_Fe;
R_FeOx_solid  = R_FeOx  .* bulk_to_solid_Fe;

% ========================================================================
% 4. Carbonate precipitation/dissolution.
% ========================================================================

solid_to_uM_carb = 1e3 .* 1e6 .* 1e-2 .* rho .* max(1 - phi, 1e-6);
uM_to_solid_carb = 1 ./ solid_to_uM_carb;

R_carb_form_solid = zeros(n,1);
R_carb_dis_solid  = zeros(n,1);

idx_form = sigma_carb > 0;
R_carb_form_solid(idx_form) = ...
    abs(sigma_carb(idx_form)).^n_form .* k_calcite .* uM_to_solid_carb(idx_form);

idx_dis1 = sigma_carb < 0 & sigma_carb > -0.2;
R_carb_dis_solid(idx_dis1) = ...
    abs(sigma_carb(idx_dis1)).^n_dis1 .* k_dis1 .* CaCO3(idx_dis1);

idx_dis2 = sigma_carb <= -0.2;
R_carb_dis_solid(idx_dis2) = ...
    abs(sigma_carb(idx_dis2)).^n_dis2 .* k_dis2 .* CaCO3(idx_dis2);

% Positive means net CaCO3 precipitation into solid phase.
R_carb_net_solid = R_carb_form_solid - R_carb_dis_solid;
R_carb_net_uM    = R_carb_net_solid .* solid_to_uM_carb;

% ========================================================================
% 5. Final explicit reaction source/sink terms for PDE RHS.
% ========================================================================

Rates = struct();

% Solids.
Rates.OM_lab = R_OM_lab;
Rates.OM_ref = R_OM_ref;
Rates.FeOOH  = -R_FeRed_solid + R_FeOx_solid;
Rates.CaCO3  = R_carb_net_solid;

% Solutes.
% For first PDE benchmark, O2 follows the current ODE-compatible primary
% OM respiration sink. Secondary O2 sinks are stored in diagnostics below.
Rates.O2  = -R_respi;

Rates.Fe2 = +R_FeRed - R_FeOx - R_FeS;

% Include sulfide oxidation as sulfate regeneration in the transient sulfur
% ledger. This is more mass-consistent than the current steady FV SO4 solver.
Rates.SO4 = -R_SRR - R_AOM + R_HSOx;

Rates.HS  = +R_SRR + R_AOM - R_FeS - R_HSOx;

Rates.CH4 = +R_Meth - R_AOM - R_CH4Ox;

% DIC/ALK ledger follows explicit reaction signs:
% carbonate precipitation consumes DIC/ALK; dissolution releases them.
Rates.DIC = +R_respi ...
            +R_FeRed ./ 4 ...
            +2 .* R_SRR ...
            +R_Meth ...
            +R_CH4Ox ...
            +R_AOM ...
            -R_carb_net_uM;

Rates.ALK = +0.5 .* R_FeRed ...
            +2 .* R_SRR ...
            +2 .* R_AOM ...
            -2 .* R_FeS ...
            -R_HSOx ...
            -2 .* R_carb_net_uM;

% ========================================================================
% 6. Diagnostics.
% ========================================================================

Diag = struct();

Diag.pH = pH;
Diag.CO3 = CO3;
Diag.H2CO3 = Carb.H2CO3;
Diag.sigma_carb = sigma_carb;

Diag.RC_mol = RC_mol;
Diag.RC_uM = RC_uM;
Diag.R_respi = R_respi;
Diag.R_FeRed = R_FeRed;
Diag.R_SRR = R_SRR;
Diag.R_AOM = R_AOM;
Diag.R_Meth = R_Meth;
Diag.R_CH4Ox = R_CH4Ox;
Diag.R_FeOx = R_FeOx;
Diag.R_FeS = R_FeS;
Diag.R_HSOx = R_HSOx;

Diag.R_carb_form_solid = R_carb_form_solid;
Diag.R_carb_dis_solid = R_carb_dis_solid;
Diag.R_carb_net_solid = R_carb_net_solid;
Diag.R_carb_net_uM = R_carb_net_uM;

Diag.HS_free = HS_free;
Diag.f_SRR = f_SRR;
Diag.f_AOM = f_AOM;
Diag.R_SRR_pot = R_SRR_pot;
Diag.R_FeRed_pot = R_FeRed_pot;
Diag.Fe_supply_scale = Fe_supply_scale;
Diag.I_FeRed_pot = I_FeRed_pot;
Diag.I_Fe_supply_ext = I_Fe_supply_ext;

Diag.O2_secondary_sink = R_FeOx + R_HSOx + R_CH4Ox;

Diag.I_RC = trapz(z, RC_uM) .* 1e-3;
Diag.I_respi = trapz(z, R_respi) .* 1e-3;
Diag.I_FeRed_C = trapz(z, R_FeRed ./ 4) .* 1e-3;
Diag.I_SRR = trapz(z, R_SRR) .* 1e-3;
Diag.I_AOM = trapz(z, R_AOM) .* 1e-3;
Diag.I_Meth = trapz(z, R_Meth) .* 1e-3;
Diag.I_CaCO3_net = trapz(z, R_carb_net_uM) .* 1e-3;

end

% ========================================================================
% Local helpers.
% ========================================================================

function v = col(x, n, name)
    v = x(:);
    if numel(v) ~= n
        error('State.%s must have length %d, but got length %d.', name, n, numel(v));
    end
    v = real(v);
end

function val = pick_scalar(S1, S2, S3, names, default_val)
    val = default_val;
    structs = {S1, S2, S3};
    for is = 1:numel(structs)
        S = structs{is};
        if isempty(S) || ~isstruct(S)
            continue
        end
        for iname = 1:numel(names)
            nm = names{iname};
            if isfield(S, nm) && ~isempty(S.(nm))
                candidate = S.(nm);
                if isscalar(candidate)
                    val = candidate;
                    return
                end
            end
        end
    end
end