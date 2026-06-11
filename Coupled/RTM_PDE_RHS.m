function dYdt = RTM_PDE_RHS(t, Y, Grid, Forcing, Params, Config) 
% RTM_PDE_RHS
% Coupled finite-volume method-of-lines RHS for transient 1-D RTM.
%
% All dynamic species are advanced together by ode15s.
% Transport is finite-volume conservative.
% Reactions are computed by RTM_Reaction_Rates.

State = Unpack_State(Y, Grid);

% Safety only for reaction calculation.
% Formal nonnegativity should be handled by ode15s NonNegative and
% depletion-aware reaction limiters.
State = floor_state(State);

[Rates, ~] = RTM_Reaction_Rates(State, Grid, Forcing, Params, Config);

dState = struct();

% ========================================================================
% Solids: conservative mixing + burial + reaction.
% Top boundary is depositional flux.
% ========================================================================

F_OM_lab_top = Forcing.F_lab_OM;               % g/cm2/yr
F_OM_ref_top = Forcing.F_ref_OM;               % g/cm2/yr
F_FeOOH_top  = Forcing.F_FeOx .* 36.5;         % umol/cm2/yr
F_CaCO3_top  = Forcing.F_CaCO3 .* 1e-4;        % g/cm2/yr

if isfield(Params, 'Fe_inventory_factor')
    F_FeOOH_top = Params.Fe_inventory_factor .* F_FeOOH_top;
end

dState.OM_lab = solid_rhs(State.OM_lab, Grid, Rates.OM_lab, F_OM_lab_top);
dState.OM_ref = solid_rhs(State.OM_ref, Grid, Rates.OM_ref, F_OM_ref_top);
dState.FeOOH  = solid_rhs(State.FeOOH,  Grid, Rates.FeOOH,  F_FeOOH_top);
dState.CaCO3  = solid_rhs(State.CaCO3,  Grid, Rates.CaCO3,  F_CaCO3_top);

% ========================================================================
% Solutes: conservative diffusion/advection + exchange + reaction.
% Top boundary is fixed concentration through face flux.
% Bottom boundary is zero diffusive gradient with advective outflow.
% ========================================================================

dState.O2 = solute_rhs(State.O2, Grid, Params.DO2, Rates.O2, Forcing.O2_top);

dState.Fe2 = solute_rhs(State.Fe2, Grid, Params.DH2S, Rates.Fe2, Forcing.Fe_top);

dState.SO4 = solute_rhs(State.SO4, Grid, Params.DSO4, Rates.SO4, Forcing.SO4_top);

dState.HS = solute_rhs(State.HS, Grid, Params.DH2S, Rates.HS, Forcing.HS_top);

dState.CH4 = solute_rhs(State.CH4, Grid, Params.DCH4, Rates.CH4, Forcing.CH4_top);

dState.DIC = solute_rhs(State.DIC, Grid, Params.DHCO3, Rates.DIC, Forcing.DIC_top);

dState.ALK = solute_rhs(State.ALK, Grid, Params.DHCO3, Rates.ALK, Forcing.ALK_top);

dState.Ca = solute_rhs(State.Ca, Grid, Params.DCa, Rates.Ca, Forcing.Ca_top);

dYdt = Pack_State(dState);
end

% ========================================================================
% Local finite-volume functions.
% ========================================================================

function dCdt = solute_rhs(C, Grid, D, R, C_top)
% d(phi*C)/dt = -div(F) + phi*R + phi*Alpha*(C_top-C)

C = C(:);
R = R(:);

n = Grid.n;
dz = Grid.dz;
phi = Grid.poros(:);
Alpha = Grid.Alpha_exchange(:);
v = Grid.v_fluid(:);

F = zeros(n+1,1);  % face flux, positive downward, unit uM*cm/yr

% Top face: Dirichlet concentration C_top.
F(1) = -phi(1) .* D .* (C(1) - C_top) ./ (0.5 .* dz) ...
       + phi(1) .* v(1) .* upwind_top(C(1), C_top, v(1));

% Internal faces.
for i = 1:n-1
    phi_f = 0.5 .* (phi(i) + phi(i+1));
    v_f   = 0.5 .* (v(i) + v(i+1));

    if v_f >= 0
        C_up = C(i);
    else
        C_up = C(i+1);
    end

    F(i+1) = -phi_f .* D .* (C(i+1) - C(i)) ./ dz ...
             + phi_f .* v_f .* C_up;
end

% Bottom face: zero diffusive gradient, advective outflow if any.
F(n+1) = phi(end) .* v(end) .* C(end);

storage_tendency = -(F(2:end) - F(1:end-1)) ./ dz ...
                   + R ...
                   + phi .* Alpha .* (C_top - C);
% d(phi*C)/dt = -div(F) + R_bulk + phi*Alpha*(C_top-C)
% R is a bulk-volume reaction term, in uM_bulk/yr.

dCdt = storage_tendency ./ max(phi, 1e-12);
end

function dCdt = solid_rhs(C, Grid, R, F_top)
% d(A*C)/dt = -div(F) + A*R
% A = rho*(1-phi), C is solid concentration per gDW.
%
% F_top is imposed depositional flux at sediment-water interface.
% Units must match A*v*C:
%   OM/CaCO3: g/cm2/yr
%   FeOOH:    umol/cm2/yr

C = C(:);
R = R(:);

n = Grid.n;
dz = Grid.dz;

rho = Grid.rho;
phi = Grid.poros(:);
A = rho .* max(1 - phi, 1e-6);

Db = Grid.D_solid_mix(:);
v = Grid.v_solid(:);

F = zeros(n+1,1);  % positive downward

% Top depositional flux.
F(1) = F_top;

% Internal faces.
for i = 1:n-1
    A_f  = 0.5 .* (A(i) + A(i+1));
    Db_f = 0.5 .* (Db(i) + Db(i+1));
    v_f  = 0.5 .* (v(i) + v(i+1));

    if v_f >= 0
        C_up = C(i);
    else
        C_up = C(i+1);
    end

    F(i+1) = -A_f .* Db_f .* (C(i+1) - C(i)) ./ dz ...
             + A_f .* v_f .* C_up;
end

% Bottom face: burial outflow, zero diffusive gradient.
F(n+1) = A(end) .* v(end) .* C(end);

storage_tendency = -(F(2:end) - F(1:end-1)) ./ dz + A .* R;

dCdt = storage_tendency ./ max(A, 1e-12);
end

function C_face = upwind_top(C_cell, C_top, v_face)
if v_face >= 0
    C_face = C_top;
else
    C_face = C_cell;
end
end

function State = floor_state(State)
names = fieldnames(State);
for i = 1:numel(names)
    name = names{i};
    if strcmp(name, 'DIC') || strcmp(name, 'ALK') || strcmp(name, 'Ca')
        State.(name) = max(real(State.(name)), 1e-12);
    else
        State.(name) = max(real(State.(name)), 0);
    end
end
end