function State = Initialize_PDE_State(Grid, Forcing, Config, Params)
% INITIALIZE_PDE_STATE
% Simple positive initial condition for PDE spin-up.
%
% This is not a formal steady state. It only gives ode15s a reasonable
% nonnegative starting point.

z = Grid.z;
n = Grid.n;

rho = Params.rho;
A_solid = rho .* max(1 - Grid.poros, 1e-6);
solid_flux = A_solid .* max(Grid.v_solid, 1e-12);

% Labile OM: depositional top value with first-order burial decay.
OM_lab_top = Forcing.F_lab_OM ./ max(solid_flux(1), 1e-12);
decay_int = cumsum((Grid.k_sed ./ max(Grid.v_solid, 1e-12)) .* Grid.dz);
State.OM_lab = max(OM_lab_top .* exp(-decay_int), 1e-12);

% Refractory OM: preserved inventory approximation.
State.OM_ref = max(Forcing.F_ref_OM ./ max(solid_flux, 1e-12), 0);

% FeOOH top inventory from external Fe flux.
% F_FeOx: mmol/m2/d -> umol/cm2/yr by *36.5.
F_FeOOH = Forcing.F_FeOx .* 36.5;
FeOOH_top = F_FeOOH ./ max(solid_flux(1), 1e-12);

if isfield(Params, 'Fe_inventory_factor')
    FeOOH_top = Params.Fe_inventory_factor .* FeOOH_top;
end

State.FeOOH = max(FeOOH_top .* exp(-z ./ 3), 1e-12);

% CaCO3 solid top inventory from depositional flux.
F_CaCO3_cm2 = Forcing.F_CaCO3 .* 1e-4;  % g/m2/yr -> g/cm2/yr
CaCO3_top = F_CaCO3_cm2 ./ max(solid_flux(1), 1e-12);
State.CaCO3 = max(CaCO3_top .* exp(-z ./ 10), 0);

% Solutes.
State.O2  = max(Forcing.O2_top .* exp(-z ./ 1.5), 0);
State.Fe2 = zeros(n,1);
State.SO4 = max(Forcing.SO4_top .* ones(n,1), 0);
State.HS  = zeros(n,1);
State.CH4 = max(Forcing.CH4_top .* ones(n,1), 0);
State.DIC = max(Forcing.DIC_top .* ones(n,1), 1e-12);
State.ALK = max(Forcing.ALK_top .* ones(n,1), 1e-12);
end