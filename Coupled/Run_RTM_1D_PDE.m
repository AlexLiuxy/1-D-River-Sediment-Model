function Result = Run_RTM_1D_PDE(Config, Params, do_plot)
% RUN_RTM_1D_PDE
% First coupled FV-MOL transient RTM driver.
%
tic

if nargin < 1 || isempty(Config)
    Config = Config_Baseline();
end
if nargin < 2 || isempty(Params)
    Params = Params_Static();
end
if nargin < 3 || isempty(do_plot)
    do_plot = true;
end

% First PDE benchmark settings.
if ~isfield(Config, 't_spinup')
    Config.t_spinup = 300;  % yr
end
if ~isfield(Config, 'dt_out_spinup')
    Config.dt_out_spinup = 2;  % yr
end

Grid = Build_Grid_1D(Config, Params);
Forcing = Build_Forcing_PDE(Config);

State0 = Initialize_PDE_State(Grid, Forcing, Config, Params);
Y0 = Pack_State(State0);

tspan = 0:Config.dt_out_spinup:Config.t_spinup;

nonnegative_idx = 1:numel(Y0);

opts = odeset( ...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-8, ...
    'NonNegative', nonnegative_idx, ...
    'MaxStep', 1);

rhs = @(t,Y) RTM_PDE_RHS(t, Y, Grid, Forcing, Params, Config);

fprintf('Starting coupled FV-MOL PDE spin-up...\n');
[tout, Yout] = ode15s(rhs, tspan, Y0, opts);
fprintf('PDE spin-up finished.\n');

State_final = Unpack_State(Yout(end,:).', Grid);
State_final = floor_output_state(State_final);

[Rates_final, Diag_final] = RTM_Reaction_Rates(State_final, Grid, Forcing, Params, Config);

Result = struct();
Result.t = tout;
Result.Y = Yout;
Result.Grid = Grid;
Result.Forcing = Forcing;
Result.Params = Params;
Result.Config = Config;
Result.State_final = State_final;
Result.Rates_final = Rates_final;
Result.Diag_final = Diag_final;

Result.Summary = Summarize_PDE_Result(Result);
% Result.Budget = RTM_Budget(Result);
% Result.Validation = RTM_Validation_Metrics(Result);

if do_plot
    Plot_PDE_Result(Result);
end

toc
end

function State = floor_output_state(State)
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