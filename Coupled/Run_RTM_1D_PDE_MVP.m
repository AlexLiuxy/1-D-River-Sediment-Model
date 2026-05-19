clear; clc;

Config = Config_Baseline();
Params = Params_Static();

Grid = Build_PDE_Grid(Config);

if Config.PDE_use_steady_IC
    load('Steady_IC_Current.mat', 'Steady');
    State0 = Initialize_PDE_State_From_Steady(Steady, Grid);
else
    State0 = Initialize_PDE_State_Generic(Config, Params, Grid);
end

Forcing = Build_Forcing_Constant(Config, Grid);

% Development spin-up / benchmark
if Config.PDE_do_spinup
    tspan_spinup = 0:Config.dt_out_spinup:Config.t_spinup;
    OutSpin = Run_PDE_Core(State0, Grid, Params, Forcing, tspan_spinup);
    StateStart = OutSpin.State_end;
else
    StateStart = State0;
end

% Transient test after spin-up
tspan_event = 0:Config.dt_out_event:Config.t_final;
Out = Run_PDE_Core(StateStart, Grid, Params, Forcing, tspan_event);

Plot_PDE_Profiles(Out, Grid);