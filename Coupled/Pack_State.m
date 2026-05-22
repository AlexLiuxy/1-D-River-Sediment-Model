function Y = Pack_State(State)
% PACK_STATE
% Convert state struct to one column vector for ode15s.

Y = [
    State.OM_lab(:)
    State.OM_ref(:)
    State.FeOOH(:)
    State.CaCO3(:)
    State.O2(:)
    State.Fe2(:)
    State.SO4(:)
    State.HS(:)
    State.CH4(:)
    State.DIC(:)
    State.ALK(:)
    State.Ca(:)
];
end