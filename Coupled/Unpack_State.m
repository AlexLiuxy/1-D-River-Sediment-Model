function State = Unpack_State(Y, Grid)
% UNPACK_STATE
% Convert ode15s vector back to state struct.

n = Grid.n;
Y = Y(:);

expected_len = 11 * n;
if numel(Y) ~= expected_len
    error('State vector length mismatch. Expected %d, got %d.', expected_len, numel(Y));
end

i1 = 1;

State.OM_lab = Y(i1:i1+n-1); i1 = i1 + n;
State.OM_ref = Y(i1:i1+n-1); i1 = i1 + n;
State.FeOOH  = Y(i1:i1+n-1); i1 = i1 + n;
State.CaCO3  = Y(i1:i1+n-1); i1 = i1 + n;

State.O2  = Y(i1:i1+n-1); i1 = i1 + n;
State.Fe2 = Y(i1:i1+n-1); i1 = i1 + n;
State.SO4 = Y(i1:i1+n-1); i1 = i1 + n;
State.HS  = Y(i1:i1+n-1); i1 = i1 + n;
State.CH4 = Y(i1:i1+n-1); i1 = i1 + n;
State.DIC = Y(i1:i1+n-1); i1 = i1 + n;
State.ALK = Y(i1:i1+n-1);
end