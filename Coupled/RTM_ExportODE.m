function Ode = RTM_ExportODE()
% RTM_EXPORTODE
% Export current ODE workspace variables into a benchmark struct.
%
% Run this from the command window after Run_RTM_1D.m finishes:
%   Ode = RTM_ExportODE();

required_vars = {'z_sed','Oxygen','Sulfate','C_Fe','C_HS','CH4', ...
                 'FeooH','C_DIC','ALK','pH','sigma_carb','CaCO3'};

caller_vars = evalin('base', 'who');

for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, caller_vars)
        error('Missing variable "%s" in base workspace. Run Run_RTM_1D.m first.', required_vars{i});
    end
end

Ode = struct();

Ode.z     = evalin('base', 'z_sed(:)');
Ode.O2    = evalin('base', 'Oxygen(:)');
Ode.SO4   = evalin('base', 'Sulfate(:)');
Ode.Fe2   = evalin('base', 'C_Fe(:)');
Ode.HS    = evalin('base', 'C_HS(:)');
Ode.CH4   = evalin('base', 'CH4(:)');
Ode.FeOOH = evalin('base', 'FeooH(:)');
Ode.DIC   = evalin('base', 'C_DIC(:)');
Ode.ALK   = evalin('base', 'ALK(:)');
Ode.pH    = evalin('base', 'pH(:)');
Ode.sigma = evalin('base', 'sigma_carb(:)');
Ode.CaCO3 = evalin('base', 'CaCO3(:)');
if ismember('C_organic_total', caller_vars)
    Corg_total = evalin('base', 'C_organic_total(:)');
elseif ismember('C_organic', caller_vars)
    Corg_total = evalin('base', 'C_organic(:)');
else
    Corg_total = NaN(size(Ode.z));
end

if ismember('C_organic_lab', caller_vars)
    Corg_lab = evalin('base', 'C_organic_lab(:)');
else
    Corg_lab = NaN(size(Ode.z));
end

if ismember('C_organic_ref', caller_vars)
    Corg_ref = evalin('base', 'C_organic_ref(:)');
else
    Corg_ref = NaN(size(Ode.z));
end

Ode.OM_total = Corg_total;
Ode.OM_lab = Corg_lab;
Ode.OM_ref = Corg_ref;

Ode.OM_total_top = Corg_total(1);
Ode.OM_total_bottom = Corg_total(end);

optional_vars = {'RC','R_respi','R_SRR','Rate_Meth', ...
                 'R_FeRed','R_FeS','R_FeOx','R_HS_Ox', ...
                 'R_AOM','R_AOM_actual','R_CH4Ox'};

for i = 1:numel(optional_vars)
    name = optional_vars{i};
    if ismember(name, caller_vars)
        Ode.(name) = evalin('base', [name, '(:)']);
    end
end

fprintf('ODE benchmark exported: %d depth nodes.\n', numel(Ode.z));
end