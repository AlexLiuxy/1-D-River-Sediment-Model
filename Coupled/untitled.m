% -------------------------------------------------------------------------
% RUN_RANGESCALED_SENSITIVITY.m
% Phase 3.5: Range-Scaled Sensitivity Analysis (Tornado Plot)
% Perturbs parameters by a fixed fraction of their feasible physical range.
% -------------------------------------------------------------------------

clc; clear; close all;

% 1. Define Parameter Space: {Name, BaseValue, MinVal, MaxVal}
% Grounded in realistic estuarine/riverine bounds
Param_Space = {
    'NPP',        200,   50,    600;   % Primary Production (g/m2/yr)
    'BE',         0.05,  0.01,  0.20;  % Burial Efficiency (fraction)
    'vbottom',    0.5,   0.1,   2.0;   % Sedimentation Rate (cm/yr)
    'F_FeOx',     1.0,   0.1,   5.0;   % Fe(III) Flux (mmol/m2/d)
    'SO4init',    200,   50,    1000;  % Boundary SO4 (uM)
    'Bioturbtop', 10,    1,     30     % Bioturbation (cm2/yr)
};

num_params = size(Param_Space, 1);
range_fraction = 0.10; % Perturb by 10% of the total physical range

% 2. Execute Baseline Run
Base_Config = Config_Baseline();
fprintf('Executing Baseline Run...\n');
try
    Base_Outputs = Run_RTM_1D(Base_Config);
catch
    error('Baseline run failed. Ensure Run_RTM_1D is stable.');
end

% Extract Baseline Targets
Base_Max_CH4   = Base_Outputs.Max_CH4;
Base_Bottom_pH = Base_Outputs.Bottom_pH;

% Preallocate
S_Max_CH4   = zeros(num_params, 1);
S_Bottom_pH = zeros(num_params, 1);
Param_Names = cell(num_params, 1);

% 3. Execute Range-Scaled Perturbation Loop
fprintf('Starting Range-Scaled Scan (%.0f%% of feasible range)...\n', range_fraction * 100);
tic;

for i = 1:num_params
    Param_Names{i} = Param_Space{i, 1};
    base_val = Param_Space{i, 2};
    min_val  = Param_Space{i, 3};
    max_val  = Param_Space{i, 4};
    
    % Calculate delta based on RANGE, not baseline
    delta_X = range_fraction * (max_val - min_val);
    perturb_val = base_val + delta_X;
    
    Run_Config = Base_Config;
    Run_Config.(Param_Names{i}) = perturb_val;
    
    fprintf('Testing %s: %.3f -> %.3f ... ', Param_Names{i}, base_val, perturb_val);
    
    try
        Outputs = Run_RTM_1D(Run_Config);
        
        % Calculate Range-Scaled Sensitivity (S)
        delta_CH4 = Outputs.Max_CH4 - Base_Max_CH4;
        % Mathematical safeguard for zero baseline methane (avoids Inf)
        denom_CH4 = max(Base_Max_CH4, 1e-6); 
        S_Max_CH4(i) = (delta_CH4 / denom_CH4) / range_fraction;
        
        delta_pH = Outputs.Bottom_pH - Base_Bottom_pH;
        S_Bottom_pH(i) = (delta_pH / Base_Bottom_pH) / range_fraction;
        
        fprintf('Done.\n');
    catch ME
        fprintf('FAILED (Stiff ODE). S assigned as NaN.\n');
        S_Max_CH4(i)   = NaN;
        S_Bottom_pH(i) = NaN;
    end
end
exec_time = toc;
fprintf('Scan Complete in %.2f seconds.\n', exec_time);

% 4. Visualization (Tornado Plots)
figure('Name', 'Range-Scaled Sensitivity Analysis', 'Color', 'w', 'Position', [150, 150, 1000, 450]);

% Subplot 1: Sensitivity of Max CH4
subplot(1,2,1);
[sorted_S_CH4, idx_CH4] = sort(S_Max_CH4, 'ascend');
sorted_Names_CH4 = Param_Names(idx_CH4);

barh(sorted_S_CH4, 'FaceColor', [0.85 0.32 0.09], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_CH4, 'TickLabelInterpreter', 'none');
xlabel('Range-Scaled Sensitivity (S_{range})');
title('Sensitivity of Max CH_4');
xline(0, 'k--', 'LineWidth', 1.5);
grid on;

% Subplot 2: Sensitivity of Bottom pH
subplot(1,2,2);
[sorted_S_pH, idx_pH] = sort(S_Bottom_pH, 'ascend');
sorted_Names_pH = Param_Names(idx_pH);

barh(sorted_S_pH, 'FaceColor', [0.0 0.44 0.74], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_pH, 'TickLabelInterpreter', 'none');
xlabel('Range-Scaled Sensitivity (S_{range})');
title('Sensitivity of Bottom pH');
xline(0, 'k--', 'LineWidth', 1.5);
grid on;

sgtitle(sprintf('Range-Scaled Sensitivity (+%.0f%% of Physical Bound)', range_fraction*100), 'FontWeight', 'bold');