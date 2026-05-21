% -------------------------------------------------------------------------
% RUN_RANGESCALED_SENSITIVITY.m
% Sensitivity Analysis (Tornado Plot)
% Perturbs parameters by a fixed fraction of their feasible physical range.
% -------------------------------------------------------------------------

clc; clear; close all;

% 1. Define Parameter Space: {Name, BaseValue, MinVal, MaxVal}
% Grounded in realistic estuarine/riverine bounds
Param_Space = {
    'NPP',        200,   50,    600;   % Primary Production (g/m2/yr)
    'BE',         0.1,  0.03,  0.15;  % Burial Efficiency (fraction)
    'vbottom',    0.5,   0.1,   2.0;   % Sedimentation Rate (cm/yr)
    'F_FeOx',     2,   0.1,   5.0;   % Fe(III) Flux (mmol/m2/d)
    'F_CaCO3',        10,    1,     30;
    'SO4init',    200,   50,    500;  % Boundary SO4 (uM)
    'HCO3init',       950,   500,   2000;
    'DICinit',        1000,  700,   2500;
    'Calcium',        1000,  200,   2000;
    'Bioturbtop', 10,    1,     30;     % Bioturbation (cm2/yr)
    'k_SO4',      20,    5,     30;     % uM
    'DSO4',      310,   100,   500;    % cm2 / yr
    'DCH4',      300,   100,   500;    % cm2 / yr
    'DH2S',      300,   429,   650;    % cm2 / yr
    'DO2',      300,   370,   730;    % cm2 / yr
    'Kreox',      500,   10,   1000;    % 1 / umol / L / yr
    'kFeOx',      10,   1,   100;    % 1 / umol / L / yr
    'kFeS',       10,   1,   100;    % 1 / umol / L / yr
    'K_CH4_SO4',  100,   10,   500;    % uM
    'k_AOM',      1,   0.1,   10;    % 1 / yr
    'k_aerobic_CH4',      6,   1,   10 ;   % 1 / yr
    % 'k_calcite',      1,     0.2,   5;
    % 'k_calcite_dis1', 0.005, 0.001, 0.05;
    'k_calcite_dis2', 10,    1,     50;
};

num_params = size(Param_Space, 1);
range_fraction = 0.10; % Perturb by 10% of the total physical range

% 2. Execute Baseline Run
Base_Config = Config_Baseline();
Base_Params = Params_Static(); 
fprintf('Executing Baseline Run...\n');
try
    Base_Outputs = Run_RTM_1D(Base_Config, Base_Params); 
catch
    error('Baseline run failed.');
end

% % Extract Baseline Targets
% Base_Max_CH4   = Base_Outputs.Max_CH4;
% Base_Bottom_pH = Base_Outputs.Bottom_pH;
% 
% % Preallocate
% S_Max_CH4   = zeros(num_params, 1);
% S_Bottom_pH = zeros(num_params, 1);
% S_Org_Burial  = zeros(num_params, 1);
% S_ALK_Bottom  = zeros(num_params, 1);
% S_Carb_Burial = zeros(num_params, 1);
% S_CH4_Flux    = zeros(num_params, 1);
% Param_Names = cell(num_params, 1);
% 
% % 3. Execute Range-Scaled Perturbation Loop
% fprintf('Starting Range-Scaled Scan (%.0f%% of feasible range)...\n', range_fraction * 100);
% tic;
% 
% for i = 1:num_params
%     Param_Names{i} = Param_Space{i, 1};
%     base_val = Param_Space{i, 2};
%     min_val  = Param_Space{i, 3};
%     max_val  = Param_Space{i, 4};
%     
%     % Calculate delta based on RANGE, not baseline
%     delta_X = range_fraction * (max_val - min_val);
%     perturb_val = base_val + delta_X;
%     
%     Run_Config = Base_Config;
%     Run_Config.(Param_Names{i}) = perturb_val;
%     
%     fprintf('Testing %s: %.3f -> %.3f ... ', Param_Names{i}, base_val, perturb_val);
%     
%     try
%         Outputs = Run_RTM_1D(Run_Config);
%         
%         % Calculate Range-Scaled Sensitivity (S)
%         delta_CH4 = Outputs.Max_CH4 - Base_Max_CH4;
%         % Mathematical safeguard for zero baseline methane (avoids Inf)
%         denom_CH4 = max(Base_Max_CH4, 1e-6); 
%         S_Max_CH4(i) = (delta_CH4 / denom_CH4) / range_fraction;
%         
%         delta_pH = Outputs.Bottom_pH - Base_Bottom_pH;
%         S_Bottom_pH(i) = (delta_pH / Base_Bottom_pH) / range_fraction;
%         
%         fprintf('Done.\n');
%     catch ME
%         fprintf('FAILED (Stiff ODE). S assigned as NaN.\n');
%         S_Max_CH4(i)   = NaN;
%         S_Bottom_pH(i) = NaN;
%     end
% end
% exec_time = toc;
% fprintf('Scan Complete in %.2f seconds.\n', exec_time);
% 
% % 4. Visualization (Tornado Plots)
% figure('Name', 'Range-Scaled Sensitivity Analysis', 'Color', 'w', 'Position', [150, 150, 1000, 450]);
% 
% % Subplot 1: Sensitivity of Max CH4
% subplot(1,2,1);
% [sorted_S_CH4, idx_CH4] = sort(S_Max_CH4, 'ascend');
% sorted_Names_CH4 = Param_Names(idx_CH4);
% 
% barh(sorted_S_CH4, 'FaceColor', [0.85 0.32 0.09], 'EdgeColor', 'k');
% set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_CH4, 'TickLabelInterpreter', 'none');
% xlabel('Range-Scaled Sensitivity (S_{range})');
% title('Sensitivity of Max CH_4');
% xline(0, 'k--', 'LineWidth', 1.5);
% grid on;
% 
% % Subplot 2: Sensitivity of Bottom pH
% subplot(1,2,2);
% [sorted_S_pH, idx_pH] = sort(S_Bottom_pH, 'ascend');
% sorted_Names_pH = Param_Names(idx_pH);
% 
% barh(sorted_S_pH, 'FaceColor', [0.0 0.44 0.74], 'EdgeColor', 'k');
% set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names_pH, 'TickLabelInterpreter', 'none');
% xlabel('Range-Scaled Sensitivity (S_{range})');
% title('Sensitivity of Bottom pH');
% xline(0, 'k--', 'LineWidth', 1.5);
% grid on;
% 
% sgtitle(sprintf('Range-Scaled Sensitivity (+%.0f%% of Physical Bound)', range_fraction*100), 'FontWeight', 'bold');

% --- 2. 预分配空间与基线提取 ---
Base_Org_Bottom = Base_Outputs.Org_Bottom;
% Base_Org_Top  = Base_Outputs.Org_Top;
Base_ALK = Base_Outputs.ALK_Bottom;
Base_pH  = Base_Outputs.pH_Bottom;
Base_CH4 = Base_Outputs.CH4_Bottom;
Base_OPD   = Base_Outputs.OPD;
Base_SO4D  = Base_Outputs.SO4_Depth;
% Base_CH4D  = Base_Outputs.CH4_Onset;
% Base_SigD  = Base_Outputs.Sigma0_Depth;
% Base_CaD   = Base_Outputs.CaCO3_Front;
B_Meth  = Base_Outputs.Integ_Meth;
B_ALK5  = Base_Outputs.ALK_Bot5;
B_Sig5  = Base_Outputs.Sigma_Top5;
B_Ca5   = Base_Outputs.CaCO3_Top5;

% 无论你有多少个 Output，都先建好全零数组
S_OrgB = zeros(num_params, 1);
% S_OrgT = zeros(num_params,1);
S_ALK = zeros(num_params, 1);
S_pH  = zeros(num_params, 1);
S_CH4 = zeros(num_params, 1);
S_OPD  = zeros(num_params,1); 
S_SO4D = zeros(num_params,1); 
% S_CH4D = zeros(num_params,1); 
% S_SigD = zeros(num_params,1); 
% S_CaD  = zeros(num_params,1);
S_ALK5 = zeros(num_params, 1);
S_Sig5 = zeros(num_params,1); 
S_Ca5  = zeros(num_params,1); 
S_Meth = zeros(num_params,1);



Param_Names = cell(num_params, 1);

% --- 3. 执行扰动循环 ---
fprintf('Starting Range-Scaled Scan...\n');
tic;
for i = 1:num_params
    Param_Names{i} = Param_Space{i, 1};
    base_val = Param_Space{i, 2};
    min_val  = Param_Space{i, 3};
    max_val  = Param_Space{i, 4};
    
    delta_X = range_fraction * (max_val - min_val);
    perturb_val = base_val + delta_X;
    
    Run_Config = Base_Config;
    Run_Params = Base_Params;
    
    % Auto-Route the parameter to the correct struct
    if isfield(Run_Config, Param_Names{i})
        Run_Config.(Param_Names{i}) = perturb_val;
    elseif isfield(Run_Params, Param_Names{i})
        Run_Params.(Param_Names{i}) = perturb_val;
    else
        error(['Parameter ', Param_Names{i}, ' not found in Config or Params.']);
    end
    
    fprintf('Testing %s: %.3f -> %.3f ... \n', Param_Names{i}, base_val, perturb_val);
    
    try
        Outputs = Run_RTM_1D(Run_Config, Run_Params); % Pass both!

        % 计算各个指标的归一化敏感度
        S_OrgB(i) = ((Outputs.Org_Bottom - Base_Org_Bottom) / max(Base_Org_Bottom, 1e-6)) / range_fraction;
%         S_OrgT(i) = ((Outputs.Org_Top    - Base_Org_Top) / max(Base_Org_Top, 1e-6)) / range_fraction;
        S_ALK(i) = ((Outputs.ALK_Bottom - Base_ALK) / max(Base_ALK, 1e-6)) / range_fraction;
        S_pH(i)  = ((Outputs.pH_Bottom  - Base_pH)  / max(Base_pH,  1e-6)) / range_fraction;
        S_CH4(i) = ((Outputs.CH4_Bottom - Base_CH4) / max(Base_CH4, 1e-6)) / range_fraction;

        S_ALK5(i) = ((Outputs.ALK_Bot5   - B_ALK5) / max(B_ALK5, 1e-6)) / range_fraction;
        S_Sig5(i) = ((Outputs.Sigma_Top5 - B_Sig5) / max(abs(B_Sig5), 1e-4)) / range_fraction;
        S_Ca5(i)  = ((Outputs.CaCO3_Top5 - B_Ca5)  / max(B_Ca5,  1e-6)) / range_fraction;
        S_Meth(i) = ((Outputs.Integ_Meth - B_Meth) / max(B_Meth, 1e-6)) / range_fraction;

        S_OPD(i)  = ((Outputs.OPD       - Base_OPD)  / max(Base_OPD,  0.1)) / range_fraction;
        S_SO4D(i) = ((Outputs.SO4_Depth - Base_SO4D) / max(Base_SO4D, 0.1)) / range_fraction;
%         S_CH4D(i) = ((Outputs.CH4_Onset - Base_CH4D) / max(Base_CH4D, 0.1)) / range_fraction;
%         S_SigD(i) = ((Outputs.Sigma0_Depth- Base_SigD) / max(Base_SigD, 0.1)) / range_fraction;
%         S_CaD(i)  = ((Outputs.CaCO3_Front- Base_CaD)  / max(Base_CaD,  0.1)) / range_fraction;


        
    catch ME
        S_OrgB(i) = NaN; S_ALK(i) = NaN; S_pH(i) = NaN; S_CH4(i) = NaN;
        fprintf('报错信息为: %s\n', ME.message);
    end
end
fprintf('Scan Complete in %.2f seconds.\n', toc);

% --- 4. 可扩展动态绘图模块 ---
% 【扩展指南】未来若要增加输出，只需在这两个 Cell Array 中添加新变量和标题即可
Targets = {S_OrgB, S_ALK, S_pH, S_CH4, S_OPD, S_SO4D, S_Meth, S_Sig5, S_Ca5}; 
Titles  = {'Bottom Organic (%)', 'Bottom ALK (\muM)', 'Bottom pH', 'Bottom CH_4 (\muM)', ...
        'OPD (O_2 < 1\muM)', 'SO_4 Depletion Depth','Integrated Methanogenesis',...
         'Mean \Omega-1 (Top 5cm)', 'Mean CaCO_3 (Top 5cm)'};

n_targets = length(Targets);

% 自动计算子图的行列数
cols = ceil(sqrt(n_targets));
rows = ceil(n_targets / cols);

% 自动调整窗口大小以适应子图数量
figure('Name', 'Multi-Output Sensitivity', 'Color', 'w', 'Position', [100, 100, cols*400, rows*350]);

for k = 1:n_targets
    subplot(rows, cols, k);
%     [sorted_S, idx] = sort(Targets{k}, 'ascend');
%     sorted_Names = Param_Names(idx);
    
    % 使用统一的配色
%      barh(sorted_S, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
barh(Targets{k}, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
set(gca, 'YTick', 1:num_params, 'YTickLabel', Param_Names, 'TickLabelInterpreter', 'none');
set(gca, 'YDir', 'reverse');
%      set(gca, 'YTick', 1:num_params, 'YTickLabel', sorted_Names, 'TickLabelInterpreter', 'none');
%     xlabel('Sensitivity Index (S_{range})');
    title(['Sensitivity of: ', Titles{k}]);
    xline(0, 'k-', 'LineWidth', 1.2);
    grid on;

    max_abs_val = max(abs(Targets{k}(~isnan(Targets{k})))); 
    
    % 2. 安全保护：如果算出全是 0 或报错，给一个默认基准 1
    if isempty(max_abs_val) || max_abs_val == 0
        max_abs_val = 1; 
    end
    
    % 3. 强制锁定坐标轴边界，乘 1.1 是为了左右两边留出 10% 的空白余量好看
    xlim([-max_abs_val * 1.1, max_abs_val * 1.1]);
end