function Plot_Validation_Obs(Result, caseName)
% Plot model profiles with extracted validation observations.
% Usage:
%   Result = Run_RTM_1D_PDE(Config, Params, false);
%   Plot_Validation_Obs(Result, Case.name);
%
% The observation data are embedded below. Depths are in cm.

if nargin < 2 || isempty(caseName)
    if isfield(Result, 'Case') && isfield(Result.Case, 'name')
        caseName = Result.Case.name;
    else
        error('caseName is required unless Result.Case.name exists.');
    end
end

obsSite = map_case_to_obs_site(caseName);
Obs = load_validation_obs();
Obs = Obs(strcmpi(Obs.site, obsSite), :);

if isempty(Obs)
    warning('No validation observations found for case: %s', caseName);
    return
end

S = Result.State_final;
D = Result.Diag_final;
z = Result.Grid.z;

plotOrder = {'OM','POC_mmol_gdw','O2','SO4','Fe2','H2S_mM','CH4','DIC','DIC_mM','ALK','ALK_mM','pH','Meth'};
vars = unique(Obs.var, 'stable');
varsOrdered = {};
for i = 1:numel(plotOrder)
    if any(strcmp(vars, plotOrder{i}))
        varsOrdered{end+1} = plotOrder{i}; %#ok<AGROW>
    end
end
for i = 1:numel(vars)
    if ~any(strcmp(varsOrdered, vars{i}))
        varsOrdered{end+1} = vars{i}; %#ok<AGROW>
    end
end

nvar = numel(varsOrdered);
ncol = min(4, nvar);
nrow = ceil(nvar / ncol);

figure('Color','w','Position',[80 80 320*ncol 330*nrow]);

hLegend = [];

hLegend = [];

for i = 1:nvar
    varName = varsOrdered{i};

    idx = strcmp(Obs.var, varName);
    obsDepth = Obs.depth_cm(idx);
    obsValue = Obs.value(idx);

    [xModel, xLabel, panelTitle] = get_model_profile_for_obs(varName, obsSite, S, D);

    subplot(nrow, ncol, i); hold on;

    hModel = plot(xModel, z, 'k-', 'LineWidth', 2);

    hObs = [];
    if ~isempty(obsValue)
        hObs = scatter(obsValue, obsDepth, 32, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    end

    if isempty(hLegend) && ~isempty(hObs)
        hLegend = [hModel hObs];
    end

    set(gca, 'YDir', 'reverse');
    ylim([0 max(z)]);

    xlabel(xLabel);
    ylabel('Depth (cm)');
    title(panelTitle);
    grid on; box on;

end

sgtitle(strrep(obsSite, '_', '\_'), 'FontWeight', 'bold');

if ~isempty(hLegend)
    lgd = legend(hLegend, {'Model', 'Observation'}, ...
        'Orientation', 'horizontal', ...
        'Box', 'off');
    lgd.Position = [0.42 0.01 0.18 0.03];
end

end

%% ========================================================================
function obsSite = map_case_to_obs_site(caseName)

c = lower(string(caseName));

if contains(c, 'geneva') && contains(c, '1')
    obsSite = 'Geneva1';
elseif contains(c, 'geneva') && contains(c, '4')
    obsSite = 'Geneva 4';
elseif contains(c, 'sap1')
    obsSite = 'Gerogia_SAP1';
elseif contains(c, 'stl2')
    obsSite = 'Gerogia_STL2';
elseif contains(c, 'michigan')
    obsSite = 'Michigan';
else
    obsSite = char(caseName);
end

end

%% ========================================================================
function [x, xLabel, panelTitle] = get_model_profile_for_obs(varName, obsSite, S, D)

switch varName
    case 'OM'
        x = 100 .* (S.OM_lab + S.OM_ref);
        xLabel = 'Organic / TOC (% gDW)';
        panelTitle = 'Organic / TOC';

    case 'POC_mmol_gdw'
        x = (S.OM_lab + S.OM_ref) .* 1000 ./ 12;
        xLabel = 'POC (mmol C gDW^{-1})';
        panelTitle = 'POC';

    case 'O2'
        x = S.O2;
        xLabel = 'O_2 (\muM)';
        panelTitle = 'O_2';

    case 'SO4'
        if strcmpi(obsSite, 'Gerogia_SAP1')
            x = S.SO4 ./ 1000;
            xLabel = 'SO_4 (mM)';
        else
            x = S.SO4;
            xLabel = 'SO_4 (\muM)';
        end
        panelTitle = 'SO_4';

    case 'Fe2'
        x = S.Fe2;
        xLabel = 'Fe^{2+} (\muM)';
        panelTitle = 'Fe^{2+}';

    
    case 'H2S'
        x = S.HS;
        xLabel = 'H_2S / HS (\muM)';
        panelTitle = 'H_2S / HS';

    case 'H2S_mM'
        x = S.HS ./ 1000;
        xLabel = 'H_2S / HS (mM)';
        panelTitle = 'H_2S / HS';

    case 'CH4'
        x = S.CH4;
        xLabel = 'CH_4 (\muM)';
        panelTitle = 'CH_4';

    case 'DIC'
        x = S.DIC;
        xLabel = 'DIC (\muM)';
        panelTitle = 'DIC';

    case 'DIC_mM'
        x = S.DIC ./ 1000;
        xLabel = 'DIC (mM)';
        panelTitle = 'DIC';

    case 'ALK'
        x = S.ALK;
        xLabel = 'ALK (\muM)';
        panelTitle = 'ALK';

    case 'ALK_mM'
        x = S.ALK ./ 1000;
        xLabel = 'ALK (mM)';
        panelTitle = 'ALK';

    case 'pH'
        x = D.pH;
        xLabel = 'pH';
        panelTitle = 'pH';

    case 'Meth'
        x = D.R_Meth ./ 365;
        xLabel = 'Methanogenesis (\muM d^{-1})';
        panelTitle = 'Methanogenesis';

    otherwise
        error('No model mapping for observation variable: %s', varName);
end

end

%% ========================================================================
function Obs = load_validation_obs()

raw = [
        "Geneva1,Fe2,0,40.7138",
        "Geneva1,Fe2,0.373507,46.5522",
        "Geneva1,Fe2,1.36257,50.3599",
        "Geneva1,Fe2,3.23658,49.8522",
        "Geneva1,Fe2,5.2147,60.5136",
        "Geneva1,Fe2,7.29693,59.4983",
        "Geneva1,Fe2,9.17094,55.6906",
        "Geneva1,Fe2,11.3573,61.2752",
        "Geneva1,Fe2,13.0751,62.5444",
        "Geneva1,Fe2,15.4176,61.2752",
        "Geneva1,Fe2,17.0834,68.6367",
        "Geneva1,Fe2,19.2698,70.1597",
        "Geneva1,Fe2,21.0397,72.4443",
        "Geneva1,Fe2,23.1739,73.9674",
        "Geneva1,Fe2,25.1,78.5366",
        "Geneva1,SO4,0,91.6667",
        "Geneva1,SO4,0,36.6667",
        "Geneva1,SO4,0,20",
        "Geneva1,SO4,0,10",
        "Geneva1,SO4,1.30045,20",
        "Geneva1,SO4,3.24365,16.6667",
        "Geneva1,SO4,4.21525,21.6667",
        "Geneva1,SO4,5.26158,6.66667",
        "Geneva1,SO4,6.15845,35",
        "Geneva1,SO4,7.27952,76.6667",
        "Geneva1,SO4,8.25112,11.6667",
        "Geneva1,SO4,10.0448,20",
        "Geneva1,SO4,11.0912,20",
        "Geneva1,SO4,13.3333,21.6667",
        "Geneva1,SO4,14.9776,21.6667",
        "Geneva1,SO4,15.426,21.6667",
        "Geneva1,SO4,17.145,11.6667",
        "Geneva1,SO4,18.1166,15",
        "Geneva1,SO4,19.3871,15",
        "Geneva1,SO4,19.9103,11.6667",
        "Geneva1,SO4,20.8072,10",
        "Geneva1,SO4,22.1525,15",
        "Geneva1,SO4,23.1988,5",
        "Geneva1,SO4,23.8714,5",
        "Geneva1,OM,0.564738,0.731959",
        "Geneva1,OM,1.57025,0.752577",
        "Geneva1,OM,2.54821,0.570447",
        "Geneva1,OM,3.53994,0.453608",
        "Geneva1,OM,4.50413,0.42268",
        "Geneva1,OM,5.53719,0.371134",
        "Geneva1,OM,6.4876,0.446735",
        "Geneva1,OM,7.50689,0.505155",
        "Geneva1,OM,8.5124,0.649485",
        "Geneva1,OM,9.50413,0.549828",
        "Geneva1,OM,10.4821,0.271478",
        "Geneva1,OM,11.4876,0.474227",
        "Geneva1,CH4,0.919811,190.476",
        "Geneva1,CH4,1.91038,285.714",
        "Geneva1,CH4,2.40566,571.429",
        "Geneva1,CH4,3.18396,857.143",
        "Geneva1,CH4,4.03302,1047.62",
        "Geneva1,CH4,5.09434,1809.52",
        "Geneva1,CH4,6.22642,2285.71",
        "Geneva1,CH4,6.86321,2888.89",
        "Geneva1,CH4,7.85377,2984.13",
        "Geneva1,CH4,8.77358,2095.24",
        "Geneva1,CH4,10.0472,4285.71",
        "Geneva1,CH4,12.1698,6000",
        "Geneva1,CH4,12.7358,4952.38",
        "Geneva1,CH4,15.7783,6380.95",
        "Geneva1,CH4,16.9811,6380.95",
        "Geneva1,CH4,17.8302,9365.08",
        "Geneva1,CH4,19.8113,5047.62",
        "Geneva1,CH4,20.8726,6952.38",
        "Geneva1,DIC,0.0707547,4955.84",
        "Geneva1,DIC,1.13208,5735.06",
        "Geneva1,DIC,2.33491,9070.13",
        "Geneva1,DIC,3.18396,9724.68",
        "Geneva1,DIC,3.82075,7324.68",
        "Geneva1,DIC,4.95283,7792.21",
        "Geneva1,DIC,6.15566,6825.97",
        "Geneva1,DIC,6.86321,7044.16",
        "Geneva1,DIC,8.06604,7418.18",
        "Geneva1,DIC,8.91509,10036.4",
        "Geneva1,DIC,9.76415,8789.61",
        "Geneva1,DIC,11.0377,8820.78",
        "Geneva1,DIC,12.1698,9818.18",
        "Geneva1,DIC,12.9481,9724.68",
        "Geneva1,DIC,14.0094,10379.2",
        "Geneva1,DIC,14.9292,7106.49",
        "Geneva1,DIC,15.9198,9288.31",
        "Geneva1,DIC,17.2642,8540.26",
        "Geneva1,DIC,18.1132,10410.4",
        "Geneva1,DIC,19.1745,8883.12",
        "Geneva1,DIC,19.8113,7761.04",
        "Geneva1,DIC,20.8726,9194.81",
        "Geneva1,DIC,21.8632,6670.13",
        "Geneva1,DIC,22.9953,8072.73",
        "Geneva1,DIC,24.0566,9070.13",
        "Geneva1,DIC,25.0472,8976.62",
        "Geneva1,DIC,25.8962,8914.29",
        "Geneva1,DIC,27.0283,11501.3",
        "Geneva1,ALK,0,4534.88",
        "Geneva1,ALK,0.638298,5232.56",
        "Geneva1,ALK,2.2695,8178.29",
        "Geneva1,ALK,3.19149,8643.41",
        "Geneva1,ALK,4.1844,6744.19",
        "Geneva1,ALK,5.10638,7093.02",
        "Geneva1,ALK,6.17021,6046.51",
        "Geneva1,ALK,7.02128,6434.11",
        "Geneva1,ALK,8.08511,6511.63",
        "Geneva1,ALK,9.14894,8992.25",
        "Geneva1,ALK,9.92908,7829.46",
        "Geneva1,ALK,10.6383,7906.98",
        "Geneva1,ALK,12.1277,8837.21",
        "Geneva1,ALK,12.9078,8488.37",
        "Geneva1,ALK,13.8298,9186.05",
        "Geneva1,ALK,14.8936,6279.07",
        "Geneva1,ALK,15.7447,8178.29",
        "Geneva1,ALK,17.0213,7558.14",
        "Geneva1,ALK,17.8723,9224.81",
        "Geneva1,ALK,18.7234,8139.53",
        "Geneva1,ALK,19.7163,6976.74",
        "Geneva1,ALK,20.6383,8255.81",
        "Geneva1,ALK,21.9149,6085.27",
        "Geneva1,ALK,22.9787,7480.62",
        "Geneva1,ALK,23.9716,8139.53",
        "Geneva1,ALK,24.8936,8062.02",
        "Geneva1,ALK,26.0993,7829.46",
        "Geneva1,ALK,27.1631,10232.6",
        "Geneva 4,CH4,0,0",
        "Geneva 4,CH4,0.738714,201.835",
        "Geneva 4,CH4,1.43639,220.183",
        "Geneva 4,CH4,1.9699,311.927",
        "Geneva 4,CH4,2.46238,311.927",
        "Geneva 4,CH4,2.7907,440.367",
        "Geneva 4,CH4,4.14501,256.881",
        "Geneva 4,CH4,4.67852,532.11",
        "Geneva 4,CH4,5.04788,587.156",
        "Geneva 4,CH4,6.03283,550.459",
        "Geneva 4,CH4,7.9617,990.826",
        "Geneva 4,CH4,9.84952,1082.57",
        "Geneva 4,CH4,11.8194,1211.01",
        "Geneva 4,CH4,13.8714,1577.98",
        "Geneva 4,CH4,15.9644,1853.21",
        "Geneva 4,CH4,17.8523,1027.52",
        "Geneva 4,CH4,19.9042,2128.44",
        "Geneva 4,CH4,21.7921,1963.3",
        "Geneva 4,CH4,23.762,2256.88",
        "Geneva 4,CH4,25.6908,3064.22",
        "Geneva 4,CH4,26.8399,2844.04",
        "Geneva 4,OM,0.412322,1.11579",
        "Geneva 4,OM,1.35071,0.803509",
        "Geneva 4,OM,2.38863,0.887719",
        "Geneva 4,OM,3.36967,0.673684",
        "Geneva 4,OM,4.47867,0.719298",
        "Geneva 4,OM,5.45972,0.778947",
        "Geneva 4,OM,6.42654,0.842105",
        "Geneva 4,OM,7.46445,0.694737",
        "Geneva 4,OM,8.40284,0.736842",
        "Geneva 4,OM,9.51185,0.894737",
        "Geneva 4,OM,10.9479,0.547368",
        "Geneva 4,OM,12.9668,0.578947",
        "Geneva 4,DIC,0.843882,2065.12",
        "Geneva 4,DIC,2.02532,3293.02",
        "Geneva 4,DIC,2.8692,3516.28",
        "Geneva 4,DIC,3.79747,2809.3",
        "Geneva 4,DIC,4.68354,3516.28",
        "Geneva 4,DIC,6.20253,3646.51",
        "Geneva 4,DIC,7.08861,3739.53",
        "Geneva 4,DIC,8.35443,3646.51",
        "Geneva 4,DIC,9.11392,3795.35",
        "Geneva 4,DIC,9.87342,3795.35",
        "Geneva 4,DIC,10.6329,4074.42",
        "Geneva 4,DIC,11.3502,3925.58",
        "Geneva 4,DIC,11.9831,3795.35",
        "Geneva 4,DIC,12.5316,4074.42",
        "Geneva 4,DIC,13.1224,4260.47",
        "Geneva 4,DIC,14.0506,4297.67",
        "Geneva 4,DIC,15.0633,4260.47",
        "Geneva 4,DIC,16.0338,4372.09",
        "Geneva 4,DIC,16.8354,4372.09",
        "Geneva 4,DIC,17.9325,4409.3",
        "Geneva 4,DIC,18.5654,4409.3",
        "Geneva 4,DIC,19.3671,4353.49",
        "Geneva 4,DIC,20.2532,4372.09",
        "Geneva 4,DIC,20.8861,4427.91",
        "Geneva 4,DIC,21.73,4632.56",
        "Geneva 4,DIC,22.9958,4855.81",
        "Geneva 4,DIC,23.8819,4911.63",
        "Geneva 4,DIC,24.9367,4874.42",
        "Geneva 4,DIC,25.8228,5823.26",
        "Geneva 4,DIC,27.8481,5079.07",
        "Geneva 4,DIC,28.9451,5469.77",
        "Geneva 4,Fe2,0.746500778,10.04854369",
        "Geneva 4,Fe2,1.866251944,17.03883495",
        "Geneva 4,Fe2,2.846034215,24.90291262",
        "Geneva 4,Fe2,3.965785381,28.54368932",
        "Geneva 4,Fe2,5.085536547,29.41747573",
        "Geneva 4,Fe2,5.925349922,31.60194175",
        "Geneva 4,Fe2,6.99844479,33.34951456",
        "Geneva 4,Fe2,7.884914463,38.44660194",
        "Geneva 4,Fe2,8.864696734,22.7184466",
        "Geneva 4,Fe2,9.9844479,11.7961165",
        "Geneva 4,Fe2,10.91757387,13.10679612",
        "Geneva 4,Fe2,11.89735614,11.94174757",
        "Geneva 4,Fe2,12.92379471,15.29126214",
        "Geneva 4,Fe2,13.90357698,14.85436893",
        "Geneva 4,Fe2,14.88335925,15",
        "Geneva 4,Fe2,15.81648523,15.72815534",
        "Geneva 4,Fe2,16.98289269,16.16504854",
        "Geneva 4,Fe2,17.96267496,15.29126214",
        "Geneva 4,Fe2,18.75583204,13.98058252",
        "Geneva 4,Fe2,19.78227061,15.29126214",
        "Geneva 4,Fe2,20.99533437,15.72815534",
        "Geneva 4,Fe2,21.74183515,16.60194175",
        "Geneva 4,Fe2,22.95489891,20.67961165",
        "Geneva 4,Fe2,24.12130638,19.36893204",
        "Geneva 4,Fe2,25.10108865,20.09708738",
        "Geneva 4,Fe2,27.06065319,22.2815534",
        "Geneva 4,Fe2,27.90046656,22.7184466",
        "Geneva 4,Fe2,29.02021773,11.06796117",
        "Geneva 4,SO4,0.322581,392.893",
        "Geneva 4,SO4,1.24424,307.614",
        "Geneva 4,SO4,2.11982,176.65",
        "Geneva 4,SO4,2.94931,76.1421",
        "Geneva 4,SO4,4.0553,46.7005",
        "Geneva 4,SO4,5.02304,46.7005",
        "Geneva 4,SO4,6.12903,48.731",
        "Geneva 4,SO4,6.95853,51.7766",
        "Geneva 4,SO4,8.06452,51.7766",
        "Geneva 4,SO4,9.17051,60.9137",
        "Geneva 4,SO4,10,68.0203",
        "Geneva 4,SO4,11.2442,83.2487",
        "Geneva 4,SO4,12.212,88.3249",
        "Geneva 4,SO4,13.0415,73.0964",
        "Geneva 4,SO4,14.0092,70.0508",
        "Geneva 4,SO4,15.1152,60.9137",
        "Geneva 4,SO4,16.2212,57.868",
        "Geneva 4,SO4,17.0507,67.0051",
        "Geneva 4,SO4,19.2627,82.2335",
        "Geneva 4,SO4,20.1843,91.3706",
        "Geneva 4,SO4,21.3364,100.508",
        "Geneva 4,SO4,22.1659,106.599",
        "Geneva 4,SO4,23.1336,115.736",
        "Geneva 4,SO4,23.9631,94.4162",
        "Geneva 4,SO4,25.0691,97.4619",
        "Geneva 4,SO4,27.0046,106.599",
        "Geneva 4,SO4,28.0645,100.508",
        "Geneva 4,SO4,29.0783,79.1878",
        "Geneva 4,CH4,0.0411523,129.63",
        "Geneva 4,CH4,0.781893,166.667",
        "Geneva 4,CH4,1.52263,222.222",
        "Geneva 4,CH4,2.22222,240.741",
        "Geneva 4,CH4,2.83951,351.852",
        "Geneva 4,CH4,3.86831,185.185",
        "Geneva 4,CH4,4.81481,555.556",
        "Geneva 4,CH4,6.04938,500",
        "Geneva 4,CH4,7.94239,907.407",
        "Geneva 4,CH4,9.79424,1111.11",
        "Geneva 4,CH4,11.8519,1240.74",
        "Geneva 4,CH4,13.8272,1555.56",
        "Geneva 4,CH4,15.9259,1888.89",
        "Geneva 4,CH4,17.9424,1000",
        "Geneva 4,CH4,20.1235,2074.07",
        "Geneva 4,CH4,21.7695,1944.44",
        "Geneva 4,CH4,24.1152,2277.78",
        "Geneva 4,CH4,26.0905,3074.07",
        "Geneva 4,CH4,26.9136,2851.85",
        "Gerogia_STL2,CH4,0,10.42345277",
        "Gerogia_STL2,CH4,0.879765396,218",
        "Gerogia_STL2,CH4,2.903225806,253",
        "Gerogia_STL2,CH4,4.046920821,301",
        "Gerogia_STL2,CH4,5.894428152,266",
        "Gerogia_STL2,CH4,7.478005865,261",
        "Gerogia_STL2,CH4,8.709677419,229",
        "Gerogia_STL2,CH4,10.38123167,225",
        "Gerogia_STL2,CH4,11.96480938,163",
        "Gerogia_STL2,CH4,13.5483871,167",
        "Gerogia_STL2,CH4,14.78005865,180",
        "Gerogia_STL2,CH4,16.71554252,180",
        "Gerogia_STL2,CH4,18.82697947,175",
        "Gerogia_STL2,CH4,20.93841642,218",
        "Gerogia_STL2,CH4,22.52199413,297",
        "Gerogia_STL2,CH4,24.36950147,315",
        "Gerogia_STL2,CH4,26.48093842,374",
        "Gerogia_STL2,CH4,28.06451613,354",
        "Gerogia_STL2,CH4,30.17595308,354",
        "Gerogia_STL2,CH4,32.55131965,409",
        "Gerogia_STL2,CH4,34.13489736,374",
        "Gerogia_STL2,CH4,36.68621701,418",
        "Gerogia_STL2,CH4,38.62170088,436",
        "Gerogia_STL2,CH4,40.73313783,532",
        "Gerogia_STL2,CH4,42.75659824,530",
        "Gerogia_STL2,CH4,44.60410557,601",
        "Gerogia_STL2,DIC,0.103092784,2630",
        "Gerogia_STL2,DIC,1.142857143,9480",
        "Gerogia_STL2,DIC,2.685714286,12900",
        "Gerogia_STL2,DIC,3.971428571,16200",
        "Gerogia_STL2,DIC,5.942857143,17500",
        "Gerogia_STL2,DIC,7.314285714,20800",
        "Gerogia_STL2,DIC,8.342857143,21000",
        "Gerogia_STL2,DIC,10.14285714,22200",
        "Gerogia_STL2,DIC,11.94285714,21900",
        "Gerogia_STL2,DIC,13.22857143,23600",
        "Gerogia_STL2,DIC,14.77142857,24100",
        "Gerogia_STL2,DIC,16.31428571,26100",
        "Gerogia_STL2,DIC,18.28571429,28400",
        "Gerogia_STL2,DIC,20.17142857,32100",
        "Gerogia_STL2,DIC,22.4,35800",
        "Gerogia_STL2,DIC,24.28571429,39300",
        "Gerogia_STL2,DIC,26.25714286,42400",
        "Gerogia_STL2,DIC,28.05714286,43800",
        "Gerogia_STL2,DIC,30.2,45300",
        "Gerogia_STL2,DIC,32.25714286,49600",
        "Gerogia_STL2,DIC,34.31428571,49300",
        "Gerogia_STL2,DIC,36.28571429,50700",
        "Gerogia_STL2,DIC,38.42857143,51500",
        "Gerogia_STL2,DIC,42.28571429,54300",
        "Gerogia_STL2,DIC,44.51428571,56100",
        "Gerogia_STL2,Fe2,2.879581152,0.483091787",
        "Gerogia_STL2,Fe2,4.45026178,0.483091787",
        "Gerogia_STL2,Fe2,7.068062827,1.449275362",
        "Gerogia_STL2,Fe2,12.04188482,0.483091787",
        "Gerogia_STL2,Fe2,13.87434555,1.449275362",
        "Gerogia_STL2,Fe2,22.68760908,0.483091787",
        "Gerogia_STL2,Fe2,30.36649215,0",
        "Gerogia_STL2,Fe2,36.64921466,2.898550725",
        "Gerogia_STL2,Fe2,38.7434555,1.93236715",
        "Gerogia_STL2,Fe2,40.4886562,1.449275362",
        "Gerogia_STL2,Fe2,42.32111693,3.381642512",
        "Gerogia_STL2,Fe2,44.5026178,1.93236715",
        "Gerogia_STL2,H2S,0.171821306,761",
        "Gerogia_STL2,H2S,1.4604811,5710",
        "Gerogia_STL2,H2S,3.092783505,12400",
        "Gerogia_STL2,H2S,4.639175258,14300",
        "Gerogia_STL2,H2S,5.841924399,15600",
        "Gerogia_STL2,H2S,7.216494845,16000",
        "Gerogia_STL2,H2S,8.934707904,15900",
        "Gerogia_STL2,H2S,10.56701031,15600",
        "Gerogia_STL2,H2S,12.02749141,17400",
        "Gerogia_STL2,H2S,13.40206186,17000",
        "Gerogia_STL2,H2S,15.20618557,18500",
        "Gerogia_STL2,H2S,16.2371134,16400",
        "Gerogia_STL2,H2S,18.29896907,17300",
        "Gerogia_STL2,H2S,20.10309278,15200",
        "Gerogia_STL2,H2S,22.16494845,16500",
        "Gerogia_STL2,H2S,24.14089347,14900",
        "Gerogia_STL2,H2S,26.28865979,14000",
        "Gerogia_STL2,H2S,28.52233677,12900",
        "Gerogia_STL2,H2S,30.32646048,14400",
        "Gerogia_STL2,H2S,32.38831615,13700",
        "Gerogia_STL2,H2S,34.27835052,13900",
        "Gerogia_STL2,H2S,35.99656357,12500",
        "Gerogia_STL2,H2S,38.40206186,13500",
        "Gerogia_STL2,H2S,40.20618557,12600",
        "Gerogia_STL2,H2S,42.43986254,12600",
        "Gerogia_STL2,H2S,44.32989691,11700",
        "Gerogia_STL2,pH,1.038062284,6.961750405",
        "Gerogia_STL2,pH,2.335640138,6.891734198",
        "Gerogia_STL2,pH,4.152249135,6.863209076",
        "Gerogia_STL2,pH,6.487889273,6.88006483",
        "Gerogia_STL2,pH,8.564013841,6.902106969",
        "Gerogia_STL2,pH,10.29411765,6.88006483",
        "Gerogia_STL2,pH,11.5916955,6.870988655",
        "Gerogia_STL2,pH,13.23529412,6.870988655",
        "Gerogia_STL2,pH,14.79238754,6.863209076",
        "Gerogia_STL2,pH,16.34948097,6.883954619",
        "Gerogia_STL2,pH,18.07958478,6.925445705",
        "Gerogia_STL2,pH,20.41522491,6.964343598",
        "Gerogia_STL2,pH,22.3183391,6.991572123",
        "Gerogia_STL2,pH,24.1349481,6.991572123",
        "Gerogia_STL2,pH,26.21107266,7.008427877",
        "Gerogia_STL2,pH,28.5467128,7.018800648",
        "Gerogia_STL2,pH,30.10380623,7.000648298",
        "Gerogia_STL2,pH,31.83391003,6.957860616",
        "Gerogia_STL2,pH,33.9100346,6.953970827",
        "Gerogia_STL2,pH,36.24567474,6.960453809",
        "Gerogia_STL2,pH,38.32179931,6.921555916",
        "Gerogia_STL2,pH,40.39792388,6.872285251",
        "Gerogia_STL2,pH,42.04152249,6.851539708",
        "Gerogia_STL2,pH,44.29065744,6.860615883",
        "Gerogia_STL2,SO4,0,24300",
        "Gerogia_STL2,SO4,1.475694444,20000",
        "Gerogia_STL2,SO4,2.777777778,16700",
        "Gerogia_STL2,SO4,4.166666667,15300",
        "Gerogia_STL2,SO4,6.25,15100",
        "Gerogia_STL2,SO4,7.552083333,13800",
        "Gerogia_STL2,SO4,9.114583333,14500",
        "Gerogia_STL2,SO4,10.41666667,13800",
        "Gerogia_STL2,SO4,12.23958333,13500",
        "Gerogia_STL2,SO4,13.45486111,13400",
        "Gerogia_STL2,SO4,15.10416667,13100",
        "Gerogia_STL2,SO4,16.40625,11800",
        "Gerogia_STL2,SO4,17.96875,9670",
        "Gerogia_STL2,SO4,20.48611111,7180",
        "Gerogia_STL2,SO4,22.39583333,5650",
        "Gerogia_STL2,SO4,24.21875,3440",
        "Gerogia_STL2,SO4,26.5625,2920",
        "Gerogia_STL2,SO4,28.55902778,1720",
        "Gerogia_STL2,SO4,30.46875,766",
        "Gerogia_STL2,SO4,32.46527778,718",
        "Gerogia_STL2,SO4,34.28819444,478",
        "Gerogia_STL2,SO4,35.67708333,431",
        "Gerogia_STL2,SO4,36.97916667,431",
        "Gerogia_STL2,SO4,38.28125,574",
        "Gerogia_STL2,SO4,39.58333333,287",
        "Gerogia_STL2,SO4,42.96875,191",
        "Gerogia_STL2,SO4,44.44444444,287",
        "Gerogia_SAP1,CH4,0,10.4235",
        "Gerogia_SAP1,CH4,0.879765,217.59",
        "Gerogia_SAP1,CH4,2.90323,252.769",
        "Gerogia_SAP1,CH4,4.04692,300.977",
        "Gerogia_SAP1,CH4,5.89443,265.798",
        "Gerogia_SAP1,CH4,7.47801,260.586",
        "Gerogia_SAP1,CH4,8.70968,229.316",
        "Gerogia_SAP1,CH4,10.3812,225.407",
        "Gerogia_SAP1,CH4,11.9648,162.866",
        "Gerogia_SAP1,CH4,13.5484,166.775",
        "Gerogia_SAP1,CH4,14.7801,179.805",
        "Gerogia_SAP1,CH4,16.7155,179.805",
        "Gerogia_SAP1,CH4,18.827,174.593",
        "Gerogia_SAP1,CH4,20.9384,217.59",
        "Gerogia_SAP1,CH4,22.522,297.068",
        "Gerogia_SAP1,CH4,24.3695,315.309",
        "Gerogia_SAP1,CH4,26.4809,373.941",
        "Gerogia_SAP1,CH4,28.0645,354.397",
        "Gerogia_SAP1,CH4,30.176,354.397",
        "Gerogia_SAP1,CH4,32.5513,409.121",
        "Gerogia_SAP1,CH4,34.1349,373.941",
        "Gerogia_SAP1,CH4,36.6862,418.241",
        "Gerogia_SAP1,CH4,38.6217,436.482",
        "Gerogia_SAP1,CH4,40.7331,531.596",
        "Gerogia_SAP1,CH4,42.7566,530.293",
        "Gerogia_SAP1,CH4,44.6041,600.651",
        "Michigan,DIC,0,2540",
        "Michigan,DIC,0.450625869,2730",
        "Michigan,DIC,1.502086231,3090",
        "Michigan,DIC,2.478442281,3080",
        "Michigan,DIC,3.454798331,2800",
        "Michigan,DIC,4.406119611,2800",
        "Michigan,DIC,5.482614743,2800",
        "Michigan,DIC,6.458970793,2890",
        "Michigan,DIC,7.485396384,2800",
        "Michigan,DIC,8.411682893,2800",
        "Michigan,DIC,9.513212796,2900",
        "Michigan,DIC,10.4394993,2880",
        "Michigan,DIC,11.34075104,2870",
        "Michigan,DIC,12.54242003,2760",
        "Michigan,DIC,13.44367177,2820",
        "Michigan,DIC,14.49513213,2740",
        "Michigan,DIC,15.47148818,2850",
        "Michigan,DIC,16.52294854,2820",
        "Michigan,DIC,17.39916551,2840",
        "Michigan,Fe2,0.531818,19.8212",
        "Michigan,Fe2,1.67727,8.26583",
        "Michigan,Fe2,2.29091,13.8561",
        "Michigan,Fe2,3.27273,15.1228",
        "Michigan,Fe2,4.58182,16.3378",
        "Michigan,Fe2,5.48182,11.9302",
        "Michigan,Fe2,6.38182,17.4752",
        "Michigan,Fe2,7.56818,20.1314",
        "Michigan,Fe2,8.50909,29.9354",
        "Michigan,Fe2,9.49091,49.6855",
        "Michigan,Fe2,10.6364,16.8031",
        "Michigan,Fe2,11.5364,32.3007",
        "Michigan,Fe2,12.4773,16.5123",
        "Michigan,Fe2,13.3773,32.0099",
        "Michigan,Fe2,14.6045,19.0198",
        "Michigan,Fe2,15.3409,55.8703",
        "Michigan,Fe2,16.2409,62.8371",
        "Michigan,Fe2,17.5909,73.9983",
        "Michigan,SO4,0,181.69",
        "Michigan,SO4,0.47191,207.042",
        "Michigan,SO4,1.48315,298.592",
        "Michigan,SO4,2.62921,230.282",
        "Michigan,SO4,3.50562,175.352",
        "Michigan,SO4,4.44944,217.606",
        "Michigan,SO4,5.57303,143.662",
        "Michigan,SO4,6.65169,116.901",
        "Michigan,SO4,7.55056,112.676",
        "Michigan,SO4,8.5618,45.0704",
        "Michigan,SO4,9.61798,11.2676",
        "Michigan,SO4,10.5843,25.3521",
        "Michigan,SO4,11.5056,8.4507",
        "Michigan,SO4,12.4045,2.11268",
        "Michigan,SO4,13.4607,16.9014",
        "Michigan,SO4,14.4719,14.7887",
        "Michigan,SO4,15.4831,12.6761",
        "Michigan,SO4,16.4944,10.5634",
        "Michigan,SO4,17.5056,10.5634",
        "Michigan,POC_mmol_gdw,0.54,4.04",
        "Michigan,POC_mmol_gdw,1.58394,3.98072",
        "Michigan,POC_mmol_gdw,2.65925,3.78204",
        "Michigan,POC_mmol_gdw,3.50302,3.84608",
        "Michigan,POC_mmol_gdw,4.40924,3.61492",
        "Michigan,POC_mmol_gdw,5.68913,3.55521",
        "Michigan,POC_mmol_gdw,6.46439,3.5538",
        "Michigan,POC_mmol_gdw,7.74069,3.28098",
        "Michigan,POC_mmol_gdw,8.57784,2.95159",
        "Michigan,POC_mmol_gdw,9.65369,2.78569",
        "Michigan,POC_mmol_gdw,10.4276,2.70232",
        "Michigan,POC_mmol_gdw,11.4373,2.61031",
        "Michigan,POC_mmol_gdw,12.6858,2.69001",
        "Michigan,POC_mmol_gdw,13.6275,2.56534",
        "Michigan,POC_mmol_gdw,14.5029,2.50637",
        "Michigan,POC_mmol_gdw,15.5818,2.5208",
        "Michigan,POC_mmol_gdw,16.4253,2.56844",
        "Michigan,POC_mmol_gdw,17.5356,2.44347",
        "Michigan,Meth,0.572614,1.36054",
        "Michigan,Meth,1.46888,0.340136",
        "Michigan,Meth,2.51452,1.02041",
        "Michigan,Meth,3.51037,1.36054",
        "Michigan,Meth,4.55602,1.02041",
        "Michigan,Meth,5.45228,1.02041",
        "Michigan,Meth,6.42324,1.02041",
        "Michigan,Meth,7.46888,1.02041",
        "Michigan,Meth,8.43983,1.02041",
        "Michigan,Meth,9.48548,3.06122",
        "Michigan,Meth,10.4564,1.36054",
        "Michigan,Meth,11.5021,1.36054",
        "Michigan,Meth,12.5228,1.36054",
        "Michigan,Meth,13.5187,1.02041",
        "Michigan,Meth,14.4647,1.02041",
        "Michigan,Meth,15.5104,1.02041",
        "Michigan,Meth,16.5062,1.02041",
        "Michigan,Meth,17.527,1.02041"
    ];

site = cell(numel(raw), 1);
var = cell(numel(raw), 1);
depth_cm = nan(numel(raw), 1);
value = nan(numel(raw), 1);

for i = 1:numel(raw)
    parts = split(raw(i), ',');
    site{i} = char(parts(1));
    var{i} = char(parts(2));
    depth_cm(i) = str2double(parts(3));
    value(i) = str2double(parts(4));
end

Obs = table(site, var, depth_cm, value);

end
