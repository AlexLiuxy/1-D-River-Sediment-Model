function Config = Config()
    % CONFIG_USGS_BASELINE: Site-specific parameters mapped from USGS data.
    % These will be randomized via Latin Hypercube Sampling in Step 4.
    
    % -- Domain & Hydrology --
    Config.Lbottom       = 30;       % Max depth [cm]
    Config.n_grid        = 101;      % Number of spatial nodes
    Config.vbottom       = 0.5;      % Sedimentation rate [cm/yr] (AUDITED: 100 -> 0.5)
    Config.vbottom_fluid = 0;        % Groundwater upwelling [cm/yr]
    
    % -- Sediment Physical Properties --
    Config.porostop      = 0.9;      % Porosity at SWI
    Config.porosbottom   = 0.7;      % Porosity at depth
    Config.porosscale    = 3;        % e-folding depth for porosity [cm]
    Config.Bioturbtop    = 10;       % SWI Bioturbation [cm2/yr]
    Config.Bioturbbottom = 1;        % Deep Bioturbation [cm2/yr]
    Config.bioturbscale  = 3;        % e-folding depth for bioturbation [cm]
    Config.Bioirrig_top  = 100;      % SWI Bioirrigation [1/yr]
    Config.Bioirrig_scale= 0.75;     % e-folding depth for bioirrigation [cm]
    
    % -- Top Boundary Conditions (SWI Concentrations) --
    Config.O2init   = 150;           % Oxygen [uM]
    Config.SO4init  = 200;           % Sulfate [uM] (AUDITED: 3000 -> 200)
    Config.DICinit  = 1000;          % DIC [uM]
    Config.HCO3init = 950;           % Alkalinity [uM]
    Config.Calcium  = 1000;          % Calcium [uM]
    Config.CH4init  = 0;             % Methane [uM]
    Config.Feinit   = 0;             % Dissolved Fe [uM]
    Config.HSinit   = 0;             % Sulfide [uM]
    
    % -- Depositional Fluxes --
    Config.NPP      = 600;           % Net Primary Production [g/m2/yr]
    Config.BE       = 0.2;           % Burial Efficiency
    Config.F_FeOx   = 10;            % FeOx flux [mmol/m2/d]
    Config.F_CaCO3  = 10;            % CaCO3 flux [g/m2/yr]
    
    % -- Climate / Environmental --
    Config.T_future = 25;            % Local water temperature [C]
    Config.ageinit  = 0.1;           % Initial age of OM at SWI [yr]
end