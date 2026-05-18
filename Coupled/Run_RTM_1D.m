%% Sediment Diagenesis Model (ODEs for organic matter, oxygen, sulfate, DIC, ALK, Carb Acid, and pH)

% -------- Process included in the model ----------------------------------
% Organic matter: Organic matter degradation using Middelburg 
% Oxygen: Oxic respiration 
% Sulfate: Anoxic respiration through sulfate reduction.
% Dissolved inorganic carbon: Organic matter degradation and carbonate precipitation
% Alkalinity: Sulfate reduction and carbonate precipitation

% function Outputs = Run_RTM_1D(Custom_Config, Custom_Params)
% RUN_RTM_1D
% Sediment Diagenesis Model - Functionized for Sensitivity Analysis
clear all
% ----------------------------- INPUT PARAMETERS ---------------------------

global v_burial Mineral_Mass z_sed Oxygen Sulfate Corg_top
global k_sed k_O2 DSO4 DH2S DO2 DPO4 k_SO4 Kreox  Bioturb Calcium DHCO3 HCO3init 
global O2init SO4init HSinit C_organic rho poros RC Alpha_Bioirrig
global R_respi R_SRR R_FeRed RC_after_Fe R_DIC_prod R_ALK_prod R_AOM R_CH4Ox Ksp_ca k_calcite DICinit R1_carb CO3_1 BE P_C_ratio Rviv1 R_FeS  R_FeOx Fe_3_init
global v_burial_Fluid CO3_activity Calcium_activity NPP kFeS FeooH Feinit R_HS_Ox kapatite
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2 CH4init Pinitial DCH4 kFeOx KFEMonod Sulfide Rapat CaCO3 F_CaCO3 O2_root
global C_HS C_Fe n_power_CaCO31 n_power_CaCO32 k_calcite_dis1 n_power_CaCO33 k_calcite_dis2 CaCO3_init Temp_factor 
global T_future Rate_Meth Salinity pH K_HS R_AOM_lag R_AOM_actual R_AOM_pot SO4_diag F_FeOx 
global F_lab_OM F_ref_OM C_organic_lab C_organic_ref C_organic_total k_ref_factor Fe_inventory_factor
%global KFe_HS Iron_conc R_iron Iron_C P_apaeq R1_carb_disso R1_carb_form
    

%     if nargin < 1 || isempty(Custom_Config)
%         Config = Config_Baseline();
%     else
%         Config = Custom_Config;
%     end
%     if nargin < 2 || isempty(Custom_Params)
%         Params = Params_Static();
%     else
%         Params = Custom_Params;
%     end
        Params = Params_Static();
    Config = Config_Baseline();
%     Hydro  = Hydro_Preprocessor(Config, Params);

    rho = Params.rho;
    Mineral_Mass = 215;   % keep as legacy until explicitly audited

    k_O2 = Params.k_O2;
    k_SO4 = Params.k_SO4;
    KFEMonod = Params.KFEMonod;

    DSO4 = Params.DSO4;
    DCH4 = Params.DCH4;
    DH2S = Params.DH2S;
    DO2  = Params.DO2;
    DHCO3 = Params.DHCO3;
    DPO4  = Params.DPO4;

    Kreox = Params.Kreox;
    kFeOx = Params.kFeOx;
    kFeS  = Params.kFeS;
    FeC_frac_max = Params.FeC_frac_max;
    Fe_inventory_factor = Params.Fe_inventory_factor;

    K_CH4_SO4 = Params.K_CH4_SO4;
    K_CH4_O2  = Params.K_CH4_O2;
    k_AOM = Params.k_AOM;
    k_aerobic_CH4 = Params.k_aerobic_CH4;

    Ksp_ca = Params.Ksp_ca;
    k_calcite = Params.k_calcite;
    k_calcite_dis1 = Params.k_calcite_dis1;
    k_calcite_dis2 = Params.k_calcite_dis2;
    n_power_CaCO31 = Params.n_power_CaCO31;
    n_power_CaCO32 = Params.n_power_CaCO32;
    n_power_CaCO33 = Params.n_power_CaCO33;
    Calcium_activity = Params.Calcium_activity;
    CO3_activity     = Params.CO3_activity;
    P_C_ratio = Params.P_C_ratio;

    kapatite = Params.kapatite;
%     KFeS = Params.KFeS;
    K_HS = Params.K_HS;
    Q10 = Params.Q10;
    T_ref = Params.T_ref;

    % Site inputs
    n = Config.n;
    Bioirrig_top    = Config.Bioirrig_top;
    Bioirrig_bottom = Config.Bioirrig_bottom;
    Bioirrig_scale  = Config.Bioirrig_scale;
    
    Lbottom         = Config.Lbottom;
    
    Bioturbtop      = Config.Bioturbtop;
    Bioturbbottom   = Config.Bioturbbottom;
    bioturbscale    = Config.bioturbscale;
    
    vbottom         = Config.vbottom;
    vbottom_fluid   = Config.vbottom_fluid;
    
    porostop        = Config.porostop;
    porosbottom     = Config.porosbottom;
    porosscale      = Config.porosscale;
    
    O2init          = Config.O2init;
    SO4init         = Config.SO4init;
    DICinit         = Config.DICinit;
    HCO3init        = Config.HCO3init;
    Calcium         = Config.Calcium;
    CH4init         = Config.CH4init;
    Feinit          = Config.Feinit;
    HSinit          = Config.HSinit;
    Pinitial        = Config.Pinitial;
    



    % Corg_top = Config.Corg_top;
    BE = Config.BE;
    NPP = Config.NPP;
    
    f_lab = Config.f_lab;
    k_ref_factor = Config.k_ref_factor;
    
    f_lab = max(0, min(1, f_lab));
    
    F_OM_total = BE * NPP * 1E-4;   % current transitional OM input
    F_lab_OM   = f_lab * F_OM_total;
    F_ref_OM   = (1 - f_lab) * F_OM_total;

    if Config.use_hydro_npp_multiplier
        NPP = NPP * Hydro.NPP_multiplier;
    end

    F_FeOx = Config.F_FeOx;
    F_CaCO3 = Config.F_CaCO3;

    T_future = Config.T_future;

    Salinity = Config.Salinity;
    ageinit  = Config.ageinit;
    age_root = Config.age_root;

    % Optional weak hydro injection to diffusion coefficients
    if Config.use_hydro_diffusion_multiplier
        DSO4 = DSO4 * Hydro.diffusion_multiplier;
        DCH4 = DCH4 * Hydro.diffusion_multiplier;
        DH2S = DH2S * Hydro.diffusion_multiplier;
        DO2  = DO2  * Hydro.diffusion_multiplier;
        DHCO3 = DHCO3 * Hydro.diffusion_multiplier;
    end

    DOC_root_1      = Config.DOC_root_1;
    O2_root_1       = Config.O2_root_1;
    POC_root_1      = Config.POC_root_1;
% --------------- Calculating initial depth profiles for input parameters ----------------

% k_sed = k_sed / 100;
% k_AOM = k_AOM / 10;
% k_aerobic_CH4 = k_aerobic_CH4 / 5;
% k_calcite = k_calcite / 10;
% kFeOx = kFeOx / 5;

    z_sed = linspace(0, Lbottom, n);
    z_biodiff = linspace(0, Lbottom, 10001);
    dz_sed = Lbottom / (n - 1);

    if Config.use_constant_porosity
        poros = Config.constant_porosity .* ones(1,n);
    elseif Config.use_hydro_phi
        poros = Hydro.phi_mix .* ones(1,n);
    else
        poros = Config.porosbottom + (Config.porostop - Config.porosbottom) .* exp(-z_sed / Config.porosscale);
    end

    Bioturb_1 = Config.Bioturbbottom + (Config.Bioturbtop - Config.Bioturbbottom) .* exp(-z_biodiff / Config.bioturbscale);
    Bioturb = interp1(z_biodiff, Bioturb_1, z_sed);

    Alpha_Bioirrig = Config.Bioirrig_bottom + (Config.Bioirrig_top - Config.Bioirrig_bottom) .* exp(-z_sed / Config.Bioirrig_scale);

    v_burial = Config.vbottom .* (1 - Config.porosbottom) ./ (1 - poros);
    v_burial_Fluid = Config.vbottom_fluid .* (1 + Config.porosbottom) ./ (1 + poros);

    age = Config.ageinit + cumsum(dz_sed ./ v_burial);
    k_sed = 10.^(-0.95 .* log10(age) - 0.81);
    Temp_factor = Q10.^((T_future - T_ref) / 10);


% ----------------------- Initial Carbonate concentration -----------------

[pH_top, CO3_top, ~] = River_Carbonate(HCO3init, DICinit, T_future, Salinity, 1); 
% CO3_top = Carb_CO3(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
CO3_1 = CO3_top*ones(1,n);
C_HS  = zeros(1,n); %intial value for sulfide
R_HS_Ox = zeros(1,n);
CaCO3 = 1E5.*zeros(1,n);
Rapat = zeros(1,n);
Rviv1 = zeros(1,n); %umol/Lsed/yr
R_FeS = ones(1,n); %umol/Lsed/yr
pH    = pH_top.*ones(1,n);

R_AOM_lag = zeros(1,n);
R_AOM_pot    = zeros(1,n);   % potential AOM from previous CH4 profile, before SO4 limitation
R_AOM_actual = zeros(1,n);   % actual SO4-supported AOM used by sulfur/carbonate ledger
R_FeOx = zeros(1,n);   % lagged Fe(II) reoxidation used to cap Fe reduction supply

% ----------------------- Correcting organic matter reactivity based on oxygen penetration depth ----------------
% After calculating the oxygen penetration depth, reactivity profiles would
% be corrected using oxic and anoxic power law by Katsev & Crowe (2015).

% ------------- ORGANIC MATTER DEGRADATION --------------------------------

hold on

% Solving ODE

if Bioturbtop == 0

        x = linspace(0,Lbottom,n);
        CorgInit  = (BE * NPP * 1E-4)./(v_burial(1) * rho * (1-poros(1)));
        C_organic = CorgInit*exp(-cumsum(k_sed./v_burial.*dz_sed));
        BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic
    
else 
    
    nmesh=1000;
    x=linspace(0,Lbottom,nmesh);
    solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
    sol = bvp4c(@organicODE,@organicbc,solinit);
    x = linspace(0,Lbottom,n);
    y = deval(sol,x);
    
    if min(y) < 0
        fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
    end
    y = max(y, 1e-12);
    
    C_organic = y(1,:);
    C_organic_lab = C_organic;
    
    solid_flux = rho .* max(1 - poros, 1e-6) .* max(v_burial, 1e-12);
    C_organic_ref = F_ref_OM ./ solid_flux;
    
    C_organic_total = C_organic_lab + C_organic_ref;
    
    C_organic = C_organic_total;  % only for plotting/output compatibility
    BEsed_org = C_organic_total ./ max(C_organic_total(1), 1e-12);

end
% if Bioturbtop == 0
%     x = linspace(0,Lbottom,n);
%     C_organic = Corg_top .* exp(-cumsum((Temp_factor .* k_sed) ./ v_burial .* dz_sed));
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% else
%     nmesh = 1000;
%     x = linspace(0,Lbottom,nmesh);
%     solinit = bvpinit(linspace(0,Lbottom,nmesh), [Corg_top 0]);
%     sol = bvp4c(@organicODE, @organicbc, solinit);
%     x = linspace(0,Lbottom,n);
%     y = deval(sol,x);
%     if min(y) < 0
%         fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
%     end
%     y = max(y, 1e-12);
%     C_organic = y(1,:);
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% end
% ---------------------------- OXYGEN -------------------------------------

RC_lab = Temp_factor .* k_sed .* C_organic_lab .* rho .* ((1-poros)./12);

RC = RC_lab; % molCorg/cm3/yr mineralization rate

O2_root = zeros(1,n);

% Solving ODE
nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@O2_ODE,@O2_bc,solinit);
x = linspace(0,Lbottom,n);
y = deval(sol,x);

if min(y) < 0
    fprintf('Initial O2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_O2 = y(1,:);
Oxygen = C_O2;

% ---------------------------- OXYGEN PENTRATION DEPTS_ --------------------

count_OPD = 0;
for i=1:n
    if Oxygen(1,i) < 1

        count_OPD = count_OPD + 1;
        OPD_1(1,count_OPD) = z_sed(1,i);
        num_OPD1(1,count_OPD) = i;

    end
end

if min(Oxygen) > 1

    OPD_1 = Lbottom;
    num_OPD1 = n;
end

OPD = min(OPD_1);
num_OPD = min(num_OPD1);

mm_count = 0;
for i=1:n

    FeOxInit_raw = (F_FeOx .* 36.5) ./ ...
    (v_burial(1) .* rho .* max(1 - poros(1), 1e-6));

    FeOxInit = Fe_inventory_factor .* FeOxInit_raw;

    Ironoxy(1,i) = 0.01.*FeOxInit.*((1/(1+exp(z_sed(1,i)-OPD)))+2*exp(-((z_sed(1,i)-OPD)^2)/2));
    mm_count=mm_count+1;
    FeOx(1,mm_count)=Ironoxy(1,i);

end
% FeooH =  FeOx;wor
% FeooH =  zeros(1,n);
FeooH = max(FeOx, 1e-12);
% ---------------------------- CORRECTING ACTIVITY PROFILES ---------------
b_oxic = 0.95;%0.977;
a_oxic = 0.81;%0.312;


b_anoxic = 0.95;%0.857;
a_anoxic = 0.81;%1.1;

% % slower reactivity decay below OPD
% b_anoxic = 0.75;
% a_anoxic = 0.60;
% k_sed(1,1:num_OPD) = 10.^(-b_oxic*log10(age(1,1:num_OPD)) - a_oxic);  % Oxic
% k_sed(1,(num_OPD+1):n) = 10.^(-b_anoxic*log10(age(1,(num_OPD+1):n)) - a_anoxic);  % Anoxic
k_oxic   = 10.^(-b_oxic  .* log10(age) - a_oxic);
k_anoxic = 10.^(-b_anoxic .* log10(age) - a_anoxic);

transition_width = 0.5;   % cm
w_oxic = 1 ./ (1 + exp((z_sed - OPD) ./ transition_width));

k_sed = w_oxic .* k_oxic + (1 - w_oxic) .* k_anoxic;

% ------------------- Blue Carbon ------------------------------------

% DOC release 

depth_rootzone = 10; % seagrass root length (cm)
z_root = 0:0.1:depth_rootzone;
mu_root = 4; % value for the center of rootzone in normal distribution
sigma_root = 6; % sigma for normal distribution of flux in the rootzone
DOC_root_2 = DOC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/year
DOC_root_22 = interp1(z_root,DOC_root_2,z_sed);
% poros_root = interp1(z_sed,poros,z_sed);

for i_root = 1:n

    if depth_rootzone >= z_sed(1,i_root) 

DOC_root(1,i_root) = 1E6 * (DOC_root_22(1,i_root)./(z_sed(1,2)-z_sed(1,1))); % (umol/l/year)

    else 

DOC_root(1,i_root) = 0;

    end

end



% O2 release 

muO2_root = 4; % value for the center of rootzone in normal distribution
sigmaO2_root = 6; % sigma for normal distribution of flux in the rootzone
O2_root_2 = O2_root_1 * 1E-4 * normpdf(z_root,muO2_root,sigmaO2_root);  % mmol/cm2/year
O2_root_22 = interp1(z_root,O2_root_2,z_sed);

for i_root = 1:n

    if depth_rootzone >= z_sed(1,i_root) 

O2_root(1,i_root) = 1E6 * (O2_root_22(1,i_root)./(z_sed(1,2)-z_sed(1,1))); % (umol/l/year)

    else 

O2_root(1,i_root) = 0;

    end

end


% Seagrass POC release 

muPOC_root = 4; % value for the center of rootzone in normal distribution
sigmaPOC_root = 6; % sigma for normal distribution of flux in the rootzone

POC_root_2 = POC_root_1 * normpdf(z_root,muPOC_root,sigmaPOC_root);  % mmol/cm2/year
POC_root = interp1(z_root,POC_root_2,z_sed);

k_sed_root = 10.^(-0.95*log10(age_root) - 0.8); % more reactive
poros2 = interp1(z_sed,poros,z_sed);


for i_root = 1:n

    if depth_rootzone >= z_sed(1,i_root) 

% RC_root(1,i_root) = Temp_factor.*k_sed_root .*POC_root(1,i_root).*rho.*((1-poros2(1,i_root))./(12)); % (umol/l/year)
RC_root(1,i_root) = 1.*k_sed_root .*POC_root(1,i_root).*rho.*((1-poros2(1,i_root))./(12)); % (umol/l/year)

    else 

RC_root(1,i_root) = 0;

    end

end

% -------------------------------------------------------------------------

iteration = 1;
count_loop = 1;

max_outer_iter = 12;
min_outer_iter = 5;
conv_tol = 0.05;     % 5% profile/rate convergence

C_Fe_prev   = [];
FeooH_prev  = [];
Sulfate_prev = [];
CH4_prev    = [];
ALK_prev    = [];
pH_prev     = [];
C_HS_prev = [];
R_FeS_prev = [];
while iteration <= max_outer_iter
    

% ------------- ORGANIC MATTER DEGRADATION --------------------------------

hold on


% Solving ODE

if Bioturbtop == 0

x = linspace(0,Lbottom,n);
CorgInit  = (BE * NPP * 1E-4)./(v_burial(1) * rho * (1-poros(1)));
C_organic = CorgInit.*exp(-cumsum((Temp_factor.*k_sed)./v_burial.*dz_sed));
BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic

else 

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@organicODE,@organicbc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);


C_organic = y(1,:);

C_organic_lab = C_organic;

solid_flux = rho .* max(1 - poros, 1e-6) .* max(v_burial, 1e-12);
C_organic_ref = F_ref_OM ./ solid_flux;

C_organic_total = C_organic_lab + C_organic_ref;

C_organic = C_organic_total;  % only for plotting/output compatibility
BEsed_org = C_organic_total ./ max(C_organic_total(1), 1e-12);


end
% if Bioturbtop == 0
%     x = linspace(0,Lbottom,n);
%     C_organic = Corg_top .* exp(-cumsum((Temp_factor .* k_sed) ./ v_burial .* dz_sed));
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% else
%     nmesh = 1000;
%     x = linspace(0,Lbottom,nmesh);
%     solinit = bvpinit(linspace(0,Lbottom,nmesh), [Corg_top 0]);
%     sol = bvp4c(@organicODE, @organicbc, solinit);
%     x = linspace(0,Lbottom,n);
%     y = deval(sol,x);
%     if min(y) < 0
%         fprintf('Initial Organic is negative in the current iteration! Minimum：%.2e\n', min(y));
%     end
%     y = max(y, 1e-12);
%     C_organic = y(1,:);
%     BEsed_org = C_organic ./ max(C_organic(1), 1e-12);
% end
% ---------------------------- OXYGEN -------------------------------------

RC_lab = Temp_factor .* k_sed .* C_organic_lab .* rho .* ((1-poros)./12);


RC = RC_lab + RC_root; % molCorg/cm3/yr mineralization rate

% Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@O2_ODE,@O2_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('O2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);


C_O2 = y(1,:);
Oxygen = C_O2;
% -------- update OPD-based k_sed every iteration --------
OPD_1 = [];
num_OPD1 = [];

for i = 1:n
    if Oxygen(i) < 1
        OPD_1(end+1) = z_sed(i);
        num_OPD1(end+1) = i;
    end
end

if isempty(OPD_1)
    OPD = Lbottom;
    num_OPD = n;
else
    OPD = OPD_1(1);
    num_OPD = num_OPD1(1);
end

R_respi = RC .* (Oxygen ./ (Oxygen + k_O2)) .* 1E9;   % umol/L/yr

% residual-carbon cascade: O2 first, then Fe, then SO4, then CH4
RC_total_uM = RC .* 1E9;                              % umol C / L / yr
RC_after_O2 = max(RC_total_uM - R_respi, 0);

% FeooH_old = FeooH;
% % cap Fe reduction so Fe branch does not consume nearly all residual carbon
% Fe_gate = FeooH_old ./ (FeooH_old + KFEMonod);
% C_to_Fe = min(FeC_frac_max .* RC_after_O2, RC_after_O2 .* Fe_gate);
% R_FeRed = 4 .* C_to_Fe;

FeooH_old = FeooH;

% Potential Fe reduction demand from carbon and local FeOOH availability.
% FeC_frac_max is no longer active.
Fe_gate = FeooH_old ./ max(FeooH_old + KFEMonod, 1e-12);
C_to_Fe_pot = RC_after_O2 .* Fe_gate;
R_FeRed_pot = 4 .* C_to_Fe_pot;

I_FeRed_pot = trapz(z_sed, R_FeRed_pot) .* 1e-3;

solid1 = rho .* max(1 - poros(1), 1e-6);

J_Fe3_burial_in = solid1 .* v_burial(1) .* FeooH_old(1);

J_Fe3_mix_in = -solid1 .* Bioturb(1) .* ...
    ((FeooH_old(2) - FeooH_old(1)) ./ dz_sed);

I_Fe_transport_supply = max(0, J_Fe3_burial_in + J_Fe3_mix_in);

I_Fe_recycle_lag = trapz(z_sed, max(R_FeOx,0)) .* 1e-3;

I_Fe_supply_ext = F_FeOx .* 36.5;  % mmol/m2/d -> umol/cm2/yr

Fe_recycle_cap_frac = 0.25;

I_Fe_recycle_eff = min(I_Fe_recycle_lag, ...
                       Fe_recycle_cap_frac * I_Fe_supply_ext);

% I_Fe_supply_cap = I_Fe_supply_ext + I_Fe_recycle_eff;

I_Fe_supply_cap = I_Fe_supply_ext ;%...
                % + I_Fe_transport_supply ...
                % + I_Fe_recycle_lag;

Fe_supply_scale = min(1, I_Fe_supply_cap ./ max(I_FeRed_pot, 1e-12));

R_FeRed = R_FeRed_pot .* Fe_supply_scale;
C_to_Fe = R_FeRed ./ 4;

C_HS_used_for_Fe2 = C_HS;
pH_used_for_Fe2 = pH;
% ------------------------ IRON(II) ---------------------------------------
 
% Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@Fe_ODE,@Fe_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('Fe2 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_Fe = y(1,:);
% Iron_C(iteration,:) = C_Fe;
% 
% for i=1:n
%   if Iron_C(iteration,i) < 0
%       Iron_C(iteration,i) = 0;
%   end 
% end

% F_diff_Fe(1,count_loop) = DH2S.*((C_Fe(1,2) - C_Fe(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr

% ------------------------ IRON(III) ---------------------------------------

% Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
% R_iron(count_loop,:) = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9; %rate of iron reduction umol/l/year
% R_FeOx = (kFeOx.*C_Fe.*Oxygen);
% % R_FeOx_1(count_loop,:) = R_FeOx;
% % Fe_3_init  = 365.*1E2.*(F_FeOx)./(v_burial(1));  %umol/l
% Fe_3_init = 36.5.*(F_FeOx).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %umol/l
% 
% % KFEMonod = 2000;
% 
% % % Solving ODE
% 
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[Fe_3_init 0]);
% sol = bvp4c(@Fe3_ODE,@Fe3_bc,solinit);
% 
% x = linspace(0,Lbottom,n);
% 
% y = deval(sol,x);
% 
% if min(y) < 0
%     fprintf('Fe3：%.2e\n', min(y));
% end
% y = max(y, 1e-12);
% 
% C_Fe_3 = y(1,:);
% FeooH = C_Fe_3;
% % FeooH = max(C_Fe_3, 0.2 .* FeOx);
% RC_after_Fe = max(RC_after_O2 - R_FeRed ./ 4, 0);   % umol C / L / yr

R_FeOx = kFeOx .* C_Fe .* Oxygen;
Fe_3_init_raw = 36.5 .* F_FeOx .* (poros(1) ./ max(1 - poros(1), 1e-6)) ./ (v_burial(1) .* rho);  % umol/g
Fe_3_init = Fe_inventory_factor .* Fe_3_init_raw;
% Fe3 top boundary is a fixed FeOOH inventory derived from Fe input using a global Fe_inventory_factor.

% Do not convert F_FeOx into a fixed FeOOH surface concentration here.

% Solve Fe(III) with the current Fe3_ODE instead of explicit forward update
nmesh = 1000;
x = linspace(0, Lbottom, nmesh);

opts_Fe3 = bvpset('NMax', 5000, 'RelTol', 1e-4);

guess_fun = @(xq) [max(interp1(z_sed, FeooH, xq, 'linear', 'extrap'), 1e-12); 0];
solinit = bvpinit(linspace(0, Lbottom, nmesh), guess_fun);

sol = bvp4c(@Fe3_ODE, @Fe3_bc, solinit, opts_Fe3);

x = linspace(0, Lbottom, n);
y = deval(sol, x);

if min(y(1,:)) < 0
    fprintf('Fe3 is negative in the current iteration! Minimum：%.2e\n', min(y(1,:)));
end

C_Fe_3_new = max(y(1,:), 1e-12);

% Under-relaxed FeOOH update for the NEXT outer iteration.
% Do not recompute R_FeRed inside the same iteration.
relax_Fe3 = 0.5;
FeooH = (1 - relax_Fe3) .* FeooH_old + relax_Fe3 .* C_Fe_3_new;

% Use the same C_to_Fe / R_FeRed that was used to solve Fe2 in this iteration.
RC_after_Fe = max(RC_after_O2 - C_to_Fe, 0);

% % ------------------------ SO4 budget pre-diagnostic ----------------------
% R_SRR_pot_max = 0.5 .* RC_after_Fe;  % maximum sulfate demand, umol SO4/L/yr
% 
% I_SO4_demand_pot = trapz(z_sed, R_SRR_pot_max) .* 1e-3;  % umol SO4/cm2/yr
% 
% % crude lower-bound/top-down supply estimates
% J_SO4_diff_ref = DSO4 .* SO4init ./ Lbottom .* 1e-3;      % umol SO4/cm2/yr
% I_SO4_irrig_max = trapz(z_sed, Alpha_Bioirrig .* SO4init) .* 1e-3;
% 
% SO4_supply_ref = J_SO4_diff_ref + I_SO4_irrig_max;
% 
% fprintf('\n--- SO4 pre-diagnostic ---\n');
% fprintf('Potential SO4 demand from RC_after_Fe: %.3f umol/cm2/yr\n', I_SO4_demand_pot);
% fprintf('Reference diffusive SO4 supply:        %.3f umol/cm2/yr\n', J_SO4_diff_ref);
% fprintf('Max irrigation SO4 supply:             %.3f umol/cm2/yr\n', I_SO4_irrig_max);
% fprintf('Supply reference total:                %.3f umol/cm2/yr\n', SO4_supply_ref);
% fprintf('Demand / supply_ref:                   %.2f\n', I_SO4_demand_pot ./ max(SO4_supply_ref, 1e-12));
% fprintf('--------------------------\n\n');

% ------------------------ SULFATE ---------------------------------------
% Use positivity-preserving finite-volume depletion-front solver.
[Sulfate, R_SRR, R_AOM_actual, SO4_diag] = Solve_SO4_FV_Front();

% fprintf('SO4 FV front: min=%.3e uM, front=%.2f cm, I_SRR=%.3f, I_AOM=%.3f, I_pot=%.3f umol/cm2/yr\n', ...
%     SO4_diag.min_SO4, SO4_diag.front_depth, SO4_diag.I_SRR, SO4_diag.I_AOM, SO4_diag.I_demand_pot);
fprintf('SO4 I_demand_pot = %.3f\n', SO4_diag.I_demand_pot);
fprintf('SO4 I_SRR        = %.3f\n', SO4_diag.I_SRR);
fprintf('SO4 I_AOM        = %.3f\n', SO4_diag.I_AOM);
fprintf('SO4 I_irrig      = %.3f\n', SO4_diag.I_irrig_source);
fprintf('SO4 J_top_down   = %.3f\n', SO4_diag.J_top_down);
fprintf('SO4 supply/demand = %.3f\n', ...
    (SO4_diag.I_irrig_source + SO4_diag.J_top_down) ./ ...
    max(SO4_diag.I_SRR + SO4_diag.I_AOM, 1e-12));
% % ------------------------ SULFATE ---------------------------------------
% 
% % Solving ODE
% 
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% % solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[SO4init 0]);
% sol = bvp4c(@SO4_ODE,@SO4_bc,solinit);
% 
% x = linspace(0,Lbottom,n);
% 
% y = deval(sol,x);
% 
% if min(y) < 0
%     fprintf('SO4 is negative in the current iteration! Minimum：%.2e\n', min(y));
% end
% y = max(y, 1e-12);
% 
% % C_SO4 = y(1,:);
% % Sulfate = C_SO4;
% % 
% % %R_SRR = RC.*Inhib.* (Sulfate./(Sulfate+k_SO4)).*1E9; %rate of sulfate reduction umol/l/year
% % 
% % R_SRR = 0.5 .* RC_after_Fe .* (Sulfate ./ (Sulfate + k_SO4));   % umol SO4/L/yr
% 
% SO4_raw = real(y(1,:));
% 
% if min(SO4_raw) < -1
%     fprintf('SO4 is substantially negative before clipping! Minimum：%.2e\n', min(SO4_raw));
% elseif min(SO4_raw) < 0
%     fprintf('SO4 has small negative numerical values before clipping. Minimum：%.2e\n', min(SO4_raw));
% end
% 
% Sulfate = max(SO4_raw, 0);
% 
% SO4_cut = 1e-6;   % uM; below this sulfate is treated as depleted
% 
% SO4_pos = max(Sulfate, 0);
% 
% f_SRR = SO4_pos ./ max(SO4_pos + k_SO4, 1e-12);
% f_SRR(SO4_pos <= SO4_cut) = 0;
% 
% R_SRR = 0.5 .* RC_after_Fe .* f_SRR;   % actual organoclastic SRR, umol SO4/L/yr
% 
% % Actual AOM supported by current sulfate, using previous CH4-derived potential.
% f_AOM = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
% f_AOM(SO4_pos <= SO4_cut) = 0;
% 
% R_AOM_actual = R_AOM_pot .* f_AOM;     % actual AOM used by sulfur/carbonate ledger

% ------------------------ SULFIDE ---------------------------------------

% Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@HS_ODE,@HS_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('HS is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_HS = max(real(y(1,:)), 0);
HS_conc = C_HS ./ (1 + ((10.^(6-pH))./K_HS));
Sulfide(iteration,:) = C_HS;   

% sigma_FeS_1 =  (Iron_C(iteration,:).*HS_conc)./((10.^(6-pH)).*KFeS);
% delta_FeS  = (sigma_FeS_1 - 1);
% 
% for i=1:n
%   if delta_FeS (1,i) > 0
%       delta_FeS1(1,i) = 1;
%   else 
%       delta_FeS1(1,i) = 0;
%   end 
% end

% R_FeS = kFeS.*C_Fe.*C_HS;
R_FeS = kFeS .* C_Fe .* HS_conc;

% R_FeS_1(count_loop,:) = R_FeS;
% R_FeS_store(iteration,:) = R_FeS;

for i=1:n
  if R_FeS(1,i) < 0
      R_FeS(1,i) = 0;
  end 
end

R_HS_Ox = max(Kreox .* HS_conc .* Oxygen, 0);   % umol/L/yr

% ---------------- Fe diagnostics ----------------

% I_FeRed = trapz(z_sed, R_FeRed) .* 1e-3;
% I_FeOx  = trapz(z_sed, R_FeOx)  .* 1e-3;
% I_FeS   = trapz(z_sed, R_FeS)   .* 1e-3;
% 
% I_Fe_irrig = trapz(z_sed, Alpha_Bioirrig .* max(C_Fe - Feinit, 0)) .* 1e-3;
% 
% J_Fe2_top_up = DH2S .* ((C_Fe(2) - C_Fe(1)) ./ dz_sed) .* 1e-3;
% 
% Fe_sink_total = I_FeOx + I_FeS + I_Fe_irrig + J_Fe2_top_up;

% fprintf('\n--- Fe diagnostics ---\n');
% fprintf('I_FeRed_pot          = %.3f umol Fe/cm2/yr\n', I_FeRed_pot);
% fprintf('I_Fe external supp = %.3f umol Fe/cm2/yr\n', I_Fe_supply_ext);
% fprintf('I_Fe transport supply= %.3f umol Fe/cm2/yr\n', I_Fe_transport_supply);
% fprintf('I_Fe recycle lag     = %.3f umol Fe/cm2/yr\n', I_Fe_recycle_lag);
% fprintf('I_Fe supply cap      = %.3f umol Fe/cm2/yr\n', I_Fe_supply_cap);
% fprintf('Fe supply scale      = %.3f\n', Fe_supply_scale);
% fprintf('I_FeRed          = %.3f umol Fe/cm2/yr\n', I_FeRed);
% fprintf('I_FeOx           = %.3f umol Fe/cm2/yr\n', I_FeOx);
% fprintf('I_FeS            = %.3f umol Fe/cm2/yr\n', I_FeS);
% fprintf('I_Fe irrigation  = %.3f umol Fe/cm2/yr\n', I_Fe_irrig);
% fprintf('J_Fe2 top up     = %.3f umol Fe/cm2/yr\n', J_Fe2_top_up);
% fprintf('Fe sink total    = %.3f umol Fe/cm2/yr\n', Fe_sink_total);
% fprintf('max Fe2          = %.1f uM\n', max(C_Fe));
% fprintf('max Fe3          = %.1f umol/g\n', max(FeooH));
% fprintf('Fe sink/source   = %.3f\n', Fe_sink_total / max(I_FeRed,1e-12));
% 
% Fe_budget_residual = Fe_sink_total - I_FeRed;
% 
% fprintf('Fe residual       = %.3f umol Fe/cm2/yr\n', Fe_budget_residual);
% fprintf('Fe residual/source= %.3f\n', Fe_budget_residual ./ max(I_FeRed,1e-12));
% Fe_sink_no_top = I_FeOx + I_FeS + I_Fe_irrig;
% 
% 
% fprintf('Fe sink no top    = %.3f umol Fe/cm2/yr\n', Fe_sink_no_top);
% fprintf('Fe no-top/source  = %.3f\n', Fe_sink_no_top ./ max(I_FeRed,1e-12));
% 
% HS_used_for_Fe2 = C_HS_used_for_Fe2 ./ (1 + ((10.^(6 - pH_used_for_Fe2)) ./ K_HS));
% R_FeS_used_for_Fe2 = kFeS .* C_Fe .* HS_used_for_Fe2;
% I_FeS_used_for_Fe2 = trapz(z_sed, R_FeS_used_for_Fe2) .* 1e-3;
% 
% fprintf('I_FeS used Fe2   = %.3f umol Fe/cm2/yr\n', I_FeS_used_for_Fe2);
% fprintf('Fe used sink/source = %.3f\n', ...
%     (I_FeOx + I_FeS_used_for_Fe2 + I_Fe_irrig + J_Fe2_top_up) ./ max(I_FeRed,1e-12));


fprintf('SO4 min = %.2f uM\n', min(Sulfate));
fprintf('SO4 bottom = %.2f uM\n', Sulfate(end));

idx_SO4_10 = find(Sulfate < 10, 1);
if isempty(idx_SO4_10)
    fprintf('SO4 depletion front <10 uM: none\n');
else
    fprintf('SO4 depletion front <10 uM: %.2f cm\n', z_sed(idx_SO4_10));
end

I_FeRed = trapz(z_sed, R_FeRed) .* 1e-3;
I_FeS   = trapz(z_sed, R_FeS) .* 1e-3;
I_FeOx_int = trapz(z_sed, R_FeOx) .* 1e-3;

fprintf('I_FeRed = %.3f umol Fe/cm2/yr\n', I_FeRed);
fprintf('I_FeS   = %.3f umol Fe/cm2/yr\n', I_FeS);
fprintf('I_FeOx  = %.3f umol Fe/cm2/yr\n', I_FeOx_int);
fprintf('FeS/FeRed = %.3f\n', I_FeS ./ max(I_FeRed,1e-12));
fprintf('max Fe2 = %.2f uM\n', max(C_Fe));
fprintf('max FeOOH = %.2f umol/g\n', max(FeooH));

fprintf('----------------------\n\n');

% ------------------------ METHANE ---------------------------------------
RC_after_SO4 = max(RC_after_Fe - 2 .* R_SRR, 0);   % umol C/L/yr

Rate_Meth = 0.5 .* RC_after_SO4;                   % umol CH4/L/yr


% Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@CH4_ODE,@CH4_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('CH4 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_CH4 = y(1,:);
CH4 = C_CH4;
SO4_pos = max(Sulfate, 0);
f_AOM_post = SO4_pos ./ max(SO4_pos + K_CH4_SO4, 1e-12);
f_AOM_post(SO4_pos <= 1e-6) = 0;

R_AOM   = k_AOM .* CH4 .* f_AOM_post;  % diagnostic current AOM based on current CH4/SO4
R_CH4Ox = k_aerobic_CH4 .* CH4 .* (Oxygen ./ max(Oxygen + K_CH4_O2, 1e-12));

% Use the same actual AOM that SO4/HS used this iteration.
% This keeps sulfur and carbonate ledger internally consistent.
R_DIC_prod = R_respi ...
           + R_FeRed ./ 4 ...
           + 2 .* R_SRR ...
           + Rate_Meth ...
           + R_CH4Ox ...
           + R_AOM_actual;

R_ALK_prod = 0.5 .* R_FeRed ...
           + 2 .* R_SRR ...
           + 2 .* R_AOM_actual ...
           - 2 .* R_FeS ...
           - R_HS_Ox;

% Store potential AOM for the next iteration's SO4 solve.
R_AOM_pot = k_AOM .* CH4;
R_AOM_lag = R_AOM;

% ---------------- Redox partition diagnostics ----------------
% I_RC     = trapz(z_sed, RC .* 1E9) .* 1e-3;        % umol C/cm2/yr
% I_respi  = trapz(z_sed, R_respi) .* 1e-3;          % umol C/cm2/yr
% I_FeC    = trapz(z_sed, R_FeRed ./ 4) .* 1e-3;     % umol C/cm2/yr
% I_SO4C   = trapz(z_sed, 2 .* R_SRR) .* 1e-3;       % umol C/cm2/yr
% I_methC  = trapz(z_sed, 2 .* Rate_Meth) .* 1e-3;   % umol C/cm2/yr
% 
% fprintf('\n--- Carbon redox partition ---\n');
% fprintf('I_RC total       = %.3f umol C/cm2/yr\n', I_RC);
% fprintf('I_O2 respiration = %.3f (%.1f%%)\n', I_respi, 100*I_respi/max(I_RC,1e-12));
% fprintf('I_Fe reduction   = %.3f (%.1f%%)\n', I_FeC,   100*I_FeC/max(I_RC,1e-12));
% fprintf('I_SO4 reduction  = %.3f (%.1f%%)\n', I_SO4C,  100*I_SO4C/max(I_RC,1e-12));
% fprintf('I_methanogenesis = %.3f (%.1f%%)\n', I_methC, 100*I_methC/max(I_RC,1e-12));
% fprintf('Closure ratio     = %.3f\n', (I_respi + I_FeC + I_SO4C + I_methC)/max(I_RC,1e-12));
% fprintf('------------------------------\n\n');

% ------------------------------- Coupled Carbonate -----------------------------------
CaCO3_init = 1E-4 .* (F_CaCO3) .* (poros(1) / (1 - poros(1))) / (v_burial(1) * rho);  % gr/grDw

solinit_coupled = bvpinit(linspace(0, Lbottom, nmesh), ...
    [DICinit, 0, HCO3init, 0, CaCO3_init, 0]);

sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled);

y_coupled = deval(sol_coupled, x);

C_DIC  = max(real(y_coupled(1,:)), 1e-12);
C_alka = max(real(y_coupled(3,:)), 1e-12);
CaCO3  = max(real(y_coupled(5,:)), 0);

ALK = C_alka;

F_diff_DIC = DHCO3 .* ((C_DIC(1,2)  - C_DIC(1,1))  ./ (x(1,2) - x(1,1))) * 1E-3;
F_diff     = DHCO3 .* ((C_alka(1,2) - C_alka(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;

% -------------- pH / CO3 / H2CO3 postprocessing -----------------------
pH_1    = zeros(1,n);
CO3_1   = zeros(1,n);
C_H2CO3 = zeros(1,n);

for i = 1:n
    [pH_1(i), CO3_1(i), C_H2CO3(i)] = River_Carbonate(ALK(i), C_DIC(i), T_future, Salinity, 1);
end

pH = pH_1;

sigma_carb = (Calcium_activity .* Calcium .* CO3_activity .* CO3_1) ./ Ksp_ca - 1;

unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho .* max(1 - poros, 1e-6));

R_carb_form = zeros(1,n);
R_carb_disso = zeros(1,n);

for i = 1:n
    if sigma_carb(i) > 0
        R_carb_form(i) = abs(sigma_carb(i))^n_power_CaCO31 .* k_calcite .* unit_conversion(i);
    else
        R_carb_form(i) = 0;
    end

    if (sigma_carb(i) < 0) && (sigma_carb(i) > -0.2)
        R_carb_disso(i) = abs(sigma_carb(i))^n_power_CaCO32 .* k_calcite_dis1 .* CaCO3(i);
    elseif sigma_carb(i) <= -0.2
        R_carb_disso(i) = abs(sigma_carb(i))^n_power_CaCO33 .* k_calcite_dis2 .* CaCO3(i);
    else
        R_carb_disso(i) = 0;
    end
end

R1_carb = R_carb_form - R_carb_disso .* (1E3 * 1E6 * 1E-2 * rho .* max(1 - poros, 1e-6));

% ------------------------ PHOSPHOROUS ------------------------------------

% % Solving ODE
% 
% nmesh=1000;
% x=linspace(0,Lbottom,nmesh);
% solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
% sol = bvp4c(@Phos_ODE,@Phos_bc,solinit);
% 
% x = linspace(0,Lbottom,n);
% 
% y = deval(sol,x);
% C_Phos = y(1,:);
% PO4 = C_Phos;
% 
% P_apaeq = (305.831-74.654.*pH + 4.583.*pH.^2).*(1-(T_apa-10).*X_apa);
% FeOxide=interp1(z_sed,FeooH./1000,x);
% Kadsp = (Kp_ads.*FeOxide.*S_ads_coef.*((1-poros)./poros).*rho.*1E3)./((10.^-pH.*1E6)+Kp_ads.*C_Phos);
% ads_P = Kadsp.*C_Phos; %Fe-bound P
% 
% PO4_diss = C_Phos - ads_P; %PO4 dissolved
% delta_apat2 = (PO4_diss-P_apaeq);
% 
% for i=1:n
%   if delta_apat2 (1,i) > 0
%       delta_apat1(1,i) = 1;
%   else 
%       delta_apat1(1,i) = 0;
%   end 
% end 
% 
% sigma1_viv = ((Iron_C.^3).*(PO4_diss.^2))./Kviv; %umol/Lsed/yr
% delta_viv1 = (sigma1_viv-1);
% 
% for i=1:n
%   if delta_viv1 (1,i) > 0
%       delta_viv1(1,i) = 1;
%   else 
%       delta_viv1(1,i) = 0;
%   end 
% end 
% 
% Rapat = kapatite.*delta_apat1.*(delta_apat2.^n_apa);
% Rviv1 = kviv.*(sigma1_viv.^alpha_viv-1).*delta_viv1.*(1-poros).*rho.*1E3; %umol/Lsed/yr

% -------------------------- Convergence coefficient ----------------------

if iteration > 1
    abs_Fe2 = max(abs(C_Fe(:) - C_Fe_prev(:)));
    rel_Fe2 = abs_Fe2 ./ max(max(abs(C_Fe_prev(:))), 1);
    
    % If Fe2 is low, absolute uM-scale changes should not dominate convergence.
    conv_Fe2 = min(rel_Fe2, abs_Fe2 ./ 25);
    conv_Fe3 = max(abs(FeooH(:) - FeooH_prev(:))) ./ max(max(abs(FeooH_prev(:))), 1);
    conv_SO4 = max(abs(Sulfate(:) - Sulfate_prev(:))) ./ max(max(abs(Sulfate_prev(:))), 1);
    conv_CH4 = max(abs(CH4(:) - CH4_prev(:))) ./ max(max(abs(CH4_prev(:))), 1);
    conv_ALK = max(abs(ALK(:) - ALK_prev(:))) ./ max(max(abs(ALK_prev(:))), 1);
    conv_pH  = max(abs(pH(:) - pH_prev(:))) ./ max(max(abs(pH_prev(:))), 1);
conv_FeS = max(abs(R_FeS(:) - R_FeS_prev(:))) ./ max(max(abs(R_FeS_prev(:))), 1);

abs_HS = max(abs(C_HS(:) - C_HS_prev(:)));
conv_HS_abs = abs_HS ./ 200;   % diagnostic only, not hard convergence

conv_all = max([conv_Fe2, conv_Fe3, conv_SO4, conv_CH4, conv_ALK, conv_pH, conv_FeS]);

fprintf('Outer iteration %d convergence: Fe2 %.3f, Fe3 %.3f, SO4 %.3f, HSabs %.3f, FeS %.3f, CH4 %.3f, ALK %.3f, pH %.3f, max %.3f\n', ...
    iteration, conv_Fe2, conv_Fe3, conv_SO4, conv_HS_abs, conv_FeS, conv_CH4, conv_ALK, conv_pH, conv_all);
else
    conv_all = Inf;
end

C_Fe_prev    = C_Fe;
FeooH_prev   = FeooH;
Sulfate_prev = Sulfate;
CH4_prev     = CH4;
ALK_prev     = ALK;
pH_prev      = pH;
C_HS_prev = C_HS;
R_FeS_prev = R_FeS;

if iteration >= min_outer_iter && conv_all < conv_tol
    break
end

iteration = iteration + 1;
count_loop = count_loop + 1;
end

% 
% % -------------------------- Sulfur budget diagnostics --------------------------
% dz = z_sed(2) - z_sed(1);
% 
% % Top boundary diffusive fluxes
% % Sign convention used here:
% %   J_SO4_top_down > 0 : downward sulfate flux INTO sediment
% %   J_HS_top_up   > 0 : upward sulfide flux OUT OF sediment
% J_SO4_top_down = -DSO4 .* ((Sulfate(2) - Sulfate(1)) ./ dz) .* 1e-3;   % umol/cm2/yr
% J_HS_top_up    =  DH2S .* ((C_HS(2)   - C_HS(1))   ./ dz) .* 1e-3;     % umol/cm2/yr
% 
% % Depth-integrated rates (convert umol/L/yr over cm to umol/cm2/yr)
% I_SRR   = trapz(z_sed, R_SRR)     .* 1e-3;
% I_AOM   = trapz(z_sed, R_AOM_lag) .* 1e-3;
% I_FeS   = trapz(z_sed, R_FeS)     .* 1e-3;
% I_HSOx  = trapz(z_sed, R_HS_Ox)   .* 1e-3;
% 
% % Irrigation terms as net source/sink diagnostics
% I_SO4_irrig = trapz(z_sed, Alpha_Bioirrig .* (SO4init - Sulfate)) .* 1e-3;   % >0 means irrigation supplies SO4
% I_HS_irrig  = trapz(z_sed, Alpha_Bioirrig .* (C_HS - HSinit))    .* 1e-3;    % >0 means irrigation removes HS
% 
% % Net sulfur source/sink summaries
% SO4_total_sink   = I_SRR + I_AOM;
% SO4_total_supply = J_SO4_top_down + I_SO4_irrig;
% 
% HS_total_source  = I_SRR + I_AOM;
% HS_total_sink    = I_FeS + I_HSOx + I_HS_irrig + J_HS_top_up;
% 
% % Where is SRR happening?
% idx_deep10 = z_sed >= 10;
% idx_deep15 = z_sed >= 15;
% 
% I_SRR_deep10 = trapz(z_sed(idx_deep10), R_SRR(idx_deep10)) .* 1e-3;
% I_SRR_deep15 = trapz(z_sed(idx_deep15), R_SRR(idx_deep15)) .* 1e-3;
% 
% frac_SRR_deep10 = I_SRR_deep10 ./ max(I_SRR, 1e-12);
% frac_SRR_deep15 = I_SRR_deep15 ./ max(I_SRR, 1e-12);
% 
% % Sulfide sink partition
% frac_FeS   = I_FeS   ./ max(HS_total_source, 1e-12);
% frac_HSOx  = I_HSOx  ./ max(HS_total_source, 1e-12);
% frac_irrig = I_HS_irrig ./ max(HS_total_source, 1e-12);
% frac_flux  = J_HS_top_up ./ max(HS_total_source, 1e-12);
% 
% fprintf('\n================ Sulfur Budget Diagnostics ================\n');
% fprintf('Top SO4 diffusive influx (downward +):   %.3f umol/cm2/yr\n', J_SO4_top_down);
% fprintf('Integrated SO4 irrigation source:        %.3f umol/cm2/yr\n', I_SO4_irrig);
% fprintf('Integrated organoclastic SRR sink:       %.3f umol/cm2/yr\n', I_SRR);
% fprintf('Integrated AOM sulfate sink:             %.3f umol/cm2/yr\n', I_AOM);
% fprintf('SO4 total supply (top+irrigation):       %.3f umol/cm2/yr\n', SO4_total_supply);
% fprintf('SO4 total sink   (SRR+AOM):              %.3f umol/cm2/yr\n', SO4_total_sink);
% 
% fprintf('\nIntegrated HS source (SRR+AOM):          %.3f umol/cm2/yr\n', HS_total_source);
% fprintf('Integrated FeS sink:                     %.3f umol/cm2/yr (%.2f%%)\n', I_FeS,   100*frac_FeS);
% fprintf('Integrated HS oxidation sink:            %.3f umol/cm2/yr (%.2f%%)\n', I_HSOx,  100*frac_HSOx);
% fprintf('Integrated HS irrigation sink:           %.3f umol/cm2/yr (%.2f%%)\n', I_HS_irrig, 100*frac_irrig);
% fprintf('Top HS diffusive efflux (upward +):      %.3f umol/cm2/yr (%.2f%%)\n', J_HS_top_up, 100*frac_flux);
% fprintf('HS total sink (FeS+Ox+irrig+flux):       %.3f umol/cm2/yr\n', HS_total_sink);
% 
% fprintf('\nDeep SRR fraction (>10 cm):              %.2f%%\n', 100*frac_SRR_deep10);
% fprintf('Deep SRR fraction (>15 cm):              %.2f%%\n', 100*frac_SRR_deep15);
% fprintf('===========================================================\n\n');

% Storing steady-state solutions as initial conditions for PDEs

Diff_fluxes = [F_diff_DIC F_diff F_diff_DIC./F_diff]';
R_ALK_DIC = F_diff./F_diff_DIC;  
F_diff_CH4 = DCH4.*((CH4(1,2) - CH4(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr


% % ----------------------------- OUTPUT PACKAGING ---------------------------
%     Outputs.z_sed = z_sed;
%     Outputs.pH_profile = pH;
%     Outputs.CH4_profile = CH4;
%     Outputs.O2_profile = Oxygen;
% %     Outputs.SO4_profile = Sulfate;
% %     Outputs.DIC_profile = C_DIC;
% 
%     % Core Diagnostics
%     Outputs.Max_CH4 = max(CH4);
% 
% 
%     Outputs.Org_Bottom = C_organic(end) * 100; % %gDw
%     Outputs.ALK_Bottom = ALK(end);             % uM
%     Outputs.pH_Bottom  = pH(end);              % 
%     Outputs.CH4_Bottom = CH4(end);             % uM
% %     Outputs.Org_Top    = C_organic(1) * 100; % %gDw
% 
%     % OPD: O2 < 1 uM
%     idx_O2 = find(Oxygen < 1, 1);
%     if isempty(idx_O2), Outputs.OPD = z_sed(end); else, Outputs.OPD = z_sed(idx_O2); end
% 
%     % SO4_Depth: SO4 降至 < 10 uM
%     idx_SO4 = find(Sulfate < 10, 1);
%     if isempty(idx_SO4), Outputs.SO4_Depth = z_sed(end); else, Outputs.SO4_Depth = z_sed(idx_SO4); end
% 
% 
%     idx_top5 = (z_sed <= 5);                   % 圈定 0-5 cm 网格
%     idx_bot5 = (z_sed >= (Lbottom - 5));       % 圈定底部 5 cm 网格
%     
%     Outputs.ALK_Bot5   = mean(ALK(idx_bot5));             % 底层 5cm 平均碱度
%     Outputs.Sigma_Top5 = mean(sigma_carb(idx_top5));      % 表层 5cm 平均饱和度 (Omega-1)
%     Outputs.CaCO3_Top5 = mean(CaCO3(idx_top5)) * 100;     % 表层 5cm 平均 CaCO3 (%gDw)
%     Outputs.Integ_Meth = trapz(z_sed, Rate_Meth);         %integrated Rate_Meth
%     
% %     % CH4_Onset_Depth: CH4 超过 10 uM 的深度
% %     idx_CH4 = find(CH4 > 10, 1);
% %     if isempty(idx_CH4), Outputs.CH4_Onset = z_sed(end); else, Outputs.CH4_Onset = z_sed(idx_CH4); end
% %     
% %     % Sigma0_Depth: 碳酸钙饱和度 Omega-1 穿过 0 的深度 (>= 0)
% %     idx_sigma = find(sigma_carb >= 0, 1);
% %     if isempty(idx_sigma), Outputs.Sigma0_Depth = z_sed(end); else, Outputs.Sigma0_Depth = z_sed(idx_sigma); end
% %     
% %     % CaCO3_Front_Depth: 碳酸钙开始显著积累的深度 (设定阈值为 1e-4，即脱离初始极小值)
% %     idx_CaCO3 = find(CaCO3 > 1e-4, 1);
% %     if isempty(idx_CaCO3), Outputs.CaCO3_Front = z_sed(end); else, Outputs.CaCO3_Front = z_sed(idx_CaCO3); end
% 
% 
% %     % Methane Appearance Depth (Depth where CH4 > 5 uM)
% %     ch4_idx = find(CH4 > 5, 1);
% %     if isempty(ch4_idx)
% %         Outputs.CH4_Depth = Lbottom; % No significant methane
% %     else
% %         Outputs.CH4_Depth = z_sed(ch4_idx);
% %     end
% % 
% %     Outputs.Convergence_Status = K_converge;
% 

% end % End of Function


% -------------------------------------------------------------------------
% ----------------------------- PLOTS -------------------------------------
clf;
n_plot = 6; % number of plots in each row
m_plot = 3; % number of total rows


% Organic

subplot(m_plot,n_plot,1);

% plot((C_organic + POC_root).*100,z_sed,'lineWidth',2); axis ij
plot((C_organic ).*100,z_sed,'lineWidth',2); axis ij
title('Organic (%gDw)')
ylabel('Depth (cm)');
box on


% Oxygen

subplot(m_plot,n_plot,2);

plot(Oxygen,z_sed,'lineWidth',2); axis ij
title('[O_2] (\muM)')
box on
grid on

ax.LineWidth = 2;

% Iron

subplot(m_plot,n_plot,3);

plot(C_Fe,z_sed,'lineWidth',2); axis ij
title('[Fe^{2+}] (\muM)')
box on
grid on

ax.LineWidth = 2;

% Sulfate

subplot(m_plot,n_plot,4);

plot(Sulfate,z_sed,'lineWidth',2); 
% xlim([0 SO4init]); 
axis ij
title('[SO_4] (\muM)')
box on
grid on

ax.LineWidth = 2;

% Sulfide

subplot(m_plot,n_plot,5);

plot(C_HS,z_sed,'lineWidth',2); axis ij
title('[H_2S] (\muM)')
box on
grid on

ax.LineWidth = 2;

% Methane

subplot(m_plot,n_plot,6);

plot(CH4,z_sed,'lineWidth',2); axis ij
title('[CH_4] (\muM)')
box on
grid on

ax.LineWidth = 2;

% CaCO3

subplot(m_plot,n_plot,7);

% plot(CaCO3.*(1E-5.*(poros2./(1-poros2)).*(1./rho)),z_sed,'lineWidth',2); axis ij
plot(CaCO3.*100,z_sed,'lineWidth',2); axis ij
title('CaCO3')
ylabel('Depth (cm)');
box on
grid on

ax.LineWidth = 2;

% DIC

subplot(m_plot,n_plot,8);

plot(C_DIC,z_sed,'lineWidth',2); axis ij
title('DIC (\muM)')
box on
grid on

ax.LineWidth = 2;

% ALK

subplot(m_plot,n_plot,9);

plot(C_alka,z_sed,'lineWidth',2); axis ij
title('ALK (\muM)')
box on
grid on

ax.LineWidth = 2;

% Carbonic Acid

subplot(m_plot,n_plot,10);

plot(C_H2CO3,z_sed,'lineWidth',2); axis ij
title('Carb Acid (\muM)')
box on
grid on

ax.LineWidth = 2;

% pH
subplot(m_plot,n_plot,11);

plot(pH,z_sed,'lineWidth',2); axis ij
title('pH')
box on
grid on

ax.LineWidth = 2;

% Organic degradation rate

subplot(m_plot,n_plot,12);

plot(RC.* 1E9,z_sed,'lineWidth',2); axis ij  %umol/l/year
title('Mineralization Rate (\mumol/l/year)')
box on
grid on

ax.LineWidth = 2;

% Burial efficiency of Organic

subplot(m_plot,n_plot,13);
plot(BEsed_org.*100,z_sed,'lineWidth',2); axis ij  %umol/l/year
title('OM Burial Efficiency')
ylabel('Depth (cm)');
box on
grid on

ax.LineWidth = 2;

% Aerobic respiration and sulfate reduction rates

% subplot(m_plot,n_plot,14);
% plot(R_SRR,z_sed,R_respi,z_sed,'lineWidth',2); axis ij %umol/l/year
% title('Rate (\mumol/l/year)')
% legend('Sulfate Red','Aerobic Resp');


subplot(m_plot,n_plot,15);
plot((0.5.*R_SRR)./365,z_sed,'lineWidth',2); axis ij %umol/l/year
title('Sulfate Reduction Rate (nmol/cm3/d)')
% legend('Sulfate Red');

box on
grid on

ax.LineWidth = 2;

% Carbonate saturation index

% subplot(m_plot,n_plot,15);
% 
% plot(sigma_carb(iteration-1,:),z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Caclite saturation (\Omega - 1)')
% 
% box on


subplot(m_plot,n_plot,17);

plot(sigma_carb(end,:),z_sed,'lineWidth',2); axis ij
title('Calciite saturation (\Omega - 1)')

box on
grid on

ax.LineWidth = 2;

subplot(m_plot,n_plot,18);

plot(FeooH(1,:),z_sed,'lineWidth',2); axis ij
title('Fe(III) (\mumol/gr)')

box on
grid on

ax.LineWidth = 2;

toc
% Mineral saturation indices
% 
% subplot(m_plot,n_plot,16);
% plot(delta_viv1,z_sed,delta_apat2,z_sed,delta_FeS,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Mineral saturation (\Omega - 1)')
% legend('Vivianite','Apatite','FeS');
% 
% box on
% grid on
% 
% ax.LineWidth = 2;


% AAA_time = toc;
% AAA_data = [Oxygen' C_DIC' C_alka' pH' sigma_carb' z_sed'];
% 
% R_iron_unit = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9.*...
%                        (poros./(1-poros)).*1E-3.*(1./rho); %rate of iron reduction umol/g/year
% R_FeOx_unit = (kFeOx.*C_Fe.*Oxygen).*...
%               (poros./(1-poros)).*1E-3.*(1./rho);
% 
% R_SRR_integ = cumsum(0.5.*R_SRR.*dz_sed.*1E-3);
% R_RC_integ = cumsum(RC.*dz_sed);
% R_iron_integ = cumsum(R_iron(count_loop-1,:).*dz_sed.*1E-3);
% R_iron_integ_unit = cumsum(R_iron_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
% R_FeOx_integ_unit = cumsum(R_FeOx_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
% R_FeOx_integ = cumsum(R_FeOx_1(count_loop-1,:).*dz_sed.*1E-3);
% R_FeS_integ = cumsum(R_FeS.*dz_sed.*1E-3);
% R_FeS_1_integ = cumsum(R_FeS_1(count_loop-1,:).*dz_sed.*1E-3);
% R_HSOX_integ = cumsum(R_Ox.*dz_sed.*1E-3);
% R_biorrig_integ = cumsum((Alpha_Bioirrig.*(HSinit-C_HS)).*dz_sed.*1E-3);
% R_biorrig_integ_iron = cumsum((Alpha_Bioirrig.*(Feinit-C_Fe)).*dz_sed.*1E-3);
% R_biorrigALK_integ = cumsum((Alpha_Bioirrig.*(HCO3init-ALK)).*dz_sed.*1E-3);
% R_ALK = RC.*1E9 - R_respi;
% R_carb_integ = cumsum(R1_carb.*dz_sed.*1E-3);
% R_ALK_integ_WITHOUT = 2.*R_FeS_integ;
% R_ALK_integ_WITH = 2.*R_FeS_integ - (R_carb_integ);
% R_CH4O2_integ = cumsum((k_aerobic_CH4.* CH4.* (Oxygen./(Oxygen+K_CH4_O2))).*dz_sed.*1E-3);
% R_CH4SO4_integ = cumsum((k_AOM.* CH4.* (Sulfate./(Sulfate+K_CH4_SO4))).*dz_sed.*1E-3);
% 
% R_net = R_SRR_integ - R_FeS_integ - R_HSOX_integ + R_biorrig_integ - F_diff_HS;
% 
% R_net_iron = R_iron_integ - R_FeS_integ - R_FeOx_integ - F_diff_Fe(1,count_loop-1) + R_biorrig_integ_iron;
% 
% R_net_Fe3  =  (v_burial(end).*C_Fe_3(end).*rho.*((1-poros(end))./poros(end))) - ...
%               (v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1))) + R_iron_integ_unit - R_FeOx_integ_unit;
% 
% R_net_Fe3_percent  =  100.*R_net_Fe3(end)./(v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1)));
% 
% R_net_org  =  (v_burial(end).*C_organic(end).*rho.*((1-poros(end))./12)) - (v_burial(1).*C_organic(1).*rho.*((1-poros(1))./12)) + ...
%               R_RC_integ;
% 
% R_net_org_percent  =  100.*R_net_org(end)./((rho.*((1-poros(1))./12)).*v_burial(1).*C_organic(1));
% 
% F_HS_tot = (R_SRR_integ - R_FeS_integ).*0.0274;
% 
% F_S_out = (R_SRR_integ - R_FeS_integ - R_HSOX_integ).*0.0274;
% 
% R_SRR_integ_store = R_SRR_integ(end).*0.0274;  %mmol/m2/d
% R_ALK_integ_WITH_store = R_ALK_integ_WITH(end).*0.0274; %mmol/m2/d
% 
% F_ox_py = R_FeS_integ./R_SRR_integ;
% 
% AAA_Store_1 = [F_FeOx R_SRR_integ_store R_ALK_integ_WITH_store];
% AAA_Store = [NPP.*BE R_ALK_integ_WITH(end) R_ALK_integ_WITHOUT(end) 2.*R_carb_integ(end) F_diff];