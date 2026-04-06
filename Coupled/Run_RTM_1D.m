%% Sediment Diagenesis Model (ODEs for organic matter, oxygen, sulfate, DIC, ALK, Carb Acid, and pH)

% -------- Process included in the model ----------------------------------
% Organic matter: Organic matter degradation using Middelburg 
% Oxygen: Oxic respiration 
% Sulfate: Anoxic respiration through sulfate reduction.
% Dissolved inorganic carbon: Organic matter degradation and carbonate precipitation
% Alkalinity: Sulfate reduction and carbonate precipitation

function Outputs = Run_RTM_1D(Custom_Config, Custom_Params)
% RUN_RTM_1D
% Sediment Diagenesis Model - Functionized for Sensitivity Analysis

% ----------------------------- INPUT PARAMETERS ---------------------------

global v_burial Mineral_Mass z_sed Oxygen Sulfate
global k_sed k_O2 DSO4 DH2S DO2 DPO4 k_SO4 Kreox  Bioturb Calcium DHCO3 HCO3init 
global O2init SO4init HSinit C_organic rho poros RC Alpha_Bioirrig
global R_respi R_SRR Ksp_ca k_calcite DICinit R1_carb CO3_1 BE P_C_ratio Rviv1 R_FeS  R_FeOx Fe_3_init
global v_burial_Fluid CO3_activity Calcium_activity NPP kFeS FeooH Feinit R_HS_Ox kapatite 
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2 CH4init Pinitial DCH4 kFeOx KFEMonod Sulfide Rapat CaCO3 F_CaCO3 O2_root
global C_HS C_Fe n_power_CaCO31 n_power_CaCO32 k_calcite_dis1 n_power_CaCO33 k_calcite_dis2 CaCO3_init Temp_factor T_future
%global KFe_HS Iron_conc R_iron Iron_C P_apaeq R1_carb_disso R1_carb_form
    

    if nargin < 1 || isempty(Custom_Config)
        Config = Config_Baseline();
    else
        Config = Custom_Config;
    end
    if nargin < 2 || isempty(Custom_Params)
        Params = Params_Static();
    else
        Params = Custom_Params;
    end
%         Params = Params_Static();
%     Config = Config_Baseline();
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
    KFeS = Params.KFeS;
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
    
    BE              = Config.BE;
    NPP             = Config.NPP;
    F_CaCO3         = Config.F_CaCO3;
    
    T_future        = Config.T_future;



    BE = Config.BE;
    NPP = Config.NPP;
    if Config.use_hydro_npp_multiplier
        NPP = NPP * Hydro.NPP_multiplier;
    end

    F_FeOx = Config.F_FeOx;
    F_CaCO3 = Config.F_CaCO3;

    T_future = Config.T_future;
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

[~, CO3_top, ~] = River_Carbonate(HCO3init, DICinit, T_future, 0.1, 1); 
% CO3_top = Carb_CO3(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
CO3_1 = CO3_top*ones(1,n);
C_HS  = zeros(1,n); %intial value for sulfide
R_HS_Ox = zeros(1,n);
CaCO3 = 1E5.*zeros(1,n);
Rapat = zeros(1,n);
Rviv1 = zeros(1,n); %umol/Lsed/yr
R_FeS = ones(1,n); %umol/Lsed/yr
pH    = 7.*ones(1,n);

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
BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic

end

% ---------------------------- OXYGEN -------------------------------------

RC = Temp_factor.*k_sed.*C_organic.*rho.*((1-poros)./(12)); % molCorg/cm3/yr mineralization rate
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

    FeOxInit  = (F_FeOx.*36.5)./(v_burial(1) * rho * (1-poros(1)));
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
k_sed(1,1:num_OPD) = 10.^(-b_oxic*log10(age(1,1:num_OPD)) - a_oxic);  % Oxic
k_sed(1,(num_OPD+1):n) = 10.^(-b_anoxic*log10(age(1,(num_OPD+1):n)) - a_anoxic);  % Anoxic


% ------------------- Blue Carbon ------------------------------------

% DOC release 

depth_rootzone = 10; % seagrass root length (cm)
z_root = 0:0.1:depth_rootzone;
mu_root = 4; % value for the center of rootzone in normal distribution
sigma_root = 6; % sigma for normal distribution of flux in the rootzone
DOC_root_2 = DOC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/year
DOC_root_22 = interp1(z_root,DOC_root_2,z_sed);
poros_root = interp1(z_sed,poros,z_sed);

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

K_converge = 1;
iteration = 1;
iteration_tolerance = 0.3;
count_loop = 1;

while abs(K_converge) > iteration_tolerance %for count_loop = 1:5 
    

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

BEsed_org = C_organic./C_organic(1);  % Burial Efficiency of Organic

end

% ---------------------------- OXYGEN -------------------------------------

RC = Temp_factor.*k_sed.*C_organic.*rho.*((1-poros)./(12)) + RC_root; % molCorg/cm3/yr mineralization rate

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

R_respi = RC.* (Oxygen./(Oxygen+k_O2)).*1E9; %rate of aerobic respiration umol/l/year

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

Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
% R_iron(count_loop,:) = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9; %rate of iron reduction umol/l/year
R_FeOx = (kFeOx.*C_Fe.*Oxygen);
% R_FeOx_1(count_loop,:) = R_FeOx;
% Fe_3_init  = 365.*1E2.*(F_FeOx)./(v_burial(1));  %umol/l
Fe_3_init = 36.5.*(F_FeOx).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %umol/l

% KFEMonod = 2000;

% % Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[Fe_3_init 0]);
sol = bvp4c(@Fe3_ODE,@Fe3_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('Fe3：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_Fe_3 = y(1,:);
FeooH = C_Fe_3;
% FeooH = max(C_Fe_3, 0.2 .* FeOx);

% ------------------------ SULFATE ---------------------------------------

% Solving ODE

nmesh=1000;
x=linspace(0,Lbottom,nmesh);
solinit = bvpinit(linspace(0,Lbottom,nmesh),[0 0]);
sol = bvp4c(@SO4_ODE,@SO4_bc,solinit);

x = linspace(0,Lbottom,n);

y = deval(sol,x);

if min(y) < 0
    fprintf('SO4 is negative in the current iteration! Minimum：%.2e\n', min(y));
end
y = max(y, 1e-12);

C_SO4 = y(1,:);
Sulfate = C_SO4;

R_SRR = RC.*Inhib.* (Sulfate./(Sulfate+k_SO4)).*1E9; %rate of sulfate reduction umol/l/year

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

C_HS = y(1,:);
Sulfide(iteration,:) = C_HS;

for i=1:n
  if Sulfide(iteration,i) < 0
      Sulfide(iteration,i) = 0;
  end 
end

HS_conc = Sulfide(iteration,:)./(1+((10.^(6-pH))./K_HS));
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

% R_HS_Ox = Sulfide(iteration,:).* Oxygen.*Kreox; %rate of sulfide reduction umol/l/year
% R_HS_Ox = Kreox .* HS_conc .* Oxygen;%XL

% F_diff_HS = DH2S.*((C_HS(1,2) - C_HS(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr

% ------------------------ METHANE ---------------------------------------

%XL calculating inhibition

Inhib_O2_meth  = (k_O2 ./ (Oxygen + k_O2));
Inhib_Fe_meth  = (KFEMonod ./ (FeooH + KFEMonod));
Inhib_SO4_meth = (k_SO4 ./ (Sulfate + k_SO4));

global Rate_Meth
 

Rate_Meth= 0.5 .* RC .* Inhib_O2_meth .* Inhib_Fe_meth .* Inhib_SO4_meth .* 1E9;

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

% ------------------------------- Coupled Carbonate -----------------------------------

CaCO3_init = 1E-4.*(F_CaCO3).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %gr/grDw

solinit_coupled = bvpinit(linspace(0, Lbottom, nmesh), [DICinit, 0, HCO3init, 0, CaCO3_init, 0]);


sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled);


% options_coupled = bvpset('NMax', 5000, 'RelTol', 1e-2); 
% sol_coupled = bvp4c(@Coupled_Carbonate_ODE, @Coupled_Carbonate_bc, solinit_coupled, options_coupled);

y_coupled = deval(sol_coupled, x);

C_DIC  = max(real(y_coupled(1,:)), 1e-12);
C_alka = max(real(y_coupled(3,:)), 1e-12);
CaCO3  = max(real(y_coupled(5,:)), 0);
ALK    = C_alka;

F_diff_DIC = DHCO3 .* ((C_DIC(1,2) - C_DIC(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;
F_diff     = DHCO3 .* ((C_alka(1,2) - C_alka(1,1)) ./ (x(1,2) - x(1,1))) * 1E-3;



  
% -------------- pH and H2CO3 (2 for 6 calculation) -----------------------

for i=1:n
    

           [pH_1(1,i), CO3_1(1,i), C_H2CO3(1,i)] = River_Carbonate(ALK(1,i), C_DIC(1,i), T_future, 0.1, 1);

end

pH = pH_1;
sigma_carb = (Calcium_activity .* Calcium .* CO3_activity .* CO3_1) ./ Ksp_ca - 1;
unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho .* (1 - poros));
R_carb_form = (sigma_carb > 0) .* abs(sigma_carb).^n_power_CaCO31 .* k_calcite .* unit_conversion;
R_carb_disso = (-0.2 < sigma_carb & sigma_carb < 0) .* abs(sigma_carb).^n_power_CaCO32 .* k_calcite_dis1 .* CaCO3 ...
             + (sigma_carb <= -0.2) .* abs(sigma_carb).^n_power_CaCO33 .* k_calcite_dis2 .* CaCO3;
             
R1_carb = R_carb_form - R_carb_disso .* (1E3 * 1E6 * 1E-2 * rho .* (1 - poros));
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

alpha_converge(1,iteration) = F_diff;

% alpha_converge(isnan(alpha_converge))=[];

     if size(alpha_converge,2) > 2
        
       K_converge = (alpha_converge(1,iteration)-alpha_converge(1,iteration-1))/alpha_converge(1,iteration-1);
  
     end

       iteration = iteration + 1;

count_loop = count_loop + 1;

end  % iteration ends here 


% Storing steady-state solutions as initial conditions for PDEs

Diff_fluxes = [F_diff_DIC F_diff F_diff_DIC./F_diff]';
R_ALK_DIC = F_diff./F_diff_DIC;  
F_diff_CH4 = DCH4.*((CH4(1,2) - CH4(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr


% ----------------------------- OUTPUT PACKAGING ---------------------------
    Outputs.z_sed = z_sed;
    Outputs.pH_profile = pH;
    Outputs.CH4_profile = CH4;
    Outputs.O2_profile = Oxygen;
%     Outputs.SO4_profile = Sulfate;
%     Outputs.DIC_profile = C_DIC;

    % Core Diagnostics
    Outputs.Max_CH4 = max(CH4);


    Outputs.Org_Bottom = C_organic(end) * 100; % %gDw
    Outputs.ALK_Bottom = ALK(end);             % uM
    Outputs.pH_Bottom  = pH(end);              % 
    Outputs.CH4_Bottom = CH4(end);             % uM
%     Outputs.Org_Top    = C_organic(1) * 100; % %gDw

    % OPD: O2 < 1 uM
    idx_O2 = find(Oxygen < 1, 1);
    if isempty(idx_O2), Outputs.OPD = z_sed(end); else, Outputs.OPD = z_sed(idx_O2); end

    % SO4_Depth: SO4 降至 < 10 uM
    idx_SO4 = find(Sulfate < 10, 1);
    if isempty(idx_SO4), Outputs.SO4_Depth = z_sed(end); else, Outputs.SO4_Depth = z_sed(idx_SO4); end


    idx_top5 = (z_sed <= 5);                   % 圈定 0-5 cm 网格
    idx_bot5 = (z_sed >= (Lbottom - 5));       % 圈定底部 5 cm 网格
    
    Outputs.ALK_Bot5   = mean(ALK(idx_bot5));             % 底层 5cm 平均碱度
    Outputs.Sigma_Top5 = mean(sigma_carb(idx_top5));      % 表层 5cm 平均饱和度 (Omega-1)
    Outputs.CaCO3_Top5 = mean(CaCO3(idx_top5)) * 100;     % 表层 5cm 平均 CaCO3 (%gDw)
    Outputs.Integ_Meth = trapz(z_sed, Rate_Meth);         %integrated Rate_Meth
    
%     % CH4_Onset_Depth: CH4 超过 10 uM 的深度
%     idx_CH4 = find(CH4 > 10, 1);
%     if isempty(idx_CH4), Outputs.CH4_Onset = z_sed(end); else, Outputs.CH4_Onset = z_sed(idx_CH4); end
%     
%     % Sigma0_Depth: 碳酸钙饱和度 Omega-1 穿过 0 的深度 (>= 0)
%     idx_sigma = find(sigma_carb >= 0, 1);
%     if isempty(idx_sigma), Outputs.Sigma0_Depth = z_sed(end); else, Outputs.Sigma0_Depth = z_sed(idx_sigma); end
%     
%     % CaCO3_Front_Depth: 碳酸钙开始显著积累的深度 (设定阈值为 1e-4，即脱离初始极小值)
%     idx_CaCO3 = find(CaCO3 > 1e-4, 1);
%     if isempty(idx_CaCO3), Outputs.CaCO3_Front = z_sed(end); else, Outputs.CaCO3_Front = z_sed(idx_CaCO3); end


%     % Methane Appearance Depth (Depth where CH4 > 5 uM)
%     ch4_idx = find(CH4 > 5, 1);
%     if isempty(ch4_idx)
%         Outputs.CH4_Depth = Lbottom; % No significant methane
%     else
%         Outputs.CH4_Depth = z_sed(ch4_idx);
%     end
% 
%     Outputs.Convergence_Status = K_converge;


end % End of Function


% % -------------------------------------------------------------------------
% % ----------------------------- PLOTS -------------------------------------
% clf;
% n_plot = 6; % number of plots in each row
% m_plot = 3; % number of total rows
% 
% 
% % Organic
% 
% subplot(m_plot,n_plot,1);
% 
% % plot((C_organic + POC_root).*100,z_sed,'lineWidth',2); axis ij
% plot((C_organic ).*100,z_sed,'lineWidth',2); axis ij
% title('Organic (%gDw)')
% ylabel('Depth (cm)');
% box on
% 
% 
% % Oxygen
% 
% subplot(m_plot,n_plot,2);
% 
% plot(Oxygen,z_sed,'lineWidth',2); axis ij
% title('[O_2] (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Iron
% 
% subplot(m_plot,n_plot,3);
% 
% plot(C_Fe,z_sed,'lineWidth',2); axis ij
% title('[Fe^{2+}] (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Sulfate
% 
% subplot(m_plot,n_plot,4);
% 
% plot(Sulfate,z_sed,'lineWidth',2); axis ij
% title('[SO_4] (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Sulfide
% 
% subplot(m_plot,n_plot,5);
% 
% plot(C_HS,z_sed,'lineWidth',2); axis ij
% title('[H_2S] (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Methane
% 
% subplot(m_plot,n_plot,6);
% 
% plot(CH4,z_sed,'lineWidth',2); axis ij
% title('[CH_4] (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % CaCO3
% 
% subplot(m_plot,n_plot,7);
% 
% % plot(CaCO3.*(1E-5.*(poros2./(1-poros2)).*(1./rho)),z_sed,'lineWidth',2); axis ij
% plot(CaCO3.*100,z_sed,'lineWidth',2); axis ij
% title('CaCO3')
% ylabel('Depth (cm)');
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % DIC
% 
% subplot(m_plot,n_plot,8);
% 
% plot(C_DIC,z_sed,'lineWidth',2); axis ij
% title('DIC (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % ALK
% 
% subplot(m_plot,n_plot,9);
% 
% plot(C_alka,z_sed,'lineWidth',2); axis ij
% title('ALK (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Carbonic Acid
% 
% subplot(m_plot,n_plot,10);
% 
% plot(C_H2CO3,z_sed,'lineWidth',2); axis ij
% title('Carb Acid (\muM)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % pH
% subplot(m_plot,n_plot,11);
% 
% plot(pH,z_sed,'lineWidth',2); axis ij
% title('pH')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Organic degradation rate
% 
% subplot(m_plot,n_plot,12);
% 
% plot(RC.* 1E9,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Mineralization Rate (\mumol/l/year)')
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Burial efficiency of Organic
% 
% subplot(m_plot,n_plot,13);
% plot(BEsed_org.*100,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('OM Burial Efficiency')
% ylabel('Depth (cm)');
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Aerobic respiration and sulfate reduction rates
% 
% % subplot(m_plot,n_plot,14);
% % plot(R_SRR,z_sed,R_respi,z_sed,'lineWidth',2); axis ij %umol/l/year
% % title('Rate (\mumol/l/year)')
% % legend('Sulfate Red','Aerobic Resp');
% 
% 
% subplot(m_plot,n_plot,15);
% plot((0.5.*R_SRR)./365,z_sed,'lineWidth',2); axis ij %umol/l/year
% title('Sulfate Reduction Rate (nmol/cm3/d)')
% % legend('Sulfate Red');
% 
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% % Carbonate saturation index
% 
% % subplot(m_plot,n_plot,15);
% % 
% % plot(sigma_carb(iteration-1,:),z_sed,'lineWidth',2); axis ij  %umol/l/year
% % title('Caclite saturation (\Omega - 1)')
% % 
% % box on
% 
% 
% subplot(m_plot,n_plot,17);
% 
% plot(sigma_carb(end,:),z_sed,'lineWidth',2); axis ij
% title('Calciite saturation (\Omega - 1)')
% 
% box on
% grid on
% 
% ax.LineWidth = 2;
% 
% subplot(m_plot,n_plot,18);
% 
% plot(FeooH(1,:),z_sed,'lineWidth',2); axis ij
% title('Fe(III) (\mumol/gr)')
% 
% box on
% grid on
% 
% ax.LineWidth = 2;

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