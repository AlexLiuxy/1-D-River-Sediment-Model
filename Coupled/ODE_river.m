%% Sediment Diagenesis Model (ODEs for organic matter, oxygen, sulfate, DIC, ALK, Carb Acid, and pH)

% -------- Process included in the model ----------------------------------
% Organic matter: Organic matter degradation using Middelburg 
% Oxygen: Oxic respiration 
% Sulfate: Anoxic respiration through sulfate reduction.
% Dissolved inorganic carbon: Organic matter degradation and carbonate precipitation
% Alkalinity: Sulfate reduction and carbonate precipitation

clear all
tic

% ----------------------------- INPUT PARAMETERS ---------------------------

global v_burial Mineral_Mass z_sed KFe_HS Oxygen Sulfate
global k_sed k_O2 DSO4 DH2S DO2 DPO4 k_SO4 Kreox Iron_conc Bioturb Calcium DHCO3 HCO3init 
global O2init SO4init HSinit C_organic rho poros RC Alpha_Bioirrig
global R_respi R_SRR Ksp_ca k_calcite DICinit R1_carb CO3_1 BE P_C_ratio Rviv1 R_FeS R_iron R_FeOx Fe_3_init
global v_burial_Fluid CO3_activity Calcium_activity NPP kFeS FeooH Feinit Iron_C R_HS_Ox kapatite P_apaeq R1_carb_disso R1_carb_form
global k_AOM k_aerobic_CH4 K_CH4_SO4 K_CH4_O2 CH4init Pinitial DCH4 kFeOx KFEMonod Sulfide Rapat CaCO3 F_CaCO3 O2_root
global C_HS C_Fe n_power_CaCO31 n_power_CaCO32 k_calcite_dis1 n_power_CaCO33 k_calcite_dis2 CaCO3_init Temp_factor T_future

Bioirrig_top = 100;%80; %1/yr
Bioirrig_bottom = 0; %1/yr
Bioirrig_scale = 0.75; %1/yr
Lbottom = 30; %cm  ''Maximum depth''
t_final = 100; %year  ''Target Year''
Bioturbtop = 10;%5; %cm2/yr bioturbation coefficient at top
Bioturbbottom = 1; %cm2/yr bioturbation coefficient at bottom
bioturbscale = 3; %cm - depth scale for decrease in bioturbation

vbottom = 1; %cm/year

vbottom_fluid = 0;%0.1; %cm/year
porosbottom = 0.7; %porosity
porostop = 0.9; 
rho = 2.73; %gram/dDw
porosscale = 3; %cm - depth scsale for decrease in porosity
ageinit = 0.1; % initial age of organic matter at the sediment water interface
age_root = 1; 
k_O2 = 2; % Oxygen half-saturation constant (uM)
k_SO4 = 20; % Sulfate half-saturation constant (uM)

KFEMonod = 200;%1E7;%500; %umol/g - Monod constant for FeOOH reduction

DSO4 = 300; % cm2/yr diffusion coefficient sulfate
DCH4 = 300; % cm2/yr diffusion coefficient sulfate
DH2S = 300; % cm2/yr diffusion coefficient sulfide
DO2 = 300;   % cm2/yr diffusion coefficient oxygen
DHCO3 = 400;   % cm2/yr diffusion coefficient oxygen
DPO4 = 400; % cm2/yr diffusion coefficient 
Kreox = 500; %1600; % 1/umol/l/year ''Sulfide oxidation rate constant''
KFe_HS = 100; %1600; % 1/umol/l/year ''Iron Sulfide oxidation rate constant''
Iron_conc = 50; % Iron concentration (uM)

Calcium = 1000; %XL 10000; % calcium concentration

Calcium_activity = 0.6;%0.2; % calcium concentration
CO3_activity = 0.6;%0.028; % calcium concentration


O2init = 50; % uM - O2 concentration at SWI

SO4init = 200; %XL28000; % uM - SO4 concentration at SWI

Pinitial = 0; % uM - PO4 concentration at SWI
Feinit = 0;  % uM - Fe concentration at SWI
HSinit = 0;  % uM - H2S concentration at SWI
Ksp_ca = 3000;%1E4; %uM2

k_calcite = 1;%1E-6;%1000; %umol/l/year

% k_calcite_dis = 10;%0.05;%0.05;%0.05;%0.005; % yr-1
Mineral_Mass = 215; %basalt 215, albite 260, forsterite 140 MgO 40


DICinit = 1000;     % [DIC] [umol/kg]


HCO3init = 950; 


CH4init = 0;      % [CH4] at SWI uM
P_C_ratio = 0.0094;
K_CH4_SO4 = 100; % AOM with sulfate half saturation (uM)
K_CH4_O2 = 1; % Aerobic methane oxidation half saturation (uM)

k_AOM = 100; % AOM rate constant (1/year)

k_aerobic_CH4 = 100; % aerobic methane oxidation rate constant (1/year)


Sed_rate = 1;  %sediment accumulation rate (gram/cm2/year)
BE = 0.2; %0.1;     % burial efficiency of organic from the water column model

NPP = 800; %200      % Net Primary Production (gram/m2/year)

F_FeOx = 100;   %mmol/m2/d


kFeOx = 10; %100; % 1/umol/l/year
kFeS = 10;%10;%0.1;%0.01;%0.08;%0.2;
kapatite = 0.05; %0.01-0.1 1/year
K_HS  = 7;
X_apa = 0.04; %0.02-0.033 Van Capellen and Berner 1988
T_apa = 25; %temperature for apatite
n_apa = 1;
Kp_ads = 6E-7;
S_ads_coef = 1E4; %umol/gr
Kviv = 3*1E6; %(umol/L)5
alpha_viv = 1.5; 
kviv = 1.7E-22;
KFeS = 2500;  %umol/l
F_CaCO3 = 10;  %XL 2000; %100; %10;%500;  % Flux of CaCO3 to sediment gram/m2/year
k_Fe_pyrite = 1000;
n_power_CaCO31 = 1.76;
n_power_CaCO32 = 0.11;
n_power_CaCO33 = 4;
k_calcite_dis1 = 0.005;
k_calcite_dis2 = 10;
Q10 = 2; 
T_ref = 25;
T_future = 30;

DOC_root_1 = 0; %500;%50;%500;;%1000;    % DOC flux release in seagrass root zone (mmol/m2/day; Eldridge & MorserMarine 2000)
O2_root_1  = 0; %500;%2;%100;%500;   % O2 flux release in seagrass root zone (mmol/m2/day; Eldridge & MorserMarine 2000)
POC_root_1 = 0; %2;%0.2;%2; %0.08;  % POC flux release in seagrass root zone (mol/m2/day; Eldridge & MorserMarine 2000)

% --------------- Calculating initial depth profiles for input parameters ----------------

% k_sed = k_sed / 100;
% k_AOM = k_AOM / 10;
% k_aerobic_CH4 = k_aerobic_CH4 / 5;
% k_calcite = k_calcite / 10;
% kFeOx = kFeOx / 5;

n=101;
MaxDepth = Lbottom;
z_sed=linspace(0,MaxDepth,n);
z_biodiff=linspace(0,MaxDepth,10001);
dz_sed=MaxDepth/(n-1);
ALK = zeros(1,n);
 
poros = porosbottom + (porostop-porosbottom)*exp(-z_sed/porosscale);
Bioturb_1 = Bioturbbottom + (Bioturbtop-Bioturbbottom)*exp(-z_biodiff/bioturbscale);
Bioturb = interp1(z_biodiff,Bioturb_1,z_sed); 

Alpha_Bioirrig = Bioirrig_bottom + (Bioirrig_top-Bioirrig_bottom)*exp(-z_sed/Bioirrig_scale);
v_burial = vbottom*(1-porosbottom)./(1-poros); %cm/year
v_burial_Fluid = vbottom_fluid*(1+porosbottom)./(1+poros); %cm/year


age = ageinit + cumsum(dz_sed./v_burial);
k_sed = 10.^(-0.95*log10(age) - 0.81);



Temp_factor = Q10.^((T_future - T_ref)/10);

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

% ---------------------------- OXYGEN PENTRATION DEPTH --------------------

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
% FeooH =  FeOx;
FeooH =  zeros(1,n);

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
Iron_C(iteration,:) = C_Fe;

for i=1:n
  if Iron_C(iteration,i) < 0
      Iron_C(iteration,i) = 0;
  end 
end

F_diff_Fe(1,count_loop) = DH2S.*((C_Fe(1,2) - C_Fe(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr

% ------------------------ IRON(III) ---------------------------------------

Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
R_iron(count_loop,:) = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9; %rate of iron reduction umol/l/year
R_FeOx = (kFeOx.*C_Fe.*Oxygen);
R_FeOx_1(count_loop,:) = R_FeOx;
% Fe_3_init  = 365.*1E2.*(F_FeOx)./(v_burial(1));  %umol/l
Fe_3_init = 36.5.*(F_FeOx).*(poros(1)/(1-poros(1)))/(v_burial(1))/rho;  %umol/l

KFEMonod = 2000;

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
sigma_FeS_1 =  (Iron_C(iteration,:).*HS_conc)./((10.^(6-pH)).*KFeS);
delta_FeS  = (sigma_FeS_1 - 1);

for i=1:n
  if delta_FeS (1,i) > 0
      delta_FeS1(1,i) = 1;
  else 
      delta_FeS1(1,i) = 0;
  end 
end

R_FeS = kFeS.*C_Fe.*C_HS;

R_FeS_1(count_loop,:) = R_FeS;
R_FeS_store(iteration,:) = R_FeS;

for i=1:n
  if R_FeS(1,i) < 0
      R_FeS(1,i) = 0;
  end 
end

R_HS_Ox = Sulfide(iteration,:).* Oxygen.*Kreox; %rate of sulfide reduction umol/l/year
F_diff_HS = DH2S.*((C_HS(1,2) - C_HS(1,1))./(x(1,2)-x(1,1)))*1E-3; %umol/cm2/yr

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

C_DIC  = y_coupled(1,:);
C_alka = y_coupled(3,:);
CaCO3  = y_coupled(5,:);
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


% -------------------------------------------------------------------------
% ----------------------------- PLOTS -------------------------------------
clf;
n_plot = 6; % number of plots in each row
m_plot = 3; % number of total rows


% Organic

subplot(m_plot,n_plot,1);

plot((C_organic + POC_root).*100,z_sed,'lineWidth',2); axis ij
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

plot(Sulfate,z_sed,'lineWidth',2); axis ij
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
title('Burial Efficiency of organic matter')
ylabel('Depth (cm)');
box on
grid on

ax.LineWidth = 2;

% Aerobic respiration and sulfate reduction rates

% subplot(m_plot,n_plot,14);
% plot(R_SRR,z_sed,R_respi,z_sed,'lineWidth',2); axis ij %umol/l/year
% title('Rate (\mumol/l/year)')
% legend('Sulfate Red','Aerobic Resp');


subplot(m_plot,n_plot,14);
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


subplot(m_plot,n_plot,16);

plot(sigma_carb(end,:),z_sed,'lineWidth',2); axis ij
title('Calciite saturation (\Omega - 1)')

box on
grid on

ax.LineWidth = 2;

subplot(m_plot,n_plot,17);

plot(FeooH(1,:),z_sed,'lineWidth',2); axis ij
title('Fe(III) (\mumol/gr)')

box on
grid on

ax.LineWidth = 2;

% Mineral saturation indices

% subplot(m_plot,n_plot,16);
% plot(delta_viv1,z_sed,delta_apat2,z_sed,delta_FeS,z_sed,'lineWidth',2); axis ij  %umol/l/year
% title('Mineral saturation (\Omega - 1)')
% legend('Vivianite','Apatite','FeS');
% 
% box on
% grid on
% 
% ax.LineWidth = 2;
% 

AAA_time = toc;
AAA_data = [Oxygen' C_DIC' C_alka' pH' sigma_carb' z_sed'];

R_iron_unit = 4.*RC.*Inhib.* (FeooH./(FeooH+KFEMonod)).*1E9.*...
                       (poros./(1-poros)).*1E-3.*(1./rho); %rate of iron reduction umol/g/year
R_FeOx_unit = (kFeOx.*C_Fe.*Oxygen).*...
              (poros./(1-poros)).*1E-3.*(1./rho);

R_SRR_integ = cumsum(0.5.*R_SRR.*dz_sed.*1E-3);
R_RC_integ = cumsum(RC.*dz_sed);
R_iron_integ = cumsum(R_iron(count_loop-1,:).*dz_sed.*1E-3);
R_iron_integ_unit = cumsum(R_iron_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
R_FeOx_integ_unit = cumsum(R_FeOx_unit.*dz_sed.*rho.*((1-poros)./poros)); %umol/cm2/yr
R_FeOx_integ = cumsum(R_FeOx_1(count_loop-1,:).*dz_sed.*1E-3);
R_FeS_integ = cumsum(R_FeS.*dz_sed.*1E-3);
R_FeS_1_integ = cumsum(R_FeS_1(count_loop-1,:).*dz_sed.*1E-3);
R_HSOX_integ = cumsum(R_HS_Ox.*dz_sed.*1E-3);
R_biorrig_integ = cumsum((Alpha_Bioirrig.*(HSinit-C_HS)).*dz_sed.*1E-3);
R_biorrig_integ_iron = cumsum((Alpha_Bioirrig.*(Feinit-C_Fe)).*dz_sed.*1E-3);
R_biorrigALK_integ = cumsum((Alpha_Bioirrig.*(HCO3init-ALK)).*dz_sed.*1E-3);
R_ALK = RC.*1E9 - R_respi;
R_carb_integ = cumsum(R1_carb.*dz_sed.*1E-3);
R_ALK_integ_WITHOUT = 2.*R_FeS_integ;
R_ALK_integ_WITH = 2.*R_FeS_integ - (R_carb_integ);
R_CH4O2_integ = cumsum((k_aerobic_CH4.* CH4.* (Oxygen./(Oxygen+K_CH4_O2))).*dz_sed.*1E-3);
R_CH4SO4_integ = cumsum((k_AOM.* CH4.* (Sulfate./(Sulfate+K_CH4_SO4))).*dz_sed.*1E-3);

R_net = R_SRR_integ - R_FeS_integ - R_HSOX_integ + R_biorrig_integ - F_diff_HS;

R_net_iron = R_iron_integ - R_FeS_integ - R_FeOx_integ - F_diff_Fe(1,count_loop-1) + R_biorrig_integ_iron;

R_net_Fe3  =  (v_burial(end).*C_Fe_3(end).*rho.*((1-poros(end))./poros(end))) - ...
              (v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1))) + R_iron_integ_unit - R_FeOx_integ_unit;

R_net_Fe3_percent  =  100.*R_net_Fe3(end)./(v_burial(1).*C_Fe_3(1).*rho.*((1-poros(1))./poros(1)));

R_net_org  =  (v_burial(end).*C_organic(end).*rho.*((1-poros(end))./12)) - (v_burial(1).*C_organic(1).*rho.*((1-poros(1))./12)) + ...
              R_RC_integ;

R_net_org_percent  =  100.*R_net_org(end)./((rho.*((1-poros(1))./12)).*v_burial(1).*C_organic(1));

F_HS_tot = (R_SRR_integ - R_FeS_integ).*0.0274;

F_S_out = (R_SRR_integ - R_FeS_integ - R_HSOX_integ).*0.0274;

R_SRR_integ_store = R_SRR_integ(end).*0.0274;  %mmol/m2/d
R_ALK_integ_WITH_store = R_ALK_integ_WITH(end).*0.0274; %mmol/m2/d

F_ox_py = R_FeS_integ./R_SRR_integ;

AAA_Store_1 = [F_FeOx R_SRR_integ_store R_ALK_integ_WITH_store];
AAA_Store = [NPP.*BE R_ALK_integ_WITH(end) R_ALK_integ_WITHOUT(end) 2.*R_carb_integ(end) F_diff]; 

toc