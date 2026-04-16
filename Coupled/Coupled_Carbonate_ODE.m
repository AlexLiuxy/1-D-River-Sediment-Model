function dYdx = Coupled_Carbonate_ODE(x, Y)
% Y(1) = DIC, Y(2) = dDIC/dx flux
% Y(3) = ALK, Y(4) = dALK/dx flux
% Y(5) = CaCO3, Y(6) = dCaCO3/dx (dummy/solid flux)

global DHCO3 RC Alpha_Bioirrig DICinit HCO3init
global v_burial_Fluid v_burial z_sed poros rho 
global k_calcite k_calcite_dis1 k_calcite_dis2 n_power_CaCO31 n_power_CaCO32 n_power_CaCO33
global Calcium Calcium_activity CO3_activity Ksp_ca
global Rate_Meth T_future R_FeS
global R_HS_Ox R_SRR R_DIC_prod R_ALK_prod Salinity
    
    v_burial_f = double(interp1(z_sed, v_burial_Fluid, x));
    v_burial_s = double(interp1(z_sed, v_burial, x)); 
    fi = double(interp1(z_sed, poros, x));
    Alpha_Bioirrig_1 = double(interp1(z_sed, Alpha_Bioirrig, x));
    
R_DIC_prod_1 = double(interp1(z_sed, R_DIC_prod, x));
R_SRR1       = double(interp1(z_sed, R_SRR, x));
R_HS_Ox_1    = double(interp1(z_sed, R_HS_Ox, x));
R_FeS_1      = double(interp1(z_sed, R_FeS, x));
R_ALK_prod_1 = double(interp1(z_sed, R_ALK_prod, x));

    
    DIC = max(real(Y(1)), 1e-12);
    ALK = max(real(Y(3)), 1e-12);
    CaCO3 = max(real(Y(5)), 0);

    % calculate real time CO3
    [~, CO3_current, ~] = River_Carbonate(ALK, DIC, T_future, Salinity, 1); 

    if isempty(CO3_current) || isnan(CO3_current)
        CO3_current = 1e-12;
    end

    % calculate real time sigma_carb (Omega)
    sigma_carb = (Calcium_activity * Calcium * CO3_activity * CO3_current) / Ksp_ca - 1;

    % calculate real time CaCO3
    unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho * (1 - fi));
    
    R_carb_form = (sigma_carb > 0) * abs(sigma_carb)^n_power_CaCO31 * k_calcite * unit_conversion;
    R_carb_disso = (-0.2 < sigma_carb & sigma_carb < 0) * abs(sigma_carb)^n_power_CaCO32 * k_calcite_dis1 * CaCO3 ...
                 + (sigma_carb <= -0.2) * abs(sigma_carb)^n_power_CaCO33 * k_calcite_dis2 * CaCO3;
                 
    
    R1_carb_total = R_carb_form - R_carb_disso * (1E3 * 1E6 * 1E-2 * rho * (1 - fi));

    
    % account for DIC loss in methanogenesis
%     R_Meth_current = double(interp1(z_sed, Rate_Meth, x)); 
    
    
    Advection_DIC = v_burial_f .* (Y(2) / (fi * DHCO3));
NR_DIC = Advection_DIC ...
       - R_DIC_prod_1 ...
       + R1_carb_total ...
       - (Alpha_Bioirrig_1 * (DICinit - DIC));
    
    Advection_ALK = v_burial_f .* (Y(4) / (fi * DHCO3));
NR_ALK = Advection_ALK ...
       - R_ALK_prod_1 ...
       + 2 * R1_carb_total ...
       - (Alpha_Bioirrig_1 * (HCO3init - ALK));

    NR_CaCO3 = R_carb_form - R_carb_disso;

    % 
    % dYdx = [ Y(2) / (fi * DHCO3);
    %          NR_DIC;
    %          Y(4) / (fi * DHCO3);
    %          NR_ALK;
    %          NR_CaCO3 / v_burial_s;
    %          0 ];

    dYdx = zeros(6,1);
    
    dYdx(1) = Y(2) / (fi * DHCO3);
    dYdx(2) = NR_DIC;
    dYdx(3) = Y(4) / (fi * DHCO3);
    dYdx(4) = NR_ALK;
    dYdx(5) = NR_CaCO3 / v_burial_s;
    dYdx(6) = 0;

end