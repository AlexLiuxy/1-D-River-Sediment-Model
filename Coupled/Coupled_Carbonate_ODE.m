function dYdx = Coupled_Carbonate_ODE(x, Y)
% State order:
% Y(1) = DIC,   Y(2) = DIC flux variable
% Y(3) = ALK,   Y(4) = ALK flux variable
% Y(5) = CaCO3, Y(6) = CaCO3 dummy flux variable

global DHCO3 DICinit HCO3init
global Alpha_Bioirrig v_burial_Fluid v_burial z_sed poros rho
global k_calcite k_calcite_dis1 k_calcite_dis2
global n_power_CaCO31 n_power_CaCO32 n_power_CaCO33
global Calcium Calcium_activity CO3_activity Ksp_ca
global T_future Salinity
global R_DIC_prod R_ALK_prod

v_burial_f = double(interp1(z_sed, v_burial_Fluid, x, 'linear', 'extrap'));
v_burial_s = double(interp1(z_sed, v_burial,       x, 'linear', 'extrap'));
fi         = double(interp1(z_sed, poros,          x, 'linear', 'extrap'));
Alpha      = double(interp1(z_sed, Alpha_Bioirrig, x, 'linear', 'extrap'));

R_DIC_prod_1 = double(interp1(z_sed, R_DIC_prod, x, 'linear', 'extrap'));
R_ALK_prod_1 = double(interp1(z_sed, R_ALK_prod, x, 'linear', 'extrap'));

DIC   = max(real(Y(1)), 1e-12);
ALK   = max(real(Y(3)), 1e-12);
CaCO3 = max(real(Y(5)), 0);

[~, CO3_current, ~] = River_Carbonate(ALK, DIC, T_future, Salinity, 1);
CO3_current = max(real(CO3_current), 1e-12);

sigma_carb = (Calcium_activity * Calcium * CO3_activity * CO3_current) / Ksp_ca - 1;

unit_conversion = 1 ./ (1E3 * 1E6 * 1E-2 * rho * max(1 - fi, 1e-6));

if sigma_carb > 0
    R_carb_form = abs(sigma_carb)^n_power_CaCO31 * k_calcite * unit_conversion;
else
    R_carb_form = 0;
end

if (sigma_carb < 0) && (sigma_carb > -0.2)
    R_carb_disso = abs(sigma_carb)^n_power_CaCO32 * k_calcite_dis1 * CaCO3;
elseif sigma_carb <= -0.2
    R_carb_disso = abs(sigma_carb)^n_power_CaCO33 * k_calcite_dis2 * CaCO3;
else
    R_carb_disso = 0;
end

R1_carb_total = R_carb_form - R_carb_disso * (1E3 * 1E6 * 1E-2 * rho * max(1 - fi, 1e-6));

Advection_DIC = v_burial_f .* (Y(2) / (fi * DHCO3));
NR_DIC = Advection_DIC ...
       - R_DIC_prod_1 ...
       + R1_carb_total ...
       - Alpha .* (DICinit - DIC);

Advection_ALK = v_burial_f .* (Y(4) / (fi * DHCO3));
NR_ALK = Advection_ALK ...
       - R_ALK_prod_1 ...
       + 2 * R1_carb_total ...
       - Alpha .* (HCO3init - ALK);

NR_CaCO3 = R_carb_form - R_carb_disso;

dYdx = zeros(6,1);
dYdx(1) = Y(2) / (fi * DHCO3);
dYdx(2) = NR_DIC;
dYdx(3) = Y(4) / (fi * DHCO3);
dYdx(4) = NR_ALK;
dYdx(5) = NR_CaCO3 / max(v_burial_s, 1e-8);
dYdx(6) = 0;
end