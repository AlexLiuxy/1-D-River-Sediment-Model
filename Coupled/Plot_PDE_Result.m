function Plot_PDE_Result(Result)
% PLOT_PDE_RESULT
% Basic diagnostic plot for the first PDE spin-up.

Grid = Result.Grid;
S = Result.State_final;
D = Result.Diag_final;

z = Grid.z;

figure('Name', 'FV-MOL PDE Result', 'Color', 'w');
clf;

n_plot = 6;
m_plot = 3;

subplot(m_plot,n_plot,1);
plot((S.OM_lab + S.OM_ref) .* 100, z, 'LineWidth', 2); axis ij
title('Organic (%gDW)');
ylabel('Depth (cm)');
grid on; box on

subplot(m_plot,n_plot,2);
plot(S.O2, z, 'LineWidth', 2); axis ij
title('O_2 (\muM)');
grid on; box on

subplot(m_plot,n_plot,3);
plot(S.Fe2, z, 'LineWidth', 2); axis ij
title('Fe^{2+} (\muM)');
grid on; box on

subplot(m_plot,n_plot,4);
plot(S.SO4, z, 'LineWidth', 2); axis ij
title('SO_4 (\muM)');
grid on; box on

subplot(m_plot,n_plot,5);
plot(S.HS, z, 'LineWidth', 2); axis ij
title('H_2S total (\muM)');
grid on; box on

subplot(m_plot,n_plot,6);
plot(S.CH4, z, 'LineWidth', 2); axis ij
title('CH_4 (\muM)');
grid on; box on

subplot(m_plot,n_plot,7);
plot(S.CaCO3 .* 100, z, 'LineWidth', 2); axis ij
title('CaCO_3 (%gDW)');
ylabel('Depth (cm)');
grid on; box on

subplot(m_plot,n_plot,8);
plot(S.DIC, z, 'LineWidth', 2); axis ij
title('DIC (\muM)');
grid on; box on

subplot(m_plot,n_plot,9);
plot(S.ALK, z, 'LineWidth', 2); axis ij
title('ALK (\muM)');
grid on; box on

subplot(m_plot,n_plot,10);
plot(S.Ca, z, 'LineWidth', 2); axis ij
title('Ca^{2+} (\muM)');
grid on; box on

subplot(m_plot,n_plot,11);
plot(D.pH, z, 'LineWidth', 2); axis ij
title('pH');
grid on; box on

subplot(m_plot,n_plot,12);
plot(D.RC_uM, z, 'LineWidth', 2); axis ij
title('Mineralization (\muM/yr)');
grid on; box on

subplot(m_plot,n_plot,13);
plot(S.FeOOH, z, 'LineWidth', 2); axis ij
title('FeOOH (\mumol/g)');
ylabel('Depth (cm)');
grid on; box on

subplot(m_plot,n_plot,14);
plot(D.R_SRR ./ 365, z, 'LineWidth', 2); axis ij
title('SRR (\muM/d)');
grid on; box on

subplot(m_plot,n_plot,15);
plot(D.R_Meth ./ 365, z, 'LineWidth', 2); axis ij
title('Methanogenesis (\muM/d)');
grid on; box on

subplot(m_plot,n_plot,16);
plot(D.R_FeS ./ 365, z, 'LineWidth', 2); axis ij
title('FeS formation (\muM/d)');
grid on; box on

subplot(m_plot,n_plot,17);
plot(D.sigma_carb, z, 'LineWidth', 2); axis ij
title('\Omega - 1');
grid on; box on

% subplot(m_plot,n_plot,18);
% plot(D.R_carb_net_solid, z, 'LineWidth', 2); axis ij
% title('Net CaCO_3 rxn (g/g/yr)');
% grid on; box on

drawnow;
end