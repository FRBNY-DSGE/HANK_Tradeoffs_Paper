nperiods_plot = 20;

plot_fig_1 = strcmp(spec_settings.fig_label{lllllll},"MonetaryPolicy");


fi1 = figure('Position', get(0, 'Screensize'));
subplot(2,3,1)
hold on
plot(1:nperiods_plot,IRFs_taylor.IRFs_Y_p(1:nperiods_plot),'-','LineWidth',1.5,'Color','red')
if ~plot_fig_1
    plot(1:nperiods_plot,IRFs_high_phi.IRFs_Y_p(1:nperiods_plot),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
title('Output','FontSize',12,'FontWeight','Bold')

grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')

subplot(2,3,2)
hold on
plot(1:nperiods_plot,400*(IRFs_taylor.pi_IRFs_p(1:nperiods_plot)-param.pi_bar),'-','LineWidth',1.5,'Color','red')
if ~plot_fig_1
    plot(1:nperiods_plot,400*(IRFs_high_phi.pi_IRFs_p(1:nperiods_plot)-param.pi_bar),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
grid on
grid minor
if ~plot_fig_1
    if spec_settings.calculate_high_phi_pi
        l = legend('Baseline',strcat('High ', texlabel('phi_pi')));
    elseif spec_settings.calculate_high_kappa
        l = legend('Baseline',strcat('High ', texlabel('kappa')));
    end
    l.Location = 'northeast';
    l.FontSize = 10;
    l.FontWeight = 'Bold';
end
title('Inflation Rate','FontSize',12,'FontWeight','Bold')
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('pp (annual)')

subplot(2,3,3)
hold on
plot(1:nperiods_plot,400*(IRFs_taylor.Rcb_IRFs_p(2:nperiods_plot+1)-param.R_cb),'-','LineWidth',1.5,'Color','red')
if ~plot_fig_1
    plot(1:nperiods_plot,400*(IRFs_high_phi.Rcb_IRFs_p(2:nperiods_plot+1)-param.R_cb),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
grid on
grid minor
title('Nominal Rate','FontSize',12,'FontWeight','Bold')
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('pp  (annual)')

subplot(2,3,4)
hold on
plot(1:nperiods_plot,IRFs_taylor.IRFs_W_p(2:nperiods_plot+1),'-','LineWidth',1.5,'Color','red')
if ~plot_fig_1
    plot(1:nperiods_plot,IRFs_high_phi.IRFs_W_p(2:nperiods_plot+1),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
hold off
title('Wage','FontSize',12,'FontWeight','Bold')
grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')

subplot(2,3,5)
hold on
plot(1:nperiods_plot,100*(IRFs_taylor.unemp_IRFs_p(1:nperiods_plot)-SS_stats.u),'-','LineWidth',1.5,'Color','red')
if ~plot_fig_1
    plot(1:nperiods_plot,100*(IRFs_high_phi.unemp_IRFs_p(1:nperiods_plot)-SS_stats.u),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
title('Unemployment rate','FontSize',12,'FontWeight','Bold')
grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('pp')


subplot(2,3,6)
hold on
plot(1:nperiods_plot,100*((IRFs_taylor.Profit_IRFs_p(1:nperiods_plot)+param.fix2)/(SS_stats.Profit+param.fix2)-1),'-','LineWidth',1.5,'Color', 'red')
if ~plot_fig_1
    plot(1:nperiods_plot,100*((IRFs_high_phi.Profit_IRFs_p(1:nperiods_plot)+param.fix2)/(SS_stats.Profit+param.fix2)-1),'-','LineWidth',1.5,'Color',[0.2 0.6 0.6 0.6])
end
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
title('Profit','FontSize',12,'FontWeight','Bold')
grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')

sgtitle([spec_settings.shock_names{lllllll}])

fig_name = sprintf('%s%s.png',spec_settings.fig_label{lllllll}, spec_settings.str_addl);
if ~exist(strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder), 'dir')
    mkdir(strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder))
end
exportgraphics(fi1, strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder,fig_name), 'Resolution',250)



