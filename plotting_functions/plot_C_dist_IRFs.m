nperiods_plot = 20;

fi1 = figure('Position', get(0, 'Screensize'));

% B10 Consumption
if spec_settings.plot_C_dist_diffs
    subplot(2,3,1)
else
    subplot(1,3,1)
end
hold on
plot(1:nperiods_plot,IRFs_taylor.IRFs_C_10_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'red')
plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_10_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', [0.2 0.6 0.6 0.6])
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
grid on
grid minor
title('C-10','FontSize',12,'FontWeight','Bold')
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')



%Q3 Consumption
if spec_settings.plot_C_dist_diffs
    subplot(2,3,2)
else
    subplot(1,3,2)
end
hold on
plot(1:nperiods_plot,IRFs_taylor.IRFs_C_Q3_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'red')
plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_Q3_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', [0.2 0.6 0.6 0.6])
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
title('C-Q3','FontSize',12,'FontWeight','Bold')
l = legend('Baseline', strcat('High ', texlabel('phi_pi')));
l.Location = 'SouthEast';
l.FontSize = 10;
l.FontWeight = 'Bold';
grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')



%Top 10 Consumption
if spec_settings.plot_C_dist_diffs
    subplot(2,3,3)
else
    subplot(1,3,3)
end
hold on
plot(1:nperiods_plot,IRFs_taylor.IRFs_C_90_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'red')
plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_90_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', [0.2 0.6 0.6 0.6])
plot(1:nperiods_plot,zeros(1,nperiods_plot),'--k','LineWidth',0.25)
hold off
title('C-90','FontSize',12,'FontWeight','Bold')
grid on
grid minor
set(gca, 'XTick', [0 5 10 15 20])
set(gca, 'FontSize',12)
set(gca,'FontWeight','Bold')
ylabel('%')



if spec_settings.plot_C_dist_diffs
    % B10 Consumption
    subplot(2,3,4)
    plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_10_p(1:nperiods_plot) - IRFs_taylor.IRFs_C_10_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'black')
    grid on
    grid minor
    title('C-10','FontSize',12,'FontWeight','Bold')
    set(gca, 'XTick', [0 5 10 15 20])
    set(gca, 'FontSize',12)
    set(gca,'FontWeight','Bold')
    ylabel('%')
    
    
    
    %Q3 Consumption
    subplot(2,3,5)
    plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_Q3_p(1:nperiods_plot) - IRFs_taylor.IRFs_C_Q3_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'black')
    hold off
    title('C-Q3','FontSize',12,'FontWeight','Bold')
    l = legend(strcat('High ', texlabel('phi_pi'), ' - Taylor'));
    l.Location = 'SouthEast';
    l.FontSize = 10;
    l.FontWeight = 'Bold';
    grid on
    grid minor
    set(gca, 'XTick', [0 5 10 15 20])
    set(gca, 'FontSize',12)
    set(gca,'FontWeight','Bold')
    ylabel('%')
    
    
    
    %Top 10 Consumption
    subplot(2,3,6)
    plot(1:nperiods_plot,IRFs_high_phi.IRFs_C_90_p(1:nperiods_plot) - IRFs_taylor.IRFs_C_90_p(1:nperiods_plot), 'LineWidth',1.5,'MarkerSize',7.5, 'Color', 'black')
    hold off
    title('C-90','FontSize',12,'FontWeight','Bold')
    grid on
    grid minor
    set(gca, 'XTick', [0 5 10 15 20])
    set(gca, 'FontSize',12)
    set(gca,'FontWeight','Bold')
    ylabel('% Difference ')

end

sgtitle([spec_settings.shock_names{lllllll}])

fig_name = sprintf('%s_C_Dist%s.png', spec_settings.fig_label{lllllll}, spec_settings.str_addl);

if ~exist(strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder), 'dir')
    mkdir(strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder))
end


exportgraphics(fi1, strcat(spec_settings.pltfolder, spec_settings.IRFpltfolder,fig_name), 'Resolution',250)


