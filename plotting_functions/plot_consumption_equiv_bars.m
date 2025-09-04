close all



nperiods_plot = 20;


% 
% figure('Position', get(0, 'Screensize'))
% hold on
% grid on
% c = {'B0.1','B1','B10','Q1','Q2','Q3','Q4','Q5','T10','T1','T0.1'};
% % b1 = bar(10000*CE_bar(2:end,:)','stacked','LineWidth',1.2);
% b1  = bar(10000*(CE_bar_compare(2:end,:)'-CE_bar_base(2:end,:)'),'stacked','LineWidth',1.2);
% 
% for jj = 1:11
%     f(jj)  = plot((jj-1)+[0.7:0.3:1.3],10000*(CE_bar_compare(1,jj)-CE_bar_base(1,jj))*ones(1,3),'-k','LineWidth',5);
% end
% hold off
% grid on
% grid minor
% b1(2).FaceColor = 'y';
% b1(2).FaceAlpha = 0.7;
% b1(2).BarWidth  = 0.7;
% b1(2).LineWidth = 1;
% 
% b1(1).FaceColor = [0.2 0.2 0.8];
% b1(1).FaceAlpha = 0.4;
% b1(1).BarWidth  = 0.7;
% b1(1).LineWidth = 1;
% 
% b1(4).FaceColor = [0.6 0.6 0.1];
% b1(4).FaceAlpha = 0.5;
% b1(4).BarWidth  = 0.7;
% b1(4).LineWidth = 1;
% 
% b1(5).FaceColor = [0 1 1];
% b1(5).FaceAlpha = 0.5;
% b1(5).BarWidth  = 0.7;
% b1(5).LineWidth = 1;
% b1(6).FaceColor = [0.2 0.8 0.1];
% b1(6).FaceAlpha = 0.8;
% b1(6).BarWidth  = 0.7;
% b1(6).LineWidth = 1;
% b1(3).FaceColor = 'm';
% b1(3).FaceAlpha = 0.6;
% b1(3).BarWidth  = 0.7;
% b1(3).LineWidth = 1;
% ylabel('bp')
% lh = legend([b1(1),b1(2),b1(3),b1(4),b1(5),f(1)],' Real rate',' Profit',' Equity price',' Real wage',' Job finding rate',' Total (High - Low \phi_\pi)');
% legend('boxoff')
% legend('Location','SouthEast')
% lh.NumColumns = 2;
% xlim([0 12])
% set(gca,'XTick',1:11,'XTickLabel',c)
% set(gca, 'FontSize',16)
% set(gca, 'FontWeight','Bold')
% 
% figure('Position', get(0, 'Screensize'))
% hold on
% grid on
% c = {'B0.1','B1','B10','Q1','Q2','Q3','Q4','Q5','T10','T1','T0.1'};
% b1  = bar(10000*(CE_bar_base(2:end,:)'),'stacked','LineWidth',1.2);
% 
% for jj = 1:11
%     f(jj)  = plot((jj-1)+[0.7:0.3:1.3],10000*(CE_bar_base(1,jj))*ones(1,3),'-k','LineWidth',5);
% end
% 
% hold off
% grid on
% grid minor
% b1(2).FaceColor = 'y';
% b1(2).FaceAlpha = 0.7;
% b1(2).BarWidth  = 0.7;
% b1(2).LineWidth = 1;
% 
% b1(1).FaceColor = [0.2 0.2 0.8];
% b1(1).FaceAlpha = 0.4;
% b1(1).BarWidth  = 0.7;
% b1(1).LineWidth = 1;
% 
% b1(4).FaceColor = [0.6 0.6 0.1];
% b1(4).FaceAlpha = 0.5;
% b1(4).BarWidth  = 0.7;
% b1(4).LineWidth = 1;
% 
% b1(5).FaceColor = [0 1 1];
% b1(5).FaceAlpha = 0.5;
% b1(5).BarWidth  = 0.7;
% b1(5).LineWidth = 1;
% b1(6).FaceColor = [0.2 0.8 0.1];
% b1(6).FaceAlpha = 0.8;
% b1(6).BarWidth  = 0.7;
% b1(6).LineWidth = 1;
% b1(3).FaceColor = 'm';
% b1(3).FaceAlpha = 0.6;
% b1(3).BarWidth  = 0.7;
% b1(3).LineWidth = 1;
% 
% ylabel('bp')
% lh = legend([b1(1),b1(2),b1(3),b1(4),b1(5),f(1)],' Real rate',' Profit',' Equity price',' Real wage',' Job finding rate',' Total');
% legend('boxoff')
% legend('Location','NorthWest')
% lh.NumColumns = 2;
% xlim([0 12])
% set(gca,'XTick',1:11,'XTickLabel',c)
% set(gca, 'FontSize',16)
% set(gca, 'FontWeight','Bold')

if strcmp(spec_settings.fig_label{lllllll}, "PriceMarkup")
figure('Position', get(0, 'Screensize'))
hold on
grid on
c = {'B10','Q1','Q2','Q3','Q4','Q5','T10'};
b1  = bar(10000*(CE_bar_compare(2:end,3:end-2)'-CE_bar_base(2:end,3:end-2)'),'stacked','LineWidth',1.2);
for jj = 1:7
    f(jj)  = plot((jj-1)+[0.7:0.3:1.3],10000*(CE_bar_compare(1,jj+2)-CE_bar_base(1,jj+2))*ones(1,3),'-k','LineWidth',5);
end
hold off
grid on
grid minor
b1(2).FaceColor = 'y';
b1(2).FaceAlpha = 0.7;
b1(2).BarWidth  = 0.7;
b1(2).LineWidth = 1;

b1(1).FaceColor = [0.2 0.2 0.8];
b1(1).FaceAlpha = 0.4;
b1(1).BarWidth  = 0.7;
b1(1).LineWidth = 1;

b1(4).FaceColor = [0.6 0.6 0.1];
b1(4).FaceAlpha = 0.5;
b1(4).BarWidth  = 0.7;
b1(4).LineWidth = 1;

b1(5).FaceColor = [0 1 1];
b1(5).FaceAlpha = 0.5;
b1(5).BarWidth  = 0.7;
b1(5).LineWidth = 1;
b1(6).FaceColor = [0.2 0.8 0.1];
b1(6).FaceAlpha = 0.8;
b1(6).BarWidth  = 0.7;
b1(6).LineWidth = 1;
b1(3).FaceColor = 'm';
b1(3).FaceAlpha = 0.6;
b1(3).BarWidth  = 0.7;
b1(3).LineWidth = 1;



ylabel('bp')
lh = legend([b1(1),b1(2),b1(3),b1(4),b1(5),f(1)],' Real rate',' Profit',' Equity price',' Real wage',' Job finding rate',' Total (High - Low \phi_\pi)');
legend('boxoff')
legend('Location','SouthEast')
lh.NumColumns = 2;
xlim([0 8])
set(gca,'XTick',1:7,'XTickLabel',c)
set(gca, 'FontSize',16)
set(gca, 'FontWeight','Bold')

if ~exist(strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder), 'dir')
    mkdir(strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder));
end
fig_name = sprintf('%s_CE_diff%s.png', spec_settings.fig_label{lllllll},spec_settings.str_addl);
exportgraphics(gca, strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder,fig_name), 'Resolution',250)
end


figure('Position', get(0, 'Screensize'))
hold on
grid on
c = {'B10','Q1','Q2','Q3','Q4','Q5','T10'};
b1  = bar(10000*(CE_bar_base(2:end,3:end-2)'),'stacked','LineWidth',1.2);

for jj = 1:7
    f(jj)  = plot((jj-1)+[0.7:0.3:1.3],10000*(CE_bar_base(1,jj+2))*ones(1,3),'-k','LineWidth',5);
end
hold off
grid on
grid minor
b1(2).FaceColor = 'y';
b1(2).FaceAlpha = 0.7;
b1(2).BarWidth  = 0.7;
b1(2).LineWidth = 1;

b1(1).FaceColor = [0.2 0.2 0.8];
b1(1).FaceAlpha = 0.4;
b1(1).BarWidth  = 0.7;
b1(1).LineWidth = 1;

b1(4).FaceColor = [0.6 0.6 0.1];
b1(4).FaceAlpha = 0.5;
b1(4).BarWidth  = 0.7;
b1(4).LineWidth = 1;

b1(5).FaceColor = [0 1 1];
b1(5).FaceAlpha = 0.5;
b1(5).BarWidth  = 0.7;
b1(5).LineWidth = 1;
b1(6).FaceColor = [0.2 0.8 0.1];
b1(6).FaceAlpha = 0.8;
b1(6).BarWidth  = 0.7;
b1(6).LineWidth = 1;
b1(3).FaceColor = 'm';
b1(3).FaceAlpha = 0.6;
b1(3).BarWidth  = 0.7;
b1(3).LineWidth = 1;

ylabel('bp')


lh = legend([b1(1),b1(2),b1(3),b1(4),b1(5),f(1)],' Real rate',' Profit',' Equity price',' Real wage',' Job finding rate',' Total');
legend('boxoff')
legend('Location','NorthWest')
lh.NumColumns = 2;
xlim([0 8]) %10
set(gca,'XTick',1:7,'XTickLabel',c) %9
set(gca, 'FontSize',16)
set(gca, 'FontWeight','Bold')

fig_name = sprintf('%s_CE%s.png', spec_settings.fig_label{lllllll},spec_settings.str_addl);

if ~exist(strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder), 'dir')
    mkdir(strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder));
end
exportgraphics(gca, strcat(spec_settings.pltfolder,"/", spec_settings.CEpltfolder,fig_name), 'Resolution',250)
