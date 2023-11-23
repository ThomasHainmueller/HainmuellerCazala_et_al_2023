function [h] = violin_figure(GFP,SEP,label)
% Make paper ready violin plot figure with overlay dots
group_labels = {'GFP', 'SEP-GluA2'};
color_order = [ 216 156 198; 159 207 125]/255;

GFPx = ones(length(GFP),1)+(rand(length(GFP),1)-.5)/2;
SEPx = 2*ones(length(SEP),1)+(rand(length(SEP),1)-.5)/2;

h = figure; hold on;
set(gcf, 'Position',  [100, 100, 400, 400])

scatter(GFPx,GFP,4,'k','filled');
scatter(SEPx,SEP,4,'k','filled');
violin({GFP,SEP},'facecolor',color_order);
legend(gca,'off');

ax = gca; % current axes
ax.FontName='Arial';
ax.FontSize = 12;
box off;
set(gca,'TickLength',[0.015, 0])

xlim([0.3 2+0.7])
xticks([1:2])
xticklabels(group_labels)
xtickangle(45)

ylabel(label,'FontSize',14)
set(gca, 'Tickdir', 'out')
set(gca,'FontName', 'Arial')

end