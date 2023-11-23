function [h, data] = plot_table_vars(table, grouping, vars_to_plot, var_labels, groups_to_plot, group_labels, color_order, output_file_prefix, bPlotScatter)
% Plot variables in table based on grouping
%
% Inputs
% table: samples x variables
% vars_to_plot: columns in table to plot (default: all columns)
% var_labels: Y-axis labels (default: table.Properties.VariableNames)
% group_labels: X-axis group labels (default: group number)
% color_order: colormap for each bar graph
% output_file_prefix: prefix for eps file saving
% bPlotScatter: plot scattered single data points (default: true)
%
% Outputs
% data{vars_to_plot,grouping}: grouped data
% Figure (h)
% Eps file
%
% Ingie Hong, Johns Hopkins Medical Institute, 2019

if nargin < 3 || isempty(vars_to_plot); vars_to_plot = 1:size(table,2); end
if nargin < 4 || isempty(var_labels); var_labels = strrep(table.Properties.VariableNames, '_', ' ');  end
if nargin < 5 || isempty(groups_to_plot); groups_to_plot = 1:size(unique(grouping),1); end
if nargin < 6 || isempty(group_labels); group_labels = num2str(groups_to_plot); end
if nargin < 7 || isempty(color_order); color_order = get(gca,'colororder'); color_order=[color_order; color_order]; end
if nargin < 9 || isempty(bPlotScatter); bPlotScatter = true; end

ExpGroup = findgroups(grouping);

for i=vars_to_plot
    figure; hold on
    set(gcf, 'Position',  [100, 100, 500, 500])
    target_gene_mean=splitapply(@nanmean, table(:,i), ExpGroup);
    target_gene_std=splitapply(@nanstd, table(:,i), ExpGroup);
    target_gene_raw=table2array(table(:,i));
    target_gene_count=splitapply(@nnz, ExpGroup , ExpGroup);
    for j=1:size(groups_to_plot,2) %  1:size(target_gene_mean,1)
        h(j)=bar(j,target_gene_mean(groups_to_plot(j)), 'FaceColor', color_order(j,:), 'FaceAlpha', 1 ); 
        if bPlotScatter
            s(j)=scatter(j+2*(rand(nnz(ExpGroup==groups_to_plot(j)),1) -1/2).^3, target_gene_raw(ExpGroup==groups_to_plot(j)), 2, '.', 'k');
        end
        data{i,j}=target_gene_raw(ExpGroup==groups_to_plot(j));
    end
    errorbar(1:size(groups_to_plot,2),target_gene_mean(groups_to_plot),zeros(size(groups_to_plot,2),1), target_gene_std(groups_to_plot)./(target_gene_count(groups_to_plot).^0.5),'.k')
    xlim([0.3 size(groups_to_plot,2)+0.7])
    xticks([1:size(groups_to_plot,2)])
    xticklabels({group_labels{groups_to_plot}})
    xtickangle(45)
        
    %ylim([0 prctile(target_gene_raw,99.5)])
    %xlabel('Experimental Group')
    ylabel(var_labels{i})
    set(gca, 'Tickdir', 'out')
    set(gca,'FontName', 'Arial')
    
    if size(groups_to_plot,2)==2 
        % Statistical test between the two groups
        disp(num2str([lillietest(target_gene_raw(ExpGroup==1)) lillietest(target_gene_raw(ExpGroup==2)) ]))
        if lillietest(target_gene_raw(ExpGroup==1)) || lillietest(target_gene_raw(ExpGroup==2)) % If any group shows non-normal data, 
            p = ranksum( target_gene_raw(ExpGroup==groups_to_plot(1)) ,target_gene_raw(ExpGroup==groups_to_plot(2)));
            disp([table.Properties.VariableNames{i} ' p-value: ' num2str(p) '   Wilcoxson ranksum test'])
        else  % If both groups show normal data, 
            [h,p,ci,stats] = ttest2( target_gene_raw(ExpGroup==1) ,target_gene_raw(ExpGroup==2));
            disp([table.Properties.VariableNames{i} ' p-value: ' num2str(p) '   Two sample t-test'])
        end
        
        title([strrep(table.Properties.VariableNames{i}, '_', ' ')  ' - p-value: ' num2str(p) ])
    elseif size(groups_to_plot,2)>2
        % Statistical test between the three or more groups
        bAllNormal = true;
        for j=1:size(groups_to_plot,2)
            try
                bAllNormal = bAllNormal * ~lillietest(target_gene_raw(ExpGroup==groups_to_plot(j)));
                disp(['Normality of group ' num2str(j) ':  ' num2str(lillietest(target_gene_raw(ExpGroup==groups_to_plot(j)))) ])
            catch
            end
        end
        if bAllNormal  % If all groups show normal data, 
            [p,tbl,stats] = anova1( target_gene_raw(any(ExpGroup==groups_to_plot,2)), ExpGroup(any(ExpGroup==groups_to_plot,2)), 'off');
            disp([table.Properties.VariableNames{i} ' p-value: ' num2str(p) '   1-way ANOVA'])
        else           % If any group shows non-normal data, 
            [p,tbl,stats] = kruskalwallis( target_gene_raw(any(ExpGroup==groups_to_plot,2)), ExpGroup(any(ExpGroup==groups_to_plot,2)), 'off');
            disp([table.Properties.VariableNames{i} ' p-value: ' num2str(p) '   Kruskal-Wallis test (non-parametric 1-way ANOVA)'])
        end     
    end
    
    drawnow
    if nargin >= 8 && ~isempty(output_file_prefix)
        print(gcf,'-depsc','-painters',[output_file_prefix num2str(i) '.eps'])
        disp(['Saved vector file as: ' output_file_prefix num2str(i) '.eps'])
    end

end

