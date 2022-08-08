function h = stripe_plot(tab, statname, masks, labels, titstr, savestr, Tpairs, varargin)
% Parameter: 
%  tab: table contain all stats
%  statname: variable / stats to plot on y axis
%  masks: cell array of masks to separate the variable into columns. 
%  labels: labels correponding to these labels. Will be plotted as
%           x tick label and legend 
%  Tpairs: a cell array of 1-by-2 arrays, indicating which pairs of columns to test between. Like {[1,2],[3,4]}
%  
% Example:
%  h = stripe_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], ...
%                    "all chan ANOVA P<1E-3", "drv_cmp", {[1,2]},'MarkerEdgeAlpha',0.3);

if nargin<7, Tpairs = {}; end
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
mean_i = mean(statscol{i});
sem_i = sem(statscol{i});
N_i = numel(statscol{i});
legstr = compose("%s %.3f(%.3f) (%d)",labels(i),mean_i,sem_i,N_i);
xjitter = 0.15*randn(N_i,1);
scatter(i+xjitter, statscol{i},'DisplayName',legstr,varargin{:});
end
xticks(1:numel(masks)); xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(statname,'interpreter','none'); title(title_str,'interpreter','none'); % xlabel(statname) 
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend('Location','best')
saveallform(figdir,compose("%s_%s_stripcmp", statname, savestr),h)
end