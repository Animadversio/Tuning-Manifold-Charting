function h=stripe_simple_plot(statscol, statname, labels, figdir, savestr, varargin)
% One array statistics stripe plot. 
% Signature:
%  h=stripe_simple_plot(statscol, statname, labels, figdir, savestr, varargin)
% 
% Parameter: 
%  statscol: an array or a cell array of multiple vectors of values
%  statname: variable / stats name to show on y axis
%  labels: labels correponding to each cell in statscol. Will be plotted as
%           x tick label and legend 
%  figdir: save dir 
%  savestr: name string to save fig
%  
% Example:
%  h = stripe_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], ...
%                    "all chan ANOVA P<1E-3", "drv_cmp", {[1,2]},'MarkerEdgeAlpha',0.3);
if nargin<=1, statname="stat";end
if nargin<=2, labels=[""]; savestr=""; figdir=""; end
if nargin<=4, savestr=""; end
% if nargin<7, Tpairs = {}; end
% global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
if ~iscell(statscol), statscol = {statscol}; end
for i = 1:numel(statscol)
mean_i = mean(statscol{i});
sem_i = sem(statscol{i});
N_i = numel(statscol{i});
legstr = compose("%s %.3f(%.3f) (%d)",labels(i),mean_i,sem_i,N_i);
xjitter = 0.15*randn(N_i,1);
scatter(i+xjitter, statscol{i},'DisplayName',legstr,varargin{:});
end
xticks(1:numel(statscol)); xticklabels(labels)
FStat = anova_cells(statscol);
% title_str = compose("Comparison of %s for\n channels %s\nANOVA F:%.3f(p=%.1e)",...
%     statname,FStat.F,FStat.F_P);
title_str = compose("Comparison of %s",statname);
for gi = 1:numel(statscol)
[~,~,sumstr,] = ttest_print(statscol{gi},labels{gi});
title_str = title_str+compose("\n")+sumstr;
% [~,P,CI,TST] = ttest(statscol{gi});
% tval = TST.tstat; pval = P;
% sumstr = sprintf("\n%s (%.1f)-0: t=%.3f(df=%d), P=%.1e, CI=[%.1f,%.1f]",labels{fi},g1_mean,tval,TST.df,pval,CI(1),CI(2));
end
% for pi = 1:numel(Tpairs)
% pair = Tpairs{pi};
% [~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
% title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
% end
Ns = cellfun(@numel,statname);
ylabel(statname,'interpreter','none'); title(title_str,'interpreter','none'); % xlabel(statname) 
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend('Location','best')
saveas(h,fullfile(figdir,compose("%s_%s_stripe.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripe.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripe.fig", statname, savestr)))

end