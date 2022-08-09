% This script is the final summary figure maker for population analysis of
% Tuning Smoothness of manifold experiment. 
%%
global sumdir
Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_MapSmooth\popstats";
sumdir = "O:\Manif_MapSmooth\summary";
mkdir(sumdir)
% load(fullfile(Matdir, "Beto_ManifPopDynamics.mat"))
% load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
% load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
Animal="Alfa"; 
alfa_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
Animal="Beto"; 
beto_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
StatsTab_sum = [alfa_StatsTab_sum;beto_StatsTab_sum];
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
%% Useful masks for analysis
% drivermsk = (StatsTab_sum.chan==StatsTab_sum.prefchan); % this is wrong.
drivermsk = zeros(size(StatsTab_sum,1),1); 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(StatsTab_sum.Animal{i})(StatsTab_sum.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (StatsTab_sum.unitnum(i) == driver_unit) & (StatsTab_sum.chan(i) == StatsTab_sum.prefchan(i));
end
%%
validmsk = ~(StatsTab_sum.unitnum==0) & ~(StatsTab_sum.Animal=="Alfa" & StatsTab_sum.Expi==10);
tunemsk = (StatsTab_sum.F_P<1E-5);
V1msk = (StatsTab_sum.chan<=48 & StatsTab_sum.chan>=33);
V4msk = (StatsTab_sum.chan>48);
ITmsk = (StatsTab_sum.chan<33);
Alfamsk = (StatsTab_sum.Animal=="Alfa");
Betomsk = (StatsTab_sum.Animal=="Beto");
%%
tunemsk = (StatsTab_sum.F_P<1E-3);
statcol = StatsTab_sum.D1E_t(drivermsk & tunemsk & validmsk);
%%
msk = drivermsk & tunemsk & validmsk;
msk_sg = msk & StatsTab_sum.D1E_t<0 & StatsTab_sum.D1E_P<1E-3;
statcol = StatsTab_sum.D1E_Dpr(msk_sg);
numel(statcol)
prctile(statcol,[0,100])
%%
msk_sg = msk & StatsTab_sum.TVE_t<0 & StatsTab_sum.TVE_P<1E-3;
statcol = StatsTab_sum.TVE_Dpr(msk_sg);
numel(statcol)
prctile(statcol,[0,100])
%% The distribution of D1E for pure noise.
figure;histogram(StatsTab_sum.D1E_Dpr((StatsTab_sum.F_P>1E-1)))
%%
figure;hold on
histogram(StatsTab_sum.D1E_Dpr(StatsTab_sum.unitnum==1&(StatsTab_sum.F_P<1E-3)))
histogram(StatsTab_sum.D1E_Dpr(StatsTab_sum.unitnum==2&(StatsTab_sum.F_P<1E-3)))
% hist_cmp_plot(StatsTab_sum, {drivermsk & tunemsk & validmsk, ~drivermsk & tunemsk & validmsk},...
%         ["Driver", "Non-Drivers"],"All Tune","driver_non","pdf");
%%
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk & validmsk),StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk & validmsk));
fprintf("d' Non driver < Driver, t=%.3f(P=%.1e) df=%d\n",TSTAT.tstat,P,TSTAT.df)
%%
tunemsk = (StatsTab_sum.F_P<1E-3);
hist_cmp_plot(StatsTab_sum, {drivermsk & tunemsk & validmsk, ~drivermsk & tunemsk & validmsk},...
        ["Driver", "Non-Drivers"],"All Tune P<1E-3","driver_non","pdf");
%% 
figure;
scatter(StatsTab_sum.chan(Betomsk & tunemsk), StatsTab_sum.D1E_Dpr(Betomsk & tunemsk), 36, StatsTab_sum.unitnum(Betomsk & tunemsk))
xlim([0,65])
%% Summary figures
hist_cmp_plot(StatsTab_sum, {drivermsk & tunemsk, ~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"All Tune","driver_non","pdf");
hist_cmp_plot(StatsTab_sum, {Alfamsk&drivermsk & tunemsk, Alfamsk&~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"Alfa Tune","driver_non_alfa","pdf");
hist_cmp_plot(StatsTab_sum, {Betomsk&drivermsk & tunemsk, Betomsk&~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"Beto Tune","driver_non_beto","pdf");
%%
tunemsk = (StatsTab_sum.F_P<1E-3);
h = stripe_minor_plot(StatsTab_sum, "D1E_Dpr", {drivermsk & tunemsk, ~drivermsk & tunemsk}, ["Driver", "Non-Drivers"],...
		{Alfamsk, Betomsk}, ["Alfa", "Beto"], "All Tuned", "driver_cmp_anim_sep", {[1,2]}, 'markeredgealpha',0.3)

function h = hist_cmp_plot(StatsTab_sum, masks, labels, titstr, namestr, normstr,h)
global sumdir
if nargin<7,h = figure;else, set(0,"CurrentFigure",h);end
clf;hold on 
for i = 1:numel(masks)
statcol = StatsTab_sum.D1E_Dpr(masks{i});
mean_i = mean(statcol);
sem_i = sem(statcol);
N_i = sum(masks{i});
legstr = compose("%s %.3f(%.3f) (%d)",labels(i),mean_i,sem_i,N_i);
histogram(statcol,20,'norm',normstr,'DisplayName',legstr)
end
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks);
% N2 = sum(masks{2});
ylabel(normstr)
xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
title(compose("Comparison of Gap of Dirichlet Energy for\n %s channels in %s Experiments\n"+...
        "t=%.2f (p=%.1e (%d))",join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
legend('Location','best')%compose("%s(%d)",labels',Ns'))%["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.png", namestr, normstr)))
saveas(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.pdf", namestr, normstr)))
savefig(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.fig", namestr, normstr)))
end

function h = stripe_minor_plot(tab, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs, varargin)
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    statcol = tab.(statname)(M);
    mean_i = mean(statcol);
	sem_i = sem(statcol);
    N_i = sum(M);
    xjitter = 0.15*randn(N_i,1);
	legstr = compose("%s %.3f(%.3f) (%d)",labels(i)+minorlabels(j),mean_i,sem_i,N_i);
    scatter(i+xjitter,statcol,'Marker',marker(j),'DisplayName',legstr,varargin{:})
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname,'interpreter','none')

statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str,'interpreter','none')
% legend(compose("%s(%d)",labelcol',Ns'),'Location','best') % ["Driver", "Non-Drivers"]
legend('Location','best')%
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end
%%
%%
% figure(3);clf;hold on 
% histogram(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),'norm','pdf')
% histogram(StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk),'norm','pdf')
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk));
% N1 = sum(drivermsk & tunemsk);
% N2 = sum(~drivermsk & tunemsk);
% ylabel("Prob Density")
% xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
% title(compose("Comparison of Gap of Dirichlet Energy for Driver\nand Non-Driver tuned channels in all Experiments\n"+...
%         "t=%.2f (p=%.1e (%d))",TSTAT.tstat,P,TSTAT.df))
% legend(["Driver", "Non-Drivers"])
% saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.png"))
% saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.pdf"))
% savefig(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.fig"))
% %%
% figure(3);hold on 
% histogram(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),'norm','count')
% histogram(StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk),'norm','count')
% ylabel("Prob Density")
% xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
% title(compose("Comparison of Gap of Dirichlet Energy for Driver\nand Non-Driver channels in all Experiments"))
% legend(["Driver", "Non-Drivers"])
% saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.png"))
% saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.pdf"))
% savefig(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.fig"))