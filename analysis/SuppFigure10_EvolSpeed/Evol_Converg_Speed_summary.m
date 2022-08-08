%% Manifold paper Supplementary figure compare evolution speed across areas
% the population statistics for evolution speed  and visualization of the trajectories
Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";

%% Load in pre computed data from Evol_Converg_Speed_Analysis.m
global figdir
figdir = "O:\EvolTraj_Cmp\summary";
A = load(fullfile(mat_dir,"Alfa"+"_EvolTrajStats.mat"),'EvolTrajStat');
B = load(fullfile(mat_dir,"Beto"+"_EvolTrajStats.mat"),'EvolTrajStat');
EvolTrajStat = [A.EvolTrajStat,B.EvolTrajStat];
StatTab = readtable(fullfile(mat_dir, "Both"+"_EvolTrajStats.csv"));

% Create masks from the table `StatTab`
StatTab.area = arrayfun(@area_map,StatTab.pref_chan);
sucsmsk = (StatTab.t_p_initmax<1E-2);
V1msk = StatTab.pref_chan <=48 & StatTab.pref_chan >= 33;
V4msk = StatTab.pref_chan <=64 & StatTab.pref_chan >= 49;
ITmsk = StatTab.pref_chan <=32 & StatTab.pref_chan >= 1;
Alfamsk = StatTab.Animal == "Alfa";
Betomsk = StatTab.Animal == "Beto";

%% Doing final population statistics for evolution speed.
%% initial paired ttest
[~,P,CI,TST] = ttest2(StatTab.step63(V4msk&sucsmsk), StatTab.step63(V1msk&sucsmsk));
fprintf("V4 - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step63(ITmsk&sucsmsk), StatTab.step63(V1msk&sucsmsk));
fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step63(ITmsk&sucsmsk), StatTab.step63(V4msk&sucsmsk));
fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
%%
[~,P,CI,TST] = ttest2(StatTab.step50(V4msk&sucsmsk), StatTab.step50(V1msk&sucsmsk));
fprintf("V4 - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&sucsmsk), StatTab.step50(V1msk&sucsmsk));
fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&sucsmsk), StatTab.step50(V4msk&sucsmsk));
fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));

%% plot the statistics as stripe across Area
h = stripe_plot(StatTab, "step63", {V1msk&sucsmsk, V4msk&sucsmsk, ITmsk&sucsmsk}, ["V1", "V4", "IT"], ...
    "all exps", "area_cmp", {[3,1],[2,1],[3,2]});
%% plot the statistics across Area and Animal 
h = stripe_minor_plot(StatTab, "step63", {V1msk&sucsmsk, V4msk&sucsmsk, ITmsk&sucsmsk}, ["V1", "V4", "IT"], ...
    {Alfamsk&sucsmsk, Betomsk&sucsmsk}, ["Alfa", "Beto"], "all exps", "area_anim_cmp", {[3,1],[2,1],[3,2]});

%% Combined statistical testing the areal progression 
sucsmsk = (StatTab.t_p_initmax<1E-2);
area_prog_cmp(StatTab, "step63", {V1msk&sucsmsk, V4msk&sucsmsk, ITmsk&sucsmsk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step85", {V1msk&sucsmsk, V4msk&sucsmsk, ITmsk&sucsmsk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step50", {V1msk&sucsmsk, V4msk&sucsmsk, ITmsk&sucsmsk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})

%% Customized testing with alternative masks
msk = (StatTab.t_p_initmax<1E-2);
area_prog_cmp(StatTab, "step63", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step85", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step50", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot trajectory summary across conditions

%% Newer version processing and API 
%% Creat masks
%% All trajectory excluding the last generation for noise reduction.
traj_col = arrayfun(@(E)E.act_traj_mean(1:end-1), EvolTrajStat,'uni',0);%{EvolTrajStat.act_traj_mean}
%%
summarize_trajs_merge_extrap(traj_col,"movmean_max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge_extrap");
%% Add single trajectory traces for visualization
summarize_trajs_merge_extrap(traj_col,"movmean_max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge_extrap_singtraj",true);
%% 
summarize_trajs_merge_extrap(traj_col,"movmean_max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge_extrap_singtraj_smooth",true);

%% Older version processing and API 
% all mean trajectory in an cell array
traj_col = arrayfun(@(E)E.act_traj_mean, EvolTrajStat,'uni',0);
% Add movemean to smooth the traj
summarize_trajs_tile(traj_col,"max",1,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim");
summarize_trajs_tile(traj_col,"max",3,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_movmean");
%%
summarize_trajs_tile(traj_col,"max",3,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_movmean");
%%
summarize_trajs_tile(traj_col,"max",1,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs");
summarize_trajs_tile(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean");

%% Merge trajs in the same tile for an monk 
summarize_trajs_merge(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_sucs_movmean_merge");
%%
figure(14);AlignAxisLimits(get(14,'Child'))
%%
summarize_trajs_merge(traj_col,"movmean_rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_SmthRngNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"movmean_rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_SmthRngNorm_scoreTraj_avg_Area_sucs_movmean_merge");
%%
summarize_trajs_merge(traj_col,"rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_RngNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_RngNorm_scoreTraj_avg_Area_sucs_movmean_merge");


function area_prog_cmp(StatTab, varnm, msks, labels, pairs)
% Statistical test for areal progression of a statistics 
% Doing 
%    pair-wise ttest, 
%    ANOVA F test, 
%    Spearman rank correlation 
%    linear model test of coef. 
% in a single line. 
fprintf("Comparison of %s among areas\n",varnm)
fullvec = [];
idxvec = [];
for i = 1:numel(msks)
varvec = StatTab.(varnm)(msks{i});
fprintf("%s: %.1f(%.1f,n=%d)  ",labels(i),mean(varvec),sem(varvec),numel(varvec));
fullvec = [fullvec; varvec];
idxvec = [idxvec; i * ones(numel(varvec), 1)];
end
fprintf("\n")
for pi = 1:numel(pairs)
i = pairs{pi}(1); j = pairs{pi}(2);
[~,P,CI,TST] = ttest2(StatTab.(varnm)(msks{i}), StatTab.(varnm)(msks{j}));
fprintf("%s - %s: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",labels(i),labels(j),P,TST.tstat,TST.df,CI(1),CI(2));
end
[cval, pval] = corr(idxvec, fullvec, 'Type','Spearman');
fprintf("Spearman corrlation of index and %s: %.3f (p=%.1e,n=%d)\n",varnm,cval,pval,numel(fullvec))
[F_p,F_tbl] = anova1(fullvec, idxvec, 'off');
Fval = F_tbl{2,5}; F_df = F_tbl{4,3};
anovastr = compose("ANOVA F=%.3f p=%.1e(df=%d,%d)\n",Fval,F_p,F_tbl{2,3},F_tbl{3,3});
lm = fitlm(idxvec, fullvec);
lmstr = compose("Linear Regres %s = %.3f + %s * %.3f \n Intercept %.3f+-%.3f, Slope %.3f+-%.3f\n Slope!=0: T=%.1f P=%.1e\n Fit Rsquare=%.3f\n",...
                varnm, lm.Coefficients.Estimate(1), "area", lm.Coefficients.Estimate(2),...
                lm.Coefficients.Estimate(1), lm.Coefficients.SE(1), ...
                lm.Coefficients.Estimate(2), lm.Coefficients.SE(2), ...
                lm.Coefficients.tStat(2), lm.Coefficients.pValue(2), ...
                lm.Rsquared.Ordinary);
fprintf(anovastr +lmstr+"\n")
% [~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V1msk&msk));
% fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
% [~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V4msk&msk));
% fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
end


function [h,score_m,score_s,blockvec] = summarize_trajs_tile(traj_col,norm_mode,movmean_N,major_msk,major_lab,minor_msk,minor_lab,figdir,fignm)
% msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
% anim_msks = {Alfamsk&validmsk, Betomsk&validmsk};
% label_col = ["V1", "V4", "IT"];
% anim_col = ["Alfa","Beto"];
% msk on the experimental session level. 
if nargin<2, norm_mode="max";end
if nargin<3, movmean_N=3;end
if nargin<4, major_msk={ones(size(traj_col,1),1,'logical')}; major_lab=["all"]; end
if nargin<6 || isempty(minor_msk)
    minor_msk={ones(size(traj_col,1),1,'logical')}; 
end
if nargin<7 || isempty(minor_lab), minor_lab = ["all"]; end 
if nargin<8, figdir=""; end
if nargin<9, fignm=compose("%sNorm_scoreTraj_avg_Area_Anim_movmean",norm_mode); end
Corder = colororder();
% Normalize the trajs 
if strcmp(norm_mode,"max") % normalize each traj to the max in that exp. 
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
elseif strcmp(norm_mode,"rng") % normalize each traj, min to 0, max to 1
normed_traj_col = cellfun(@(traj) (traj-min(traj))/(max(traj)-min(traj)),traj_col,'uni',0);
elseif strcmp(norm_mode,"movmean_rng") % normalize each traj to min max of smooth traj
normed_traj_col = [];
for i = 1:numel(traj_col)
smooth_traj = movmean(traj_col{i},3);
MIN = min(smooth_traj); 
MAX = max(smooth_traj);
normed_traj_col{i} = (traj_col{i}-MIN)/(MAX-MIN);
end
normed_traj_col = reshape(normed_traj_col, size(traj_col));
else
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
end
% block cell array
block_col = cellfun(@(traj)1:numel(traj),traj_col,'uni',0);
% [score_m,score_s,blockvec] = sort_scoreblock(block_col, normed_traj_col)
nrow = numel(minor_msk); ncol = numel(major_msk);
h=figure; %fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean", Animal);
set(h,'pos',[125,   258,   310*ncol, 100+310*nrow])
T = tiledlayout(nrow,ncol,"pad",'compact',"tilespac",'compact');
for animi=1:nrow
for mski=1:ncol
nexttile(T,mski+ncol*(animi-1))
msk = major_msk{mski} & minor_msk{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_col(msk),normed_traj_col(msk));
score_C_col_mov_m = movmean(score_C_col_m,movmean_N);
score_C_col_mov_s = movmean(score_C_col_s,movmean_N);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory\n%s %s (n=%d)",minor_lab(animi),major_lab(mski),sum(msk)))
% legend(["Full","50D"])
end
end
title(T, compose("Summary of mean evol trajectory (%d block movmean) for each area",movmean_N), 'FontSize',16)
saveallform(figdir,fignm);
for animi=1:nrow
for mski=1:ncol
    nexttile(T,mski+ncol*(animi-1))
    msk = major_msk{mski} & minor_msk{animi};
    xlim([1,prctile(cellfun(@numel,block_col(msk)),80)]);
end
end
saveallform(figdir,fignm+"_Xlim");
end
