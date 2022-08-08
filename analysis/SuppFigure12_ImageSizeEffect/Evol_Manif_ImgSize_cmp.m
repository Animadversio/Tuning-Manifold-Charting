%% Target analysis and comparison of 1deg and 3deg images. 
Animal = "Beto";
Set_Path;
mat_dir = "O:\Mat_Statistics"; 
%% 
load(fullfile(mat_dir, Animal+"_Evol_stats.mat"))
load(fullfile(mat_dir, Animal+"_Manif_stats.mat"))
load(fullfile(mat_dir, Animal+"_Manif_stats.mat"))
sucstab = readtable(fullfile(mat_dir, Animal+"_EvolTrajStats.csv"));
% Get pairs of experiments
[Expi2cmp, prefchan_arr] = get_cmp_pairs(sucstab);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evolution Comparison 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evolution Success Ratio comparison
fprintf("1 deg Evolution: %d/9 success per P<0.01\n", sum(sucstab.t_p_initend(Expi2cmp(:,1))<0.01))
fprintf("3 deg Evolution: %d/9 success per P<0.01\n", sum(sucstab.t_p_initend(Expi2cmp(:,2))<0.01))

%% Final stage score comparison 
% Population level 
ttest2_print(sucstab.fina_act(Expi2cmp(:,1)),...
             sucstab.fina_act(Expi2cmp(:,2)),...
            "1 deg","3 deg",true) % paired comparison 
%% Process final activations of Experiments 
fina_act_mean = [];
fina_act_sem = [];
fina_act_std = [];
fina_Z_mean = [];
fina_Z_std = [];
fina_Z_sem = [];
for Expi = 1:numel(Stats)
    psth_all = cat(3,EStats(Expi).evol.psth{:});
    act_all_vec = squeeze(mean(psth_all(:,51:200,:),[1,2]));
    act_all_mean = mean(act_all_vec);
    act_all_std = std(act_all_vec);
    
    psth_mat = cat(3,EStats(Expi).evol.psth{end-2:end-1});
    act_vec = squeeze(mean(psth_mat(:,51:200,:),[1,2]));
    bsl_vec = squeeze(mean(psth_mat(:,1:50,:),[1,2]));
    fina_act_mean(Expi) = mean(act_vec) - mean(bsl_vec);
    fina_act_std(Expi)  = std(act_vec);
    fina_act_sem(Expi)  = sem(act_vec);
    fina_Z_mean(Expi) = mean((act_vec - act_all_mean) / act_all_std); 
    fina_Z_std(Expi)  = std(act_vec / act_all_std); 
    fina_Z_sem(Expi)  = sem(act_vec / act_all_std); 
end
%% 
figure(1);set(1,'pos',[ 680   631   374   347])
errorbar(fina_act_mean(Expi2cmp(:,1)), fina_act_mean(Expi2cmp(:,2)),...
        fina_act_sem(Expi2cmp(:,1)), fina_act_sem(Expi2cmp(:,1)),...
        fina_act_sem(Expi2cmp(:,2)), fina_act_sem(Expi2cmp(:,2)),'o')
axis equal;xlim([50,350]);ylim([50,350])
xlabel("1 deg Evol");ylabel("3 deg Evol");
title("Final Activation Comparison across Image size")
add_diagonal()
saveallform(outdir, "Evol_final_act_cmp", 1)
%% 
figure(2);set(2,'pos',[ 680   631   374   347])
errorbar(fina_Z_mean(Expi2cmp(:,1)), fina_Z_mean(Expi2cmp(:,2)),...
        fina_Z_sem(Expi2cmp(:,1)), fina_Z_sem(Expi2cmp(:,1)),...
        fina_Z_sem(Expi2cmp(:,2)), fina_Z_sem(Expi2cmp(:,2)),'o')
axis equal;%xlim([50,350]);ylim([50,350])
xlabel("1 deg Evol");ylabel("3 deg Evol");
title("Final Activation z score Comparison across Image size")
add_diagonal();xlim([0,0.85]);ylim([0,0.85])
saveallform(outdir, "Evol_final_act_zscore_cmp", 2)
%% 
ttest2_print(fina_Z_mean(Expi2cmp(:,1)),...
             fina_Z_mean(Expi2cmp(:,2)),...
            "zscore 1 deg","zscore 3 deg",true) % paired comparison 
%% std of manifold experiment
manif_actstd = arrayfun(@(Expi) std(cellfun(@(P)mean(P(1,51:200,:),'all'),...
        Stats(Expi).manif.psth{1}),0,'all'), Expi2cmp);
ttest2_print(manif_actstd(:,1),manif_actstd(:,2),...
            "manif std 1 deg","manif std 3 deg",true)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manifold Comparison 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demo of a Manifold Map 
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"),'MapVarStats')
outdir = "H:\Manuscript_Manifold\FigureS3C";
Manif_Map_show_fun(MapVarStats, Animal, 44, [44])
set(3,'pos',[976   128   500   520])
saveallform(outdir,"manif_map_Beto_Exp44_ch43_1deg",3)
Manif_Map_show_fun(MapVarStats, Animal, 45, [45], 4)
set(4,'pos',[976   128   500   520])
saveallform(outdir,"manif_map_Beto_Exp45_ch43_3deg",4)

%% Load Kent fitting statistics 
poptabdir = "O:\Manif_Fitting\popstats";
poptab = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));%
% Create the masks
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
driver_fittab = poptab(drivermsk, :);
driverPC23_fittab = driver_fittab(driver_fittab.space==1,:);
assert(size(driverPC23_fittab,1) == 45) % 45 Beto Exp with PC23 space. 
fittab_B = driverPC23_fittab;
%%
ttest2_print(fittab_B.R2(Expi2cmp(:,1)),fittab_B.R2(Expi2cmp(:,2)),...
            "R2 1 deg","3 deg",true); % paired comparison 
ttest2_print(fittab_B.kappa(Expi2cmp(:,1)),fittab_B.kappa(Expi2cmp(:,2)),...
            "Kappa 1 deg","3 deg",true); % paired comparison 
ttest2_print(fittab_B.beta(Expi2cmp(:,1)),fittab_B.beta(Expi2cmp(:,2)),...
            "Beta 1 deg","3 deg",true); % paired comparison 
fprintf("R2 table: 1deg 3deg\n")
fittab_B.R2(Expi2cmp)
fprintf("kappa table: 1deg 3deg\n")
fittab_B.kappa(Expi2cmp)
%% Get the Monkey B's explained variance and normalized R2 data here.
popfitdir = "O:\Manif_Fitting\popstats";
expableVar = [];
expvartab = readtable(fullfile(popfitdir,"Both_Driver_ExpVar.csv"));
expvartab(expvartab.Animal=="Beto" & expvartab.space==1);
BPC12tab = expvartab(expvartab.Animal=="Beto" & expvartab.space==1,:);
%%
BPC12tab.expvar(Expi2cmp)
BPC12tab.normR2(Expi2cmp)
%% Plot the explained variance against the explanable variance. 
figure(10);clf;hold on
scatter(fittab_B.R2(Expi2cmp(:,1)),BPC12tab.expvar(Expi2cmp(:,1)),"Display","1 deg Map");
scatter(fittab_B.R2(Expi2cmp(:,2)),BPC12tab.expvar(Expi2cmp(:,2)),"Display","3 deg Map")
xlabel("Explained Var (R2)")
ylabel("Explainable Var (noise ceiling)")
axis equal;xlim([0,1]);ylim([0,1])
add_diagonal()
legend('location','best')
title(["Comparison of Kent Fitting R2 and"," Noise Ceiling of 1deg 3deg evolutions"])
saveallform(outdir,"manif_R2_expvar_imsize_cmp",10)
%%
msk = (fittab_B.R2(Expi2cmp(:,1))>0.4)&(fittab_B.R2(Expi2cmp(:,2))>0.4);
paired_stripe_plot({fittab_B.kappa(Expi2cmp(:,1)),fittab_B.kappa(Expi2cmp(:,2))},...
            ["Kappa 1 deg","Kappa 3 deg"],{msk},["fit succes"])
%%

%% Load in non parametric statistics 
nonpardir = "O:\Manif_NonParam\summary";
NonParTab = readtable(fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"));
NonParTab = NonParTab(NonParTab.Animal=="Beto",:);
drivermsk = zeros(size(NonParTab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats(NonParTab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (NonParTab.unitnum(i) == driver_unit) & (NonParTab.chan(i) == NonParTab.prefchan(i));
end
driver_nonpartab = NonParTab(drivermsk, :);
driverPC23_nonpartab = driver_nonpartab(driver_nonpartab.space==1,:);
assert(size(driverPC23_nonpartab,1) == 45) % 45 Beto Exp with PC23 space. 
nonpartab_B = driverPC23_nonpartab;
%%
[~,~,sumstr] = ttest2_print(nonpartab_B.normAUS_bsl(Expi2cmp(:,1)),...
             nonpartab_B.normAUS_bsl(Expi2cmp(:,2)),...
            "normVUS bsl 1 deg","normVUS bsl 3 deg",true);
paired_stripe_plot({nonpartab_B.normAUS_bsl(Expi2cmp(:,1)),...
                    nonpartab_B.normAUS_bsl(Expi2cmp(:,2))},...
            ["normVUS bsl 1 deg","normVUS bsl 3 deg"],{},["fit succes"])
figh = gcf; set(figh,'pos',[ 1000         566         306         416])
title(["NormVUS comparison",sumstr]);xlim([0.5,2.5])
ylabel("normVUS")
saveallform(outdir,"normVUS_imsize_cmp",figh)


%%  Radial tuning curve AUC 
summarydir = "O:\ImMetricTuning\summary";
tab = readtable(fullfile(summarydir,"Beto"+"_RadialTuningStatsTab_squ.csv"));
%%
normAUC_cmp = tab.normAUC_mani(Expi2cmp);
ttest2_print(normAUC_cmp(:,1),normAUC_cmp(:,2),...
            "normAUC bsl 1 deg","normAUC bsl 3 deg",true)
paired_stripe_plot({normAUC_cmp(:,1), normAUC_cmp(:,2)},...
            ["normAUC bsl 1 deg","normAUC bsl 3 deg"],{},["fit succes"])
%%
ttest2_print(tab.normAUC_gab(Expi2cmp(:,1)),tab.normAUC_gab(Expi2cmp(:,2)),...
    "normAUC gabor bsl 1 deg","normAUC gabor bsl 3 deg",true)
%%

normAUC_cmp = tab.corr_mani(Expi2cmp);
ttest2_print(normAUC_cmp(:,1),normAUC_cmp(:,2),...
            "Manif corr 1 deg","Manif corr 3 deg",true)


%% Visualize pairs of Evolution and Manifold pairs. 
for ipair = 1:9
vis_EvolManif_pair(Expi2cmp(ipair,:), EStats, Stats, 1);
fprintf("%.3f  %.3f \n",normAUC_cmp(ipair,1),normAUC_cmp(ipair,2))
pause
end

%% Experiment on that day of Beto Manifold Exp 44,45
load('N:\Data-Ephys-MAT\Beto-05122019-001_formatted.mat')
RFStats = RF_Calc_Stats_fun({meta}, {rasters}, {Trials});
%% Experiment of the previous day
load('N:\Data-Ephys-MAT\Beto64chan-04122019-003_formatted.mat')
RFStats = RF_Calc_Stats_fun({meta}, {rasters}, {Trials});

%%
figdir = "E:\OneDrive - Harvard University\Manuscript_Manifold\FigureS3C";
% figh = RF_contour_plot(RFStats, "edge");
figure(2)
RFStats_indiv_chan_plot(RFStats,44)
saveallform(figdir, "Beto_Exp44-45_chan43_RFMapa", 2);
%%
S = RFStats_indiv_chan_gen_mask(RFStats,44,-1.5:0.05:1.5,-1.5:0.05:1.5);
% saveallform(figdir, "Beto_Exp44-45_chan43_RFMap_small", 2);
%%
figure()
imagesc(S(44).Xq,S(44).Yq,S(44).pred_rfmat{1}) % Gaussian Process regression. 
axis image; set(gca,'ydir','normal')
colormap("gray")
%%
figh = figure('pos',[801   579   335   280]);
imagesc(S(44).Xq,S(44).Yq,S(44).interp_rfmat{1}) % spline interpolation
axis image; set(gca,'ydir','normal')
maxval = max(S(44).interp_rfmat{1},[],'all');
caxis([maxval/20,maxval])
colormap("gray")
saveallform(figdir,"rfmap_spline_half_active",figh)



function vis_EvolManif_pair(idxs, EStats, Stats, fignum)
if nargin == 3
figure('pos',[680   444   670   530])
else
figh = figure(fignum);
set(figh,'pos',[680   444   670   530])
end
T = tiledlayout(2,2,'TileSp','compact','padding','compact');
psthaxs = {};
mapaxs = {};
for i = 1:2
Expi = idxs(i);
actcol = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),...
                    EStats(Expi).evol.psth,'uni',0);
actmean = cellfun(@mean, actcol);
actstd = cellfun(@sem, actcol);
refactcol = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),...
                    EStats(Expi).ref.psth,'uni',0);
refactmean = cellfun(@mean, refactcol);
refactstd = cellfun(@sem, refactcol);
psthaxs{i} = nexttile(T,i);
shadedErrorBar(1:numel(actmean)-1,actmean(1:end-1),actstd(1:end-1))
shadedErrorBar(1:numel(refactmean)-1,refactmean(1:end-1),refactstd(1:end-1),...
    'lineprops',{'Color','g'})
axis tight
mapaxs{i} = nexttile(T,i+2);
ui = 1;
actmap = cellfun(@(P)squeeze(mean(P(ui,51:200,:),[1,2,3])),...
                    Stats(Expi).manif.psth{1},'uni',1);
imagesc(actmap)
colorbar();
axis image
end
AlignAxisLimits(psthaxs);
AlignAxisCLimits(mapaxs);
end

function [Expi2cmp,prefchan_arr] = get_cmp_pairs(sucstab)
% Utils for getting pairs of 1deg 3deg experiments paired in 9-by-2 array 
Expi2cmp = [];
prefchan_arr = [];
prefchanlast = 0;
rowi = 0;
for Expi = 28:45
    if ~(sucstab.pref_chan(Expi) == prefchanlast)
        rowi = rowi + 1;
        prefchanlast = sucstab.pref_chan(Expi);
        prefchan_arr(end+1) = prefchanlast;
    end
    if sucstab.imgsize(Expi) == 1
        Expi2cmp(rowi,1) = Expi;
    elseif sucstab.imgsize(Expi) == 3
        Expi2cmp(rowi,2) = Expi;
    else
        error("")
    end
end
assert(all(size(Expi2cmp)==[9,2]))
assert(length(prefchan_arr)==9)
end