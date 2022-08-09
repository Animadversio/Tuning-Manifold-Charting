%% Manif_Tune_Smooth_Exemplar
%  Find exemplar for tuning maps of a given d' or other stats. 
global sumdir
Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_MapSmooth\popstats";
sumdir = "O:\Manif_MapSmooth\summary";
mkdir(sumdir)
%% Prepare the tables and the structure.
Animal="Alfa"; 
alfa_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
Animal="Beto"; 
beto_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
StatsTab_sum = [alfa_StatsTab_sum;beto_StatsTab_sum];
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
%% Get the masks 
drivermsk = zeros(size(StatsTab_sum,1),1); 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(StatsTab_sum.Animal{i})(StatsTab_sum.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (StatsTab_sum.unitnum(i) == driver_unit) & (StatsTab_sum.chan(i) == StatsTab_sum.prefchan(i));
end
validmsk = ~(StatsTab_sum.unitnum==0) & ~(StatsTab_sum.Animal=="Alfa" & StatsTab_sum.Expi==10);
tunemsk = (StatsTab_sum.F_P<1E-5);
V1msk = (StatsTab_sum.chan<=48 & StatsTab_sum.chan>=33);
V4msk = (StatsTab_sum.chan>48);
ITmsk = (StatsTab_sum.chan<33);
Alfamsk = (StatsTab_sum.Animal=="Alfa");
Betomsk = (StatsTab_sum.Animal=="Beto");
%%
Animal="Alfa"; 
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
%% Sample tuning maps based on d' or other comparison. 
Expi=3; si=1; ui=1;
pref_chan_id = MapVarStats(Expi).units.pref_chan_id;
iCh = pref_chan_id(ui);
unitlab = MapVarStats(Expi).units.unit_name_arr(iCh);
act_col = cellfun(@(A)A(iCh,:)',MapVarStats(Expi).manif.act_col{si},'uni',0);
calc_visualize_map_smoothness(act_col,Animal,Expi,unitlab,si,sumdir);

function [h,h2]=calc_visualize_map_smoothness(act_col,Animal,Expi,unitlab,si,figdir,h,h2)
% Inputs:
%  act_col: cell array of activation vectors corresponding to trial response to each stimuli.
%  Others: meta info to label the title and savename
% Demo:
%  actmap = cellfun(@(P)mean(P(ui,51:200,:),'all'),Stats(Expi).manif.psth{si});
%  act_col = cellfun(@(P)squeeze(mean(P(ui,51:200,:),2)),Stats(Expi).manif.psth{si},'uni',0);
%  calc_visualize_map_smoothness(act_col)
%  
if nargin == 1, doplot = false;
else, doplot = true; end

spacelabels = ["PC23","PC4950","RND12"];
splab = spacelabels(si);
if nargin<=6, h=figure();h2=figure(); end
clf(h);clf(h2);
set(h, 'pos',[275         423        1086         435]);
set(h2,'pos',[386         135        1500         766]);
%% Calculate shuffled/Boostrapped statistics
actmap = cellfun(@mean,act_col);
[D2E_sp, TVE_sp, D1E_sp, D2map_sp, GPmap_sp, detJac_map] = sph_map_energy(actmap);
btrp_Scol_sp = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col); % Resample trial responses
[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp, D2map_sp_btrp, GPmap_sp_btrp, detJac_map] = sph_map_energy(actmap_btrp);
btrp_Scol_sp = [btrp_Scol_sp;[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp]];
end
shfl_Scol_sp = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl, D2map_sp_shfl, GPmap_sp_shfl, detJac_map] = sph_map_energy(actmap_shfl);
shfl_Scol_sp = [shfl_Scol_sp;[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl]];
end
[~,P_D2E_sp,~,tSt_D2E_sp] = ttest2(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
Dpr_D2E_sp = computeCohen_d(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
[~,P_TVE_sp,~,tSt_TVE_sp] = ttest2(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
Dpr_TVE_sp = computeCohen_d(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
[~,P_D1E_sp,~,tSt_D1E_sp] = ttest2(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));
Dpr_D1E_sp = computeCohen_d(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));

%% Spherical version statistics comparison Plot
set(0,'CurrentFigure',h);clf;
T = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile(T,1);hold on
histogram(btrp_Scol_sp(:,1),'Norm','pdf');
histogram(shfl_Scol_sp(:,1),'Norm','pdf');
vline(D2E_sp)
title(compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E_sp.tstat,P_D2E_sp,Dpr_D2E_sp))
nexttile(T,2);hold on
histogram(btrp_Scol_sp(:,2),'Norm','pdf');
histogram(shfl_Scol_sp(:,2),'Norm','pdf');
vline(TVE_sp)
title(compose("Total Variation Energy\n t=%.2f(%.1e) d'=%.2f",tSt_TVE_sp.tstat,P_TVE_sp,Dpr_TVE_sp))
nexttile(T,3);hold on
histogram(btrp_Scol_sp(:,3),'Norm','pdf');
histogram(shfl_Scol_sp(:,3),'Norm','pdf');
vline(D1E_sp)
title(compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E_sp.tstat,P_D1E_sp,Dpr_D1E_sp))
legend(["Trial Resampled", "Matrix Shuffled"])
if doplot
title(T,["Map Smoothness Energy (Sphere) Compared",compose("%s Exp %d Space%d(%s) Unit %s",Animal,Expi,si,splab,unitlab)])
saveallform(figdir, compose("Smooth_DE_sph_hist_cmp_%s_Exp%d_%s_%s",Animal,Expi,splab,unitlab), h)
end
%% Spherical version map comparison Plot
xarr = -90:18:90;
yarr = -90:18:90; % x ticks for the map
set(0,'CurrentFigure',h2);clf;
T = tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(T,1);
imagesc(xarr,yarr,actmap)
axis image;colorbar
title("Neural Activation")
ylabel("Trial Resampling")
xlabel("PC3 Phi (deg)")
nexttile(T,2);
imagesc(xarr,yarr,GPmap_sp_btrp)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi (deg)")
nexttile(T,3);
imagesc(xarr,yarr,GPmap_sp_btrp.*detJac_map)
axis image;colorbar
title(["Gradient Power * Grid Weights",compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E_sp.tstat,P_D1E_sp,Dpr_D1E_sp)])
xlabel("PC3 Phi (deg)")
nexttile(T,4);
imagesc(xarr,yarr,abs(D2map_sp_btrp))
axis image;colorbar
title(["Laplacian Filtered Amplitude",compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E_sp.tstat,P_D2E_sp,Dpr_D2E_sp)])
xlabel("PC3 Phi (deg)")
nexttile(T,5);
imagesc(xarr,yarr,actmap_shfl)
axis image;colorbar
title("Neural Activation")
ylabel("Shuffled Control")
xlabel("PC3 Phi (deg)")
nexttile(T,6);
imagesc(xarr,yarr,GPmap_sp_shfl)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi (deg)")
nexttile(T,7);
imagesc(xarr,yarr,GPmap_sp_shfl.*detJac_map)
axis image;colorbar
title("Gradient Power * Grid Weights")
xlabel("PC3 Phi (deg)")
nexttile(T,8);
imagesc(xarr,yarr,abs(D2map_sp_shfl))
axis image;colorbar
title("Laplacian Filtered Amplitude")
xlabel("PC3 Phi (deg)")
if doplot
title(T,compose("Gradient map Compared (Sphere)\n%s Exp %d Space%d(%s) Unit %s",Animal,Expi,si,splab,unitlab))
saveallform(figdir, compose("Smooth_Map_sph_cmp_%s_Exp%d_%s_%s",Animal,Expi,unitlab,splab), h2)
end
end