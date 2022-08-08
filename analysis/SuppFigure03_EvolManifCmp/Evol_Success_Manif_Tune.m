%% Evolution Successfulness and Manifold Tuning Width
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
Kent_path = "E:\OneDrive - Washington University in St. Louis\Manif_SUHash\summary";
load(fullfile(Kent_path, Animal+"_Manif_Kent_Fit.mat"),"Kent_stats")

SuccCorr_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Succ_Manif_Kappa";
%% Evolution Statistics Extraction
EvolSuccStats = repmat(struct(),1,length(EStats));
for Expi = 1:length(Stats)
blockn = EStats(Expi).evol.block_n;
block_mean_score = cellfun(@(psth)mean(psth(1,51:200,:),'all'),EStats(Expi).evol.psth(1:blockn-1));
[~,peakBlock]=max(block_mean_score);
if peakBlock == blockn-1, peakBlock = blockn-2; end % avoid the last generation.
endspsths = cell2mat(reshape(EStats(Expi).evol.psth(peakBlock:peakBlock+1),1,1,[]));
initpsths = cell2mat(reshape(EStats(Expi).evol.psth(2:3),1,1,[]));
window = [51:200];
initacts = squeeze(mean(initpsths(1,window,:),2));
endsacts = squeeze(mean(endspsths(1,window,:),2));
[H,P,CI,STATS] = ttest2(endsacts, initacts);
EvolSuccStats(Expi).tstat = STATS.tstat;
EvolSuccStats(Expi).DAOA = (mean(endsacts) - mean(initacts)) / mean(initacts);
end
DAOA_arr = arrayfun(@(E)E.DAOA, EvolSuccStats)';
tstat_arr = arrayfun(@(E)E.tstat, EvolSuccStats)';

%% Get the Kappa fit for the Manifold data.
drivechan_tune_kappa = zeros(length(Stats),1);
for Expi = 1:length(Stats)
% ui = find(Stats(Expi).units.pref_chan_id ==
% EStats(Expi).units.pref_chan_id); % this can fail since the numbering
% could change from Evol to Manif
ui = EStats(Expi).units.unit_num_arr(EStats(Expi).units.pref_chan_id);
si = 1; % number of subspace, in case there are multiple subspace there.
kappa = Kent_stats(Expi).act_fit{ui,si}.coef(4);
drivechan_tune_kappa(Expi) = kappa;
end
%%
corr(DAOA_arr,  drivechan_tune_kappa) 
corr(tstat_arr, drivechan_tune_kappa) 
%%
prefchan_arr = arrayfun(@(E) E.units.pref_chan,EStats);
prefname_arr = arrayfun(@(E) E.units.unit_name_arr(E.units.pref_chan_id),EStats);
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
ColorArr = [V1msk',V4msk',ITmsk']; % R for V1, G for V4, B for IT
% ColorArr = 'b'; 
%%
figure(1);
subtightplot(1,2,1,0.05,0.08,0.05)
S1 = scatter(tstat_arr, drivechan_tune_kappa, 25, ColorArr);
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
ylabel("Manif Tuning Kappa");xlabel("t stat for Evol")
title(compose("Corr Coef %.3f",corr(tstat_arr,drivechan_tune_kappa)))
subtightplot(1,2,2,0.05,0.10,0.05)
S2 = scatter(DAOA_arr,  drivechan_tune_kappa, 25, ColorArr);
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
msk = DAOA_arr < 2 & drivechan_tune_kappa <2 ;
xlabel("Delta Activation / Activation_0")
title(compose("Corr Coef %.3f (< 2 subregion Corr Coef %.3f)",...
    corr(DAOA_arr,drivechan_tune_kappa),corr(DAOA_arr(msk),drivechan_tune_kappa(msk))))
suptitle(compose("%s Correlation between Evol successfulness and Manif tuning width\n(R,G,B for V1 V4 IT)",Animal))
%%
saveas(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr.jpg"))
savefig(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr.fig"))
saveas(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr_area.jpg"))

%% Evolution Ref Signal to Noise Ratio
for Expi = 1:length(Stats)
refimgmean = cellfun(@(psth)mean(psth(1,window,:),'all'),EStats(Expi).ref.psth_arr);
refimgstd = cellfun(@(psth)std(mean(psth(1,window,:),2)),EStats(Expi).ref.psth_arr);
figure(3);
scatter(refimgmean,refimgstd)
ylabel("std");xlabel("mean")
title(compose("Exp %d pref chan %s", Expi, prefname_arr(Expi)));
pause
end


%% Newer version Plotting Script
SuccCorr_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Succ_Manif_Kappa";
Animal="Both"; Set_Path;
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
V1msk = (ExpTab_cmb.pref_chan <=48) & (ExpTab_cmb.pref_chan >= 33);
V4msk =  ExpTab_cmb.pref_chan > 48;
ITmsk =  ExpTab_cmb.pref_chan < 33;
Alfamsk = (ExpTab_cmb.Animal=="Alfa");
Betomsk = (ExpTab_cmb.Animal=="Beto");
%% Load the population stats
sumdir = "O:\Manif_Fitting\summary";
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];
% Create the masks
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
driver_fittab = poptab(drivermsk, :);
driverPC23_fittab = driver_fittab(driver_fittab.space==1,:);
%% 
figure;
scatter(ExpTab_cmb.DAOA_initmax(driverPC23_fittab.R2>0.5), driverPC23_fittab.kappa(driverPC23_fittab.R2>0.5))
%%
valmsk = driverPC23_fittab.R2>0.5;
tmptab = [ExpTab_cmb(:,"DAOA_initend"), driverPC23_fittab(:,"kappa")];
xscatter_plot(tmptab, "kappa", "DAOA_initend", {V1msk&valmsk, V4msk&valmsk, ITmsk&valmsk}, ["V1","V4","IT"],...
              "Evol Sucess ~ Manifold Sharpness", "Both", SuccCorr_dir);
%%
valmsk = driverPC23_fittab.R2>0.5;
tmptab = [ExpTab_cmb(:,"DAOA_initend"), driverPC23_fittab(:,"kappa")];
figh = xscatter_minor_plot(tmptab, "kappa", "DAOA_initend", {V1msk&valmsk, V4msk&valmsk, ITmsk&valmsk}, ["V1","V4","IT"],...
              {Alfamsk, Betomsk}, ["A","B"],"Evol Sucess ~ Manifold Sharpness", "Both_anim_area_sep", SuccCorr_dir);
%%
xlim([0,5]);ylim([0,8]);
saveallform(SuccCorr_dir, compose("%s_%s_xscat_%s", "kappa", "DAOA_initend", "Both_anim_area_sep_XYLIM"), figh)
%%
valmsk = driverPC23_fittab.R2>0.5;
diary(fullfile(SuccCorr_dir,"stat_summary.log"))
msks = {valmsk, valmsk&V1msk, valmsk&V4msk, valmsk&ITmsk, valmsk&Alfamsk, valmsk&Betomsk};
labels = ["All", "V1","V4","IT","Alfa","Beto"];
summary_corr_wmasks(ExpTab_cmb.DAOA_initend,"DAOA",driverPC23_fittab.kappa,"kappa",labels,msks)
summary_corr_wmasks(ExpTab_cmb.t_initend,"tstat",driverPC23_fittab.kappa,"kappa",labels,msks)
diary off
winopen(fullfile(SuccCorr_dir,"stat_summary.log"))
%%
summary_corr_wmasks(-ExpTab_cmb.t_initend,"tstat_initend",driverPC23_fittab.kappa,"kappa",labels,msks)%,"all",{[]})
function summary_corr_wmasks(vec1,name1,vec2,name2,labels,msks)
for i = 1:numel(msks)
msk = msks{i};
if isempty(msk), msk = ones(length(vec1),1,'logical');end
label = labels(i);
N = sum(msk);
[cval, pval] = corr(vec1(msk),vec2(msk),'type','spearman');
stat1_m = mean(vec1(msk));
stat1_s = sem(vec1(msk));
stat2_m = mean(vec2(msk));
stat2_s = sem(vec2(msk));
fprintf("For %s exps, Spearman Corr %.3f (%.1e) N=%d\n%s %.1f+-%.1f %s %.2f+-%.2f\n", label, cval, pval, N, name1, stat1_m, stat1_s, name2, stat2_m, stat2_s)
end
end