%% Manifold Paper Figure 4: Tuning Width Progression along the Ventral Stream
%  Post processing of the manifold fitting data. 
%  Majorly test the trend of Kappa and Beta. 
%  Write the diary to the summary folder 
poptabdir = "O:\Manif_Fitting\popstats";
sumdir = "O:\Manif_Fitting\Kent_summary";

poptab = readtable(fullfile(poptabdir,"Both_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
%% Create masks for the population. 
valmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10)); % excluding one exp with fewer trials
Fmsk = poptab.F_P<0.001;
V1msk = poptab.area == "V1";
V4msk = poptab.area == "V4";
ITmsk = poptab.area == "IT";
Alfamsk = poptab.Animal == "Alfa";
Betomsk = poptab.Animal == "Beto";
drivermsk = poptab.is_driver;

%% create array of masks and public mask 
Vmsk = drivermsk&poptab.R2>.5; % is driver and the fitting is good. 
msks = {Vmsk&V1msk,Vmsk&V4msk,Vmsk&ITmsk}; 
labels = ["V1","V4","IT"];

%%
varnm = "kappa";
diary(fullfile(sumdir,"sumstat.log"))
summary_stat(poptab, varnm, {Vmsk&V1msk,Vmsk&V4msk,Vmsk&ITmsk}, ["V1","V4","IT"])
summary_stat(poptab, varnm, {Vmsk&Alfamsk,Vmsk&Betomsk}, ["Alfa","Beto"])
fprintf("\nProgression of %s value in Both monkey\n",varnm)
lm = test_progression(poptab, varnm, {Vmsk&V1msk,Vmsk&V4msk,Vmsk&ITmsk}, ["V1","V4","IT"]);
Vmsk = drivermsk&poptab.R2>.5;%&(poptab.space==1);
fprintf("\nProgression of %s value in Alfa\n",varnm)
lm_A = test_progression(poptab, varnm, {Vmsk&V1msk&Alfamsk,Vmsk&V4msk&Alfamsk,Vmsk&ITmsk&Alfamsk}, ["V1","V4","IT"]);
fprintf("\nProgression of %s value in Beto\n",varnm)
lm_B = test_progression(poptab, varnm, {Vmsk&V1msk&Betomsk,Vmsk&V4msk&Betomsk,Vmsk&ITmsk&Betomsk}, ["V1","V4","IT"]);
diary off
%
varnm = "beta";
diary(fullfile(sumdir,"sumstat.log"))
summary_stat(poptab, varnm, {Vmsk&V1msk,Vmsk&V4msk,Vmsk&ITmsk}, ["V1","V4","IT"])
summary_stat(poptab, varnm, {Vmsk&Alfamsk,Vmsk&Betomsk}, ["Alfa","Beto"])
fprintf("\nProgression of %s value in Both monkey\n",varnm)
lm = test_progression(poptab, varnm, {Vmsk&V1msk,Vmsk&V4msk,Vmsk&ITmsk}, ["V1","V4","IT"]);
Vmsk = drivermsk&poptab.R2>.5;%&(poptab.space==1);
fprintf("\nProgression of %s value in Alfa\n",varnm)
lm_A = test_progression(poptab, varnm, {Vmsk&V1msk&Alfamsk,Vmsk&V4msk&Alfamsk,Vmsk&ITmsk&Alfamsk}, ["V1","V4","IT"]);
fprintf("\nProgression of %s value in Beto\n",varnm)
lm_B = test_progression(poptab, varnm, {Vmsk&V1msk&Betomsk,Vmsk&V4msk&Betomsk,Vmsk&ITmsk&Betomsk}, ["V1","V4","IT"]);
diary off

function summary_stat(tab,varnm,msks,labels)
% Akin to pandas `groupby` `agg` syntax. 
% Summarize scalar statistics under the masks with labels.
var_cols = {};
fprintf("\n%s value\n",varnm)
for i = 1:numel(msks)
msk = msks{i};
var_arr = tab.(varnm)(msk);
fprintf("%s: mean %.3f +- std %.3f sem %.3f med %.3f [5-95]Prct [%.3f-%.3f] (N=%d)\n",...
    labels(i),mean(var_arr),std(var_arr),sem(var_arr),median(var_arr),prctile(var_arr,5),prctile(var_arr,95),numel(var_arr))
var_cols{end+1} = var_arr;
end
S = anova_cells(var_cols);
fprintf("ANOVA F=%.3f p=%.1e (df=%d)\n",S.F,S.F_P,S.STATS.df)
end

function lm = test_progression(tab,varnm,msks,labels)
% Test the progression 
%     via linear model fitting for the masks and labels. 
var_vec = []; idx_vec = [];
for i = 1:numel(msks)
msk = msks{i};
var_arr = tab.(varnm)(msk);
idx_arr = i*ones(numel(var_arr),1);
var_vec = [var_vec; var_arr];
idx_vec = [idx_vec; idx_arr];
end
lm = fitlm(idx_vec, var_vec);
disp(lm)
end
