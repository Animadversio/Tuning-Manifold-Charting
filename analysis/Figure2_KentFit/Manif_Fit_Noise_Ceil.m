%% Manif Fit Noise Ceiling 

Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_Fitting\popstats";
% load(fullfile(Matdir, "Beto_ManifPopDynamics.mat"))
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')

%% 
expvar_arr_all = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
expvar_arr = nan(numel(Stats,1),1);
for Expi = 1:numel(Stats)
    if Expi == 10 && Animal=="Alfa", continue;end
    actcolmap = cellfun(@(P)squeeze(mean(P(1,50:200,:),[1,2])), Stats(Expi).manif.psth{1}, 'unif', 0);
    expvar_arr(Expi) = bstrp_expvar(actcolmap);
end
nanmean(expvar_arr)
expvar_arr_all = [expvar_arr_all; expvar_arr];
end

%% 
Scol = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
for Expi = 1:numel(Stats)
    S = struct();
    S.Animal = Animal;
    S.Expi = Expi;
    S.prefchan = Stats(Expi).units.pref_chan;
    S.area = area_map(S.prefchan);
    if Expi == 10 && Animal=="Alfa"
        S.space = 1;
        S.expvar = nan;
        Scol = [Scol;S]; 
        continue; 
    end
    for spi = 1:numel(Stats(Expi).manif.psth)
    actcolmap = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])), Stats(Expi).manif.psth{spi}, 'unif', 0);
    S.space = spi;
    S.expvar = bstrp_expvar_fixtarg(actcolmap);
    Scol = [Scol;S];
    csr = csr+1;
    end
end
end
tab = struct2table(Scol);
%%
writetable(tab, fullfile(tabdir,"Both_Driver_ExpVar.csv"))
%%
fprintf("Bootstrapped Exp Var. mean %.3f (N=%d) median %.3f\n",mean(expvar_arr_all),numel(expvar_arr_all),median(expvar_arr_all))
msk_col = {tab.area=="V1", tab.area=="V4", tab.area=="IT", tab.Animal=="Alfa", tab.Animal=="Beto"};
lab_col = ["V1","V4","IT","Alfa","Beto"];
summarize_masks_print(tab, "expvar", msk_col, lab_col, {[]},[""]);
for mi = 1:numel(msk_col)
    fprintf("ExpVar in %s %.3f (N=%d)\n", lab_col(mi), nanmean(tab.expvar(msk_col{mi})),sum(~isnan(tab.expvar(msk_col{mi}))))
end
%%
summarize_masks_print(tab, "expvar", msk_col, lab_col, {drivertab.F_P<0.001},["Fsig"]);
%% 
drivertab = readtable(fullfile(tabdir,"Both_Exp_Driver_KentStat_bsl_pole.csv"),'Format','auto');
%% Load up Kent Stats and get the table containing the drivers. 
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
alfatab_pop = readtable(fullfile(tabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
betotab_pop = readtable(fullfile(tabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
poptab = [alfatab_pop;betotab_pop];
poptab.area = arrayfun(@area_map, poptab.chan);
for i = 1:size(poptab,1)
    poptab.imgsize(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgsize;
    poptab.imgposX(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(1);
    poptab.imgposY(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(2);
end
% %% Creat masks for analysis
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = poptab.area == "V1";%(poptab.chan<=48 & poptab.chan>=33);
V4msk = poptab.area == "V4";%(poptab.chan>48);
ITmsk = poptab.area == "IT";%(poptab.chan<33);
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
prefchmsk = poptab.chan==poptab.prefchan;
Fsigmsk = poptab.F_P<0.001;
R2msk = poptab.R2>0.5;
drivertab = poptab(drivermsk,:);
writetable(drivertab,fullfile(tabdir,"Both_Exp_Driver_KentStat_bsl_pole.csv"))
%% Get the normalized explained variance
tab.normR2 = drivertab.R2 ./ tab.expvar;
summarize_masks_print(tab, "normR2", msk_col, lab_col, {[]}, ["all"])
summarize_masks_print(tab, "normR2", msk_col, lab_col, {drivertab.F_P<0.001}, ["Fsig"])


%% Test expvar 
expvar_both = bstrp_expvar(actcolmap);
expvar_fixtrg = bstrp_expvar_fixtarg(actcolmap,500);
fprintf("Usually with fixed target, exp var is larger. %.3f < %.3f\n",expvar_both,expvar_fixtrg)
%  Testing expvar
test_expvarfunc()
function test_expvarfunc(actcolmap)
if nargin == 0
    actcolmap = reshape(arrayfun(@(i)rand(9,1),1:121,'unif',0),[11,11]);
end
rng(0);expvar1 = bstrp_expvar(actcolmap);
rng(0);expvar2 = bstrp_expvar_forloop(actcolmap);
assert(abs(expvar1-expvar2)<eps, "Implementation error in explained variance. ")
fprintf("Implementation match!\n")
end

function mean_expvar = bstrp_expvar_fixtarg(actcolmap, Nbtrp)
if nargin == 1, Nbtrp=500; end
actmap_mean = cellfun(@mean, actcolmap);
actcolmap_btrp = cellfun(@(P)reshape(bootstrp(Nbtrp, @mean, P),1,1,[]), actcolmap, 'unif', 0);
acttsr_btrp = cell2mat(actcolmap_btrp);
actmat_btrp = reshape(acttsr_btrp, 11*11, Nbtrp);
actvec = reshape(actmap_mean, 11*11, 1);
L2mat = pdist2(actmat_btrp',actvec','euclidean'); % row needs to be a actmap pattern. 
resvarmat = L2mat .^2 / 121;
varval = var(actmap_mean,1,'all');%var(actmat_btrp,1,1);
expvarmat = 1 - resvarmat ./ varval;
mean_expvar = nanmean(expvarmat,'all');
% disp(mean_expvar)
end

function mean_expvar = bstrp_expvar(actcolmap, Nbtrp)
if nargin == 1, Nbtrp=100; end
actcolmap_btrp = cellfun(@(P)reshape(bootstrp(Nbtrp, @mean, P),1,1,[]), actcolmap, 'unif', 0);
acttsr_btrp = cell2mat(actcolmap_btrp);
actmat_btrp = reshape(acttsr_btrp,11*11, Nbtrp);
L2mat = squareform(pdist(actmat_btrp','euclidean')); % row needs to be a actmap pattern. 
resvarmat = L2mat .^2 / 121;
varvec = var(actmat_btrp,1,1);
expvarmat = 1 - resvarmat ./ varvec;
expvarmat = expvarmat + diag(nan(1,100));
mean_expvar = nanmean(expvarmat,'all');
% disp(mean_expvar)
end

function mean_expvar = bstrp_expvar_forloop(actcolmap, Nbtrp)
	if nargin == 1, Nbtrp=100; end
actcolmap_btrp = cellfun(@(P)reshape(bootstrp(Nbtrp, @mean, P),1,1,[]), actcolmap, 'unif', 0);
acttsr_btrp = cell2mat(actcolmap_btrp);
actmat_btrp = reshape(acttsr_btrp,11*11, Nbtrp);
for i = 1:100
	for j = 1:100 
	expvarmat(i,j) = 1 - sum((actmat_btrp(:,i)-actmat_btrp(:,j)).^2) / sum((actmat_btrp(:,j) - mean(actmat_btrp(:,j))).^2);
	end
end 
expvarmat = expvarmat + diag(nan(1,100));
mean_expvar = nanmean(expvarmat,'all');
% disp(mean_expvar)
end

% %% 
% actcolmap = cellfun(@(P)squeeze(mean(P(1,50:200,:),[1,2])), Stats(1).manif.psth{1}, 'unif', 0);
% %%
% actcolmap_btrp = cellfun(@(P)reshape(bootstrp(100, @mean, P),1,1,[]), actcolmap, 'unif', 0);
% acttsr_btrp = cell2mat(actcolmap_btrp);
% actmat_btrp = reshape(acttsr_btrp,11*11,100);
% %% 1 - corr^2
% corrmat_actbtrp = corr(actmat_btrp);
% corrmat_actbtrp = corrmat_actbtrp + diag(nan(1,100));
% mean_corrsq = nanmean(corrmat_actbtrp.^2,'all');
% disp(mean_corrsq)
% %% 1 - res / var
% L2mat = squareform(pdist(actmat_btrp','euclidean')); % row needs to be a actmap pattern. 
% resvarmat = L2mat .^2 / 121;
% varvec = var(actmat_btrp,1,1);
% expvarmat = 1 - resvarmat ./ varvec;
% expvarmat = expvarmat + diag(nan(1,100));
% mean_expvar = nanmean(expvarmat,'all');
% disp(mean_expvar)
% %% 