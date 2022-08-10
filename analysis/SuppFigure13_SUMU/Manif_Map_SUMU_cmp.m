%% Manifold Paper Supplementary Figure. Single unit, Multi unit comparison. Compared across areas
%  Compare all aspects of SU, MU. Kent fitting, Non parametrics, Corr, 
%  Code to compare paired Single and Multi unit in the dataset.

sumdir = "O:\CortiDistCorr\summary";
SUMUdir = "O:\Manif_SUHash\summary";
nonpardir = "O:\Manif_NonParam\summary";
popfitdir = "O:\Manif_Fitting\popstats";
Animal = "Both";Set_Path;
% NPStatTab = readtable(fullfile(nonpardir,Animal+"_Driver_NonParamStat.csv"));
% NPStatTab = readtable(fullfile(nonpardir,Animal+"_Driver_NonParamWidth.csv"));
nonpardir = "O:\Manif_NonParam\summary";
NonParTab = readtable(fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"),'Format','auto');
%% Load the Kent fitting data for the Population
popfittab = readtable(fullfile(popfitdir,"Both_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
%% Load the tuning map similarity matrix 
Animal = "Both";Set_Path;
mat_dir = "O:\Mat_Statistics";
if strcmp(Animal, "Both") % load stats
A = load(fullfile(mat_dir, "Alfa"+'_CortiDisCorr.mat'), 'CortiDisCorr');
B = load(fullfile(mat_dir, "Beto"+'_CortiDisCorr.mat'), 'CortiDisCorr');
CortiDisCorr = [A.CortiDisCorr, B.CortiDisCorr];
clear A B
else
load(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr') 
end

%%  Supplementary Figure 13
%   Find examples of SU-MU Tune map pair and show examples montage of their maps
Animal="Alfa"; iTr = 15;
unit_num_arr = CortiDisCorr(iTr).units.unit_num_arr;
spikeID = CortiDisCorr(iTr).units.spikeID;
Fmsk = struct2table(CortiDisCorr(iTr).FStats).F_P<1E-2;
valmsk = unit_num_arr>0;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
ITmsk =  spikeID < 33;
pairs = get_SUMUpair_ids(unit_num_arr', spikeID, valmsk&Fmsk);%&ITmsk
for rowi = 1:size(pairs,1)
Manif_Map_show_fun(MapVarStats, "Alfa", iTr, pairs(rowi,:))
pause
end
%%
saveallform(SUMUdir, compose("SUMUpair_cmp_%s_Exp%02d_%s",Animal,iTr,strrep(num2str(pairs(rowi,:)),"  ","_")))

%% Compare the max activation for SU and MU
SU_statvec = [];
MU_statvec = [];
% for Animal = ["Alfa", "Beto"]
% load(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr');
% for Expi = 1:max(popfittab.Expi(popfittab.Animal==Animal))
for iTr = 1:numel(CortiDisCorr)
Animal = CortiDisCorr(iTr).Animal;
Expi = CortiDisCorr(iTr).Expi;
unit_num_arr = CortiDisCorr(iTr).units.unit_num_arr;
spikeID = CortiDisCorr(iTr).units.spikeID;
% Get the stats table for this experiment
exptab = popfittab(popfittab.Animal==Animal & popfittab.Expi==Expi & popfittab.space==1,:);
% expstatvec = exptab.kappa;
NPexptab = NonParTab(NonParTab.Animal==Animal & NonParTab.Expi==Expi & NonParTab.space==1,:);
expstatvec = NPexptab.maxAct;
Fmsk = struct2table(CortiDisCorr(iTr).FStats).F_P<1E-2; % Well modulated
R2msk = true;%exptab.R2>0.5; % Well fit 
valmsk = unit_num_arr>0; % unit 1,2,3 not 0
assert(numel(expstatvec)==numel(spikeID))
pair_ids = get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&R2msk&valmsk);
for j = 1:size(pair_ids,1)
SU_statvec(end+1,:) = expstatvec(pair_ids(j,1));
MU_statvec(end+1,:) = expstatvec(pair_ids(j,2));
end
end
ttest2_print(SU_statvec, MU_statvec, "SU maxAct", "MU maxAct", true)

%% SU,MU kappa value (parametric tuning width) of all pairs
SU_statvec = [];
MU_statvec = [];
for iTr = 1:numel(CortiDisCorr)
Animal = CortiDisCorr(iTr).Animal;
Expi = CortiDisCorr(iTr).Expi;
unit_num_arr = CortiDisCorr(iTr).units.unit_num_arr;
spikeID = CortiDisCorr(iTr).units.spikeID;
% Get the stats table for this experiment
exptab = popfittab(popfittab.Animal==Animal & popfittab.Expi==Expi & popfittab.space==1,:);
expstatvec = exptab.kappa;
Fmsk = struct2table(CortiDisCorr(iTr).FStats).F_P<1E-2; % Well modulated
R2msk = exptab.R2>0.5; % Well fit 
valmsk = unit_num_arr>0; % unit 1,2,3 not 0
assert(numel(expstatvec)==numel(spikeID))
pair_ids = get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&R2msk&valmsk);
for j = 1:size(pair_ids,1)
SU_statvec(end+1,:) = expstatvec(pair_ids(j,1));
MU_statvec(end+1,:) = expstatvec(pair_ids(j,2));
end
end
ttest2_print(SU_statvec, MU_statvec, "SU kappa", "MU kappa", true)

%% Compare the similarity of SU-MU
%% Collect pairs of SU MU correlation into 
areacorrvec_col = {}; % nExp by 5 cell array, containing vectorized corr coef
SUMUcorrvec_col = {}; % nExp by 5 cell array, containing vectorized corr coef for SU-MU
for iTr=1:numel(CortiDisCorr)
% Fmsk = struct2table(CortiDisCorr(1).FStats).F_P<1E-2;
unit_num_arr = CortiDisCorr(iTr).units.unit_num_arr;
spikeID = CortiDisCorr(iTr).units.spikeID;
Fmsk = struct2table(CortiDisCorr(iTr).FStats).F_P<1E-2;
valmsk = unit_num_arr>0;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
ITmsk =  spikeID < 33;
corrmat = CortiDisCorr(iTr).avgsph_corrmat;
% corrmat = CortiDisCorr(1).res_corrmat;
msk_col = {V1msk, V4msk, ITmsk, true};
label_col = ["V1","V4","IT","all"];
for mi = 1:numel(msk_col)
areamsk = msk_col{mi}; arealab = label_col{mi};
pair_ids = get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&areamsk);
SUMUpair_linidx = sub2ind( size(corrmat), pair_ids(:,1), pair_ids(:,2));
SUMUcorrvec = corrmat(SUMUpair_linidx);
in_area_pair_corrvec = get_submat_value(corrmat, areamsk & valmsk & Fmsk, SUMUpair_linidx);
ttest2corr_print(SUMUcorrvec, in_area_pair_corrvec, compose("%s SU-MU",label_col(mi)), compose("%s val pairs",label_col(mi)));
SUMUcorrvec_col{iTr,mi} = SUMUcorrvec;
areacorrvec_col{iTr,mi} = in_area_pair_corrvec;
end

pair_ids = get_SUMUpair_ids(unit_num_arr',spikeID, Fmsk);
SUMUpair_linidx = sub2ind( size(corrmat), pair_ids(:,1), pair_ids(:,2));
SUMUcorrvec = corrmat(SUMUpair_linidx);
all_area_pair_corrvec = cat(1, get_submat_value(corrmat, V1msk & valmsk, SUMUpair_linidx)...
                             , get_submat_value(corrmat, V4msk & valmsk, SUMUpair_linidx)...
                             , get_submat_value(corrmat, ITmsk & valmsk, SUMUpair_linidx));
ttest2corr_print(SUMUcorrvec, all_area_pair_corrvec, "All SU-MU", "All same area pairs");
SUMUcorrvec_col{iTr,5} = SUMUcorrvec;
areacorrvec_col{iTr,5} = all_area_pair_corrvec;
end
%
ttest2corr_print(cat(1,SUMUcorrvec_col{:,1}), cat(1,areacorrvec_col{:,1}), "All Exp SU-MU (in V1)", "All Exp same area pairs (in V1)");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,2}), cat(1,areacorrvec_col{:,2}), "All Exp SU-MU (in V4)", "All Exp same area pairs (in V4)");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,3}), cat(1,areacorrvec_col{:,3}), "All Exp SU-MU (in IT)", "All Exp same area pairs (in IT)");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,5}), cat(1,areacorrvec_col{:,5}), "All Exp SU-MU", "All Exp same area pairs");
%
label_col = ["V1","V4","IT","all","all_same_array"];
save(fullfile(sumdir,Animal+"SUMU_MapCorr_cmp.mat"), "SUMUcorrvec_col", "areacorrvec_col", "label_col");
save(fullfile(matdir,Animal+"SUMU_MapCorr_cmp.mat"), "SUMUcorrvec_col", "areacorrvec_col", "label_col");


%% Visualize and Summary
%cat(1,SUMUcorrvec_col{:,1}), cat(1,areacorrvec_col{:,1})
h=stripe_simple_plot({cat(1,SUMUcorrvec_col{:,5}), cat(1,areacorrvec_col{:,5})},...
        "corr",["SUMU","All other"], sumdir, "SUMU_bsl_cmp", "markeredgealpha", 0.2);
%%
h=figure();set(h,'pos',[1000  412   400   600]);
h = violin({cat(1,SUMUcorrvec_col{1:46,5}), ...
            cat(1,areacorrvec_col{1:46,5}), ...
            cat(1,SUMUcorrvec_col{47:end,5}), ...
            cat(1,areacorrvec_col{47:end,5})},'facealpha',0.4);
xticks(1:4); box off
xticklabels(["SU-MU A","All A","SU-MU B","All B"]);
ylabel("Correlation")
title(compose("Compare Tuning Map Similarity\n in SU-MU vs Other Pairs"))
saveallform(sumdir,"SUMU_other_cmp_violin",h)

%%
h=figure;set(h,'pos',[1000  412   484   555]);
cc_col = {cat(1,SUMUcorrvec_col{1:46,5}),...
          cat(1,areacorrvec_col{1:46,5}),...
          cat(1,SUMUcorrvec_col{47:end,5}),...
          cat(1,areacorrvec_col{47:end,5})};
label_arr = ["SU-MU A", "All A", "SU-MU B", "All B"];
violinplot_cell(cc_col, label_arr,'showData',true)
ylabel("Correlation")
saveallform(sumdir,"SUMU_other_cmp_violin_scatter",h)


function corrvec = get_submat_value(corrmat, msk, exclude_linidx)
linidxmat = reshape(1:numel(corrmat),size(corrmat));
linidx_submat = linidxmat(msk,msk);
linidx_submat_tril = tril(linidx_submat,-1);
linidx_final = nonzeros(linidx_submat_tril);
if nargin > 2
    linidx_final = setdiff(linidx_final, exclude_linidx);
end
corrvec = corrmat(linidx_final);
end

function pairs = get_SUMUpair_ids(unit_num_arr,chan_num_arr,msk)
U1idxs = strfind(unit_num_arr,[1,2]);
pairs = [U1idxs',U1idxs'+1];
chan_ids = chan_num_arr(pairs);
assert(all(chan_num_arr(U1idxs')==chan_num_arr(U1idxs'+1)),"the spike id doesn't match ")
if nargin >2
validunit = msk(pairs);
if size(pairs,1)==1, validunit=validunit'; end
valid_row = find(all(validunit,2));
pairs = pairs(valid_row,:);
assert(all(validunit(valid_row,:),'all'),"the spike id doesn't match ")
end
end
