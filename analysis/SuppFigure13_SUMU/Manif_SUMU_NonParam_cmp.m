%% Manifold paper Supplementary
% Compare non parametric statistics for SU and MU esp. activation heights, normVUS tuning width
Set_Path;
mat_dir = "O:\Mat_Statistics";
nonpardir = "O:\Manif_NonParam\summary";
%% Do population statistics of SU-MU 
NonParTab = readtable(fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"));
%%
SU_statvec = []; % the statistics, normAUS_bsl for SU 
MU_statvec = [];
SU_stattab = []; % other meta info for that chanel. 
MU_stattab = [];
for Animal = ["Alfa", "Beto"]
% load(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr');
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
for Expi = 1:max(NonParTab.Expi(NonParTab.Animal==Animal))
% section out the table for this animal + Exp, at PC23 space. 
exptab = NonParTab(NonParTab.Animal==Animal & NonParTab.Expi==Expi & NonParTab.space==1,:);
unit_num_arr = MapVarStats(Expi).units.unit_num_arr;
spikeID = MapVarStats(Expi).units.spikeID;
expstatvec = exptab.normAUS_bsl;
Fmsk = exptab.F_P<1E-2;
% R2msk = exptab.R2>0.5;
valmsk = spikeID>0;
assert(numel(expstatvec)==numel(spikeID))
pair_ids = get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&valmsk);
for j = 1:size(pair_ids,1)
SU_statvec(end+1,:) = expstatvec(pair_ids(j,1));
MU_statvec(end+1,:) = expstatvec(pair_ids(j,2));
SU_stattab = cat(1,SU_stattab,exptab(pair_ids(j,1),:));
MU_stattab = cat(1,MU_stattab,exptab(pair_ids(j,2),:));
end
end
end
%
ttest2_print(SU_statvec, MU_statvec, "SU normVUS", "MU normVUS", true)
%%
ttest2_print(SU_statvec, MU_statvec, "SU normVUS", "MU normVUS", true);
msk=SU_stattab.area=="V1"; ttest2_print(SU_statvec(msk), MU_statvec(msk), "V1 SU normVUS", "V1 MU normVUS", true);
msk=SU_stattab.area=="V4"; ttest2_print(SU_statvec(msk), MU_statvec(msk), "V4 SU normVUS", "V4 MU normVUS", true);
msk=SU_stattab.area=="IT"; ttest2_print(SU_statvec(msk), MU_statvec(msk), "IT SU normVUS", "IT MU normVUS", true);
msk=SU_stattab.prefchan==SU_stattab.chan; ttest2_print(SU_statvec(msk), MU_statvec(msk), "Driver SU normVUS", "Driver MU normVUS", true);

%% Plot the comparison of normVUS between SU MU
drv_msk = SU_stattab.prefchan==SU_stattab.chan;
V1msk = SU_stattab.area=="V1";
V4msk = SU_stattab.area=="V4";
ITmsk = SU_stattab.area=="IT";
Alfamsk = SU_stattab.Animal=="Alfa";
Betomsk = SU_stattab.Animal=="Beto";
%% Separated by Driver Channel or not
paired_stripe_plot({SU_statvec, MU_statvec}, ["SU", "MU"], ...
    {drv_msk,~drv_msk}, ["Driver","non-driver"], 'MarkerEdgeAlpha',0.65)
ylabel("normVUS bsl"); ylim([0,2*pi])
title("SU-MU Non-Parametric Tuning Width Comparison")
saveallform(nonpardir,"SUMU-normVUS_bsl_cmp_drv_sep")
%% Separated by Area 
paired_stripe_plot({SU_statvec, MU_statvec}, ["SU", "MU"], ...
    {V1msk,V4msk,ITmsk}, ["V1","V4","IT"], 'MarkerEdgeAlpha',0.65)
ylabel("normVUS bsl"); ylim([0,2*pi])
title("SU-MU Non-Parametric Tuning Width Comparison")
saveallform(nonpardir,"SUMU-normVUS_bsl_cmp_area_sep")
%% Separated by Animal 
paired_stripe_plot({SU_statvec, MU_statvec}, ["SU", "MU"], ...
    {Alfamsk,Betomsk}, ["Alfa","Beto"], 'MarkerEdgeAlpha',0.65)
ylabel("normVUS bsl"); ylim([0,2*pi])
title("SU-MU Non-Parametric Tuning Width Comparison")
saveallform(nonpardir,"SUMU-normVUS_bsl_cmp_anim_sep")

%% Tuning Map Range Comprison separated by Driver or not
paired_stripe_plot({SU_stattab.Act_range, MU_stattab.Act_range}, ["SU", "MU"], ...
    {drv_msk,~drv_msk}, ["Driver","non-driver"], 'MarkerEdgeAlpha',0.65)
ylabel("Activation Range"); %ylim([0,2*pi])
title("SU-MU Non-Parametric Amplitude Comparison")
saveallform(nonpardir,"SUMU-ActRange_cmp_drv_sep")
%% Tuning Map Range Comprison
paired_stripe_plot({SU_stattab.Act_range, MU_stattab.Act_range}, ["SU", "MU"], ...
    {V1msk,V4msk,ITmsk}, ["V1","V4","IT"], 'MarkerEdgeAlpha',0.65)
ylabel("Activation Range"); %ylim([0,2*pi])
title("SU-MU Non-Parametric Amplitude Comparison")
saveallform(nonpardir,"SUMU-ActRange_cmp_area_sep")

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