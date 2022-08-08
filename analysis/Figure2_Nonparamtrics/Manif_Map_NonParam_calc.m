%% Compute the nonparametric statistics for tuning maps of all exp all channels pop. 
%  Similar content to Manif_NonParametrics_Pop.m 
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];

%% Compute and curate non parametric stats about single channel tuning map
tic
S_col = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
for Expi = 1:numel(Stats)
for spi = 1:MapVarStats(Expi).manif.subsp_n
exptab = poptab(poptab.Animal==Animal & poptab.Expi==Expi & poptab.space==spi,:);
assert(size(exptab,1)==numel(MapVarStats(Expi).units.spikeID))
for iCh = 1:numel(MapVarStats(Expi).units.spikeID)
S = struct();
S.Animal = Animal;
S.Expi = Expi;
S.chan = MapVarStats(Expi).units.spikeID(iCh);
S.unitnum = MapVarStats(Expi).units.unit_num_arr(iCh);
S.prefchan = MapVarStats(Expi).units.pref_chan;
S.unitstr = MapVarStats(Expi).units.unit_name_arr(iCh);
S.area = area_map(S.chan);
S.iCh = iCh;
S.F = exptab.F(iCh);
S.F_P = exptab.F_P(iCh);
S.space = spi;
% for nm = ["Animal", "Expi", "unitstr", "unitnum", "chan", "prefchan", "space","F", "F_P"] % copy stats.
% S.(nm) = Tab.(nm);
% end
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2), MapVarStats(Expi).manif.act_col{spi}, 'uni',1);
bslmat = cell2mat(reshape(MapVarStats(Expi).manif.bsl_col{spi},1,[]));
bslmean = mean(bslmat(iCh,:));
% bslstd = std(bslmat(iCh,:));
% Major workhorse here 
S = Manif_Map_NonParam_stat_fun(actmap_mean, bslmean, S);
S_col = [S_col;S];
end
end
toc
end
end % 74 sec for computation
%%
nonpardir = "O:\Manif_NonParam\summary";
NonParTab = struct2table(S_col);
writetable(NonParTab,fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"))
% Check `fullfile(nonpardir,"Both"+"_Pop_NonParamStat.csv"))`


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

function S = Manif_Map_NonParam_stat_fun(actmap_mean, bslmean, S)
if nargin==2, S = struct(); end
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
XX = cosd(theta_grid).* cosd(phi_grid);
YY = sind(theta_grid) .* cosd(phi_grid);
ZZ = sind(phi_grid);
%  Integration Weight Matrix: 
%  Integrate the area of the patch that each node governs
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9, -90) /180 *pi;
theta2_grid = min(theta_grid + 9,  90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);

[maxAct, maxIdx] = max(actmap_mean,[],'all','linear');
[minAct, ~] = min(actmap_mean,[],'all','linear');
S.maxAct = maxAct;
S.minAct = minAct;
S.Act_range = maxAct - minAct;
S.bslAct = bslmean;
thresh = 0.9 * maxAct; 
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_90 = CoMtheta;
S.CoMphi_90 = CoMphi;
S.CoMrho_90 = CoMrho;
S.Mtheta_90 = Mtheta;
S.Mphi_90 = Mphi;
S.Mrho_90 = Mrho;
% Core of the loop, use the Spherical Integration to get the activation maps. 
integ = SphIntegration(actmap_mean,0,theta_grid,phi_grid,Wgrid);
S.AUS = integ; % Area under tuning map
S.normAUS = integ / maxAct; % Normalized area under tuning map.
integ = SphIntegration(actmap_mean,bslmean,theta_grid,phi_grid,Wgrid);
S.AUS_bsl = integ;
S.normAUS_bsl = integ / (maxAct - bslmean);
end

function [CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap, thresh, XX, YY, ZZ, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
msk = actmap > thresh;
Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
Actvec = actmap(msk);
CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
% [TH,PHI,R] = cart2sph(CoMvec(),CoMvec(),CoMvec())
[CoMtheta, CoMphi, CoMrho] = cart2sph(CoMvec(1),CoMvec(2),CoMvec(3));
[Mtheta, Mphi, Mrho] = cart2sph(Mvec(1),Mvec(2),Mvec(3));
end

function [Integ,meanInteg] = SphIntegration(actmap, baseline, theta_grid, phi_grid, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
mask = (actmap - baseline) > 0;
Integ = sum(mask .* (actmap - baseline) .* Weight,'all');
meanInteg = sum(mask .* (actmap - baseline) .* Weight,'all') / sum(Weight,'all');
% Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
% Actvec = actmap(msk);
% CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
% Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
end