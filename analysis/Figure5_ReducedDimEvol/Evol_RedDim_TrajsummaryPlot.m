%% Create the summary plots for Compare evolution scores for Alfa/Beto.
%  See if the Evolution works 
%  (Newer Refactored version using new plotting API comparing to Evol_RedDim_summary.m. May. 2021)
%  Plot the evolution trajectory comparison
%  With handy stats function at bottom. 

%% Set path and load most basic exp summary
Animal = "Both"; Set_Path;
global figdir
figdir = "O:\Evol_ReducDim\summary";

ExpType = "RDEvol";
if strcmp(Animal,"Both") % load stats
A = load(fullfile(matdir, "Alfa"+"_RDEvol_stats.mat"), 'RDStats');
B = load(fullfile(matdir, "Beto"+"_RDEvol_stats.mat"), 'RDStats');
RDStats = [A.RDStats, B.RDStats];
else
load(fullfile(matdir, Animal+"_RDEvol_stats.mat"), 'RDStats')
end
Corder = colororder; % default color seq

%% Collect all traces into a cell array from the `RDStats`.
%  For future plotting purpose.
RDEvol_Stats = []; % struct list, higher level summary, easy for plotting and testing.
unitnum_arr = zeros(length(RDStats),1);
prefchan_arr = zeros(length(RDStats),1);
validmsk = ones(numel(RDStats), 1, 'logical'); % whether the exp should be excluded. (unconventional optimizer setup)
score_traces = cell(numel(RDStats), 2); % activation vector 
zscore_traces = cell(numel(RDStats), 2); % activation vector zscored for each session 
block_traces = cell(numel(RDStats), 2); % vector of block num for plotting
% ref_traces = cell(numel(RDStats), 2);
% ref_block_traces = cell(numel(RDStats), 2); % vector of block num for plotting
for Expi = 1:length(RDStats)
thread_n = RDStats(Expi).evol.thread_num;
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == RDStats(Expi).evol.unit_in_pref_chan(1)) ...
          & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
ui = find(prefchan_id==chid); % unit index in the evol.psth array
unitnum_arr(Expi) = RDStats(Expi).evol.unit_in_pref_chan(1);
prefchan_arr(Expi) = RDStats(Expi).evol.pref_chan(1);
% Collect scaler score of each trial and the block num into cell array
window = [51:200]; bslwdw = [1:50];
evoke_vec_col = cellfun(@(psth)squeeze(mean(psth(ui,window,:),[1,2])),RDStats(Expi).evol.psth,'uni',0);
bsl_vec_col =  cellfun(@(psth)squeeze(mean(psth(ui,bslwdw,:),[1,2])),RDStats(Expi).evol.psth,'uni',0);
bslvec = cat(1, bsl_vec_col{:}); % Debug and update the baseline calculation method @July 4th. 
bslmean = mean(bslvec); 
score_vec_col = cellfun(@(act) act-bslmean, evoke_vec_col, 'uni',0);
zscore_vec_col = zscore_cellarr(score_vec_col);
block_vec_col = cellfun(@(idx)RDStats(Expi).evol.block_arr(idx),...
                              RDStats(Expi).evol.idx_seq,'uni',0);

S = struct();
S.Animal = RDStats(Expi).Animal;
S.Expi = RDStats(Expi).Expi;
S.pref_chan = RDStats(Expi).evol.pref_chan(1);
S.pref_unit = RDStats(Expi).evol.unit_in_pref_chan(1);
% Generate stats on the cell array of scores. 
if all(RDStats(Expi).evol.optim_names == ["ZOHA Sphere lr euclid", "ZOHA Sphere lr euclid ReducDim"]) 
% make sure the thread order optimizer order is correct
for thr_i = [1,2] % collect the vectorized scores into the summary stats
score_traces{Expi, thr_i} = cat(1,score_vec_col{thr_i,1:end-1}); % Full vector of score
zscore_traces{Expi, thr_i} = cat(1,zscore_vec_col{thr_i,1:end-1}); % Full vector of zscore
block_traces{Expi, thr_i} = cat(1,block_vec_col{thr_i,1:end-1});
% Collect Stats
end
midgen = round((size(score_vec_col,2) - 1)/2);
S = optimtraj_integral(score_vec_col, S);
S = Dprime_integral(score_vec_col,S);
S = score_cmp_stats(score_vec_col(:,2:3), "init23", S);
S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23", S);
S = score_cmp_stats(score_vec_col(:,midgen:midgen+1), "middle", S);
% quantify the successfulness of evolution. 
[t_full,p_full,sumstr_full]=ttest2_print(cat(1,score_vec_col{1,1:2}),cat(1,score_vec_col{1,end-1:end}),"init12","end12");
[t_50,p_50,sumstr_50]=ttest2_print(cat(1,score_vec_col{2,1:2}),cat(1,score_vec_col{2,end-1:end}),"init12","end12");
S.t_succ = [t_full,t_50];
S.t_p_succ = [p_full,p_50];
S.sumstr_succ = [sumstr_full,sumstr_50];
else
fprintf("Exp %02d, Skip Non standard optimnames",Expi)
disp(RDStats(Expi).evol.optim_names)
validmsk(Expi) = 0; % Mask out this exp when computing stats. 
S = optimtraj_integral({}, S);
S = Dprime_integral({}, S);
S = score_cmp_stats({}, "init23", S);
S = score_cmp_stats({}, "last23", S);
S = score_cmp_stats({}, "middle", S);
S.t_succ = [nan,nan];
S.t_p_succ = [nan,nan];
S.sumstr_succ = ["",""];
end
RDEvol_Stats = [RDEvol_Stats, S];
end
RDEvolTab = struct2table(RDEvol_Stats);
%% save the collected stats as mat and csv
writetable(RDEvolTab, fullfile(figdir, Animal+"_RDEvol_trajcmp.csv"))
save(fullfile(figdir, Animal+"_RDEvol_summaryStats.mat"), "RDEvol_Stats")

%% Newer API verions: Process the collected data from all the experiment pairs. 
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["Full", "50D"], "max12", 55);

%% Plot the summary of trajectories.
figdir = "O:\Evol_ReducDim\summary";
outdir = "O:\Manuscript_Manifold\Figure3\RedDimEffectProg";
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),RDStats');
area_arr = arrayfun(@area_map,prefchan_arr);
anim_arr = [RDStats.Animal]';
%% creat masks to filter array
V1msk = area_arr=="V1";
V4msk = area_arr=="V4";
ITmsk = area_arr=="IT";
Alfamsk = anim_arr=="Alfa";
Betomsk = anim_arr=="Beto";
anysucsmsk = any(RDEvolTab.t_p_succ<0.01,2);
allsucsmsk = all(RDEvolTab.t_p_succ<0.01,2);

%% Plot optimization trajectory separated by visual area, with individual optim traces 
figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
	{V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {}, [], ["Full", "50D"], true);

%% Plot optimization trajectory separated by visual area and animal, with individual optim traces 
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {Alfamsk,Betomsk}, ["Alfa","Beto"], ["Full", "50D"], true);
title(T, "Reduced Dimension Evolution Trajectory Comparison (End gen success P<0.01)")
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth_windiv_NOYLIM",figh,["png","pdf"])
for i=1:6,nexttile(T,i);ylim([0,1]);end
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth_windiv",figh,["png","pdf"])

%% Plot optimization trajectory separated by visual area and animal, summary only, no individual traces.
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {Alfamsk,Betomsk}, ["Alfa","Beto"], ["Full", "50D"], false);
title(T, "Reduced Dimension Evolution Trajectory Comparison (End gen success P<0.01)")
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth",figh,["png","pdf"])


function S = optimtraj_integral(score_traj2cmp, S)
% Integrate the area under the optimization trajectory.
% Return a structure. Collect various kinds of stats for 2 trajectories
% comparison.
% 
% score_traj2cmp: cell array of scores, 2-by-blocknum. 
%      in each cell is a 1-by-N array of single trial image scores. 
% this function will clip the last generation for stability. 
% No need to clip the score trajectory before entering.
if nargin <= 1, S=struct(); end
if isempty(score_traj2cmp), % use this to generate nan default dict
    S.("traj_int") = nan(1, 2);
    S.("traj_int_norm") = nan(1, 2);
    S.("traj_int_ratio") = nan;
    S.("traj_incr_int") = nan(1, 2);
    S.("traj_incr_int_norm") = nan(1, 2);
    S.("traj_incr_int_ratio") = nan;
    return; 
end
score_m_arr = cellfun(@mean, score_traj2cmp(:,1:end-1));
traj_int = sum(score_m_arr,2)'; % AUC without subtracting 1st gen scores
traj_int_norm = mean(score_m_arr,2)';
traj_incr_int = sum(score_m_arr-score_m_arr(:,1),2)'; % AUC with subtracting 1st gen scores 
traj_incr_int_norm = mean(score_m_arr-score_m_arr(:,1),2)'; % Area for only increased activation.
S.("traj_int") = traj_int;
S.("traj_int_norm") = traj_int_norm;
S.("traj_int_ratio") = traj_int(2) / traj_int(1); % ratio between 2 AUC
S.("traj_incr_int") = traj_incr_int;
S.("traj_incr_int_norm") = traj_incr_int_norm;
S.("traj_incr_int_ratio") = traj_incr_int(2) / traj_incr_int(1);
end

function S = Dprime_integral(score_traj2cmp, S)
% Integrate D prime between 2 threads. 
% 
% score_traj2cmp: cell array of scores, 2-by-blocknum. 
%      in each cell is a 1-by-N array of single trial image scores. 
% this function will clip the last generation for stability. 
% No need to clip the score trajectory before entering.
if nargin <= 1, S=struct(); end
if isempty(score_traj2cmp), % use this to generate nan default dict
    S.("Dpr_int") = nan;
    S.("Dpr_int_norm") = nan;
    return; 
end
Dpr_arr = [];
for blocki = 1:size(score_traj2cmp,2)-1 
Dpr_arr(blocki) = computeCohen_d(score_traj2cmp{1,blocki}, score_traj2cmp{2,blocki});
end
Dpr_int = sum(Dpr_arr);
Dpr_int_norm = sum(Dpr_arr)/numel(Dpr_arr);
S.("Dpr_int") = Dpr_int;
S.("Dpr_int_norm") = Dpr_int_norm;
end

function S = score_cmp_stats(score2cmp, prefix, S)
% Comparing cell array of firing rates and write stats in S
% 
% score2cmp: a cell array of scores, 2 rows, each row is a thread / optimizer
% prefix: prefix to name the fields of the struct.
% S: Struct containing the stats, if given then write the stats into it; if not, create a new one.  
%
% Example: 
% S = score_cmp_stats(score_vec_col(:,2:3), "init23");
% S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23");
if nargin == 1, prefix=""; end
if nargin <= 2, S=struct(); end
if isempty(score2cmp), % use this to generate nan default dict
    S.(prefix+"_cmp_t") = nan;
    S.(prefix+"_cmp_p") = nan;
    S.(prefix+"_cmp_dpr") = nan;
    S.(prefix+"_mean") = nan(1,2);
    S.(prefix+"_sem") = nan(1,2);
    S.(prefix+"_m_ratio") = nan;
return; 
end
score_thr1 = cat(1,score2cmp{1,:});
score_thr2 = cat(1,score2cmp{2,:});
[~,P,~,STS] = ttest2(score_thr1, score_thr2);
S.(prefix+"_cmp_t") = STS.tstat;
S.(prefix+"_cmp_p") = P;
S.(prefix+"_cmp_dpr") = computeCohen_d(score_thr1, score_thr2);
S.(prefix+"_mean") = [mean(score_thr1), mean(score_thr2)];
S.(prefix+"_sem") = [sem(score_thr1), sem(score_thr2)];
S.(prefix+"_m_ratio") = mean(score_thr2) / mean(score_thr1);
end

function [zscore_vec_col] = zscore_cellarr(score_vec_col)
score_all_vec = cat(1,score_vec_col{:});
scoreM = mean(score_all_vec);
scoreS = std(score_all_vec);
zscore_vec_col = cellfun(@(vec)(vec-scoreM) / scoreS, score_vec_col, 'uni', 0);
end
