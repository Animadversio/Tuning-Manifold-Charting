%% Create the summary plots for Compare evolution scores for Alfa/Beto.
%  See if the Evolution works 
%  (Newer Refactored version using new plotting API comparing to Evol_RedDim_summary.m. May. 2021)
%  Plot the evolution trajectory comparison
%  With handy stats function below. 

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
%% 
RDEvol_Stats = []; % struct list, higher level summary, easy for plotting and testing.
unitnum_arr = zeros(length(RDStats),1);
prefchan_arr = zeros(length(RDStats),1);
validmsk = ones(numel(RDStats), 1, 'logical'); % whether the exp should be excluded. (unconventional optimizer setup)
score_traces = cell(numel(RDStats), 2); 
zscore_traces = cell(numel(RDStats), 2); 
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
% make sure the thread order
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
%% Collect stats and save
writetable(RDEvolTab, fullfile(figdir, Animal+"_RDEvol_trajcmp.csv"))
save(fullfile(figdir, Animal+"_RDEvol_summaryStats.mat"), "RDEvol_Stats")

%% Newer API verions
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["Full", "50D"], "max12", 55);
%%
figdir = "O:\Evol_ReducDim\summary";
outdir = "O:\Manuscript_Manifold\Figure3\RedDimEffectProg";
anim_arr = [RDStats.Animal]';
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),RDStats');
area_arr = arrayfun(@area_map,prefchan_arr);
V1msk = area_arr=="V1";
V4msk = area_arr=="V4";
ITmsk = area_arr=="IT";
Amsk = anim_arr=="Alfa";
Bmsk = anim_arr=="Beto";
anysucsmsk = any(RDEvolTab.t_p_succ<0.01,2);
allsucsmsk = all(RDEvolTab.t_p_succ<0.01,2);
%%
figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
	{V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {}, [], ["Full", "50D"], true);
%%
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["Full", "50D"], true);
title(T, "Reduced Dimension Evolution Trajectory Comparison (End gen success P<0.01)")
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth_windiv_NOYLIM",figh,["png","pdf"])
for i=1:6,nexttile(T,i);ylim([0,1]);end
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth_windiv",figh,["png","pdf"])
%%
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["Full", "50D"], false);
title(T, "Reduced Dimension Evolution Trajectory Comparison (End gen success P<0.01)")
saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth",figh,["png","pdf"])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trajectory Summary plot: 
%% Re-normalize the trajectories to show together.
score_C_m_trajs = {};
score_G_m_trajs = {};
score_C_s_trajs = {};
score_G_s_trajs = {};
block_trajs = {};
for Expi=1:size(score_traces,1)
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},score_traces{Expi,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},score_traces{Expi,2});
    normmin = 0; % (score_C_m(1)+score_G_m(1))/2;%0; %
    normmax = max(max(score_C_m),max(score_G_m));
    scaling = normmax;%abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_C_m_trajs{Expi} = score_C_m_norm;
    score_G_m_trajs{Expi} = score_G_m_norm;
    score_C_s_trajs{Expi} = score_C_s/scaling;
    score_G_s_trajs{Expi} = score_G_s/scaling;
    block_trajs{Expi} = blockvec;
end
%% Filtering Array
Animal_arr = struct2table(RDEvol_Stats).Animal;
Alfamsk = Animal_arr=="Alfa";
Betomsk = Animal_arr=="Beto";
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
anysucsmsk = any(RDEvolTab.t_p_succ<0.01,2);
allsucsmsk = all(RDEvolTab.t_p_succ<0.01,2);
%% Plot trajectory comparison for each individual session.
h=figure;hold on;fignm=compose("%s_MaxNorm_scoreTraj_indiv_all", Animal);
set(h,'pos',[1000         315         765         663])
for Expi=1:numel(block_trajs)
    if isempty(block_trajs{Expi}), continue;end
    shadedErrorBar(block_trajs{Expi},score_C_m_trajs{Expi},score_C_s_trajs{Expi},'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1},'patchSaturation',0.05)
    shadedErrorBar(block_trajs{Expi},score_G_m_trajs{Expi},score_G_s_trajs{Expi},'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1},'patchSaturation',0.05)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("%s All Sessions Evol Trajectory Comparison", Animal),'FontSize',14)
legend(["Full","50D"])
saveallform(figdir,fignm);
xlim([0,35]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);

%% Plot trajectory comparison averaging all sessions.
h=figure;hold on;fignm=compose("%s_MaxNorm_scoreTraj_avg_all", Animal);
set(h,'pos',[1000         315         765         663])
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{:}),cat(2,score_C_m_trajs{:}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{:}),cat(2,score_G_m_trajs{:}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("%s All Sessions Averaged Optimization Trajectory",Animal))
legend(["Full","50D"])
saveallform(figdir,fignm);
xlim([0,35]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT
msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT
msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_movmean", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
score_C_col_mov_m = movmean(score_C_col_m,3);
score_G_col_mov_m = movmean(score_G_col_m,3);
score_C_col_mov_s = movmean(score_C_col_s,3);
score_G_col_mov_s = movmean(score_G_col_s,3);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_mov_m,score_G_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3 % relimit x axis by 80 percentile of block number
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT in A B (Final version)
msk = validmsk&anysucsmsk;
msk_col = {V1msk&msk, V4msk&msk, ITmsk&msk};
anim_msks = {Alfamsk&msk, Betomsk&msk};
label_col = ["V1", "V4", "IT"];
anim_col = ["Alfa","Beto"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean_anysucs", Animal);
set(h,'pos',[125   258   935   715])
T = tiledlayout(2,numel(msk_col),"pad",'compact',"tilespac",'compact');
for animi=1:2
for mski=1:3
nexttile(T,mski+3*(animi-1))
msk = msk_col{mski} & anim_msks{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
score_C_col_mov_m = movmean(score_C_col_m,3);
score_G_col_mov_m = movmean(score_G_col_m,3);
score_C_col_mov_s = movmean(score_C_col_s,3);
score_G_col_mov_s = movmean(score_G_col_s,3);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_mov_m,score_G_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory\n%s %s (n=%d)",anim_col(animi),label_col(mski),sum(msk)))
legend(["Full","50D"])
end
end
title(T, compose("%s Summary of mean evol trajectory (3 block movmean) for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for animi=1:2
for mski=1:3
    nexttile(T,mski+3*(animi-1))
    msk = msk_col{mski} & anim_msks{animi};
    xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
end
saveallform(figdir,fignm+"_Xlim");

%% Area average with individual curves plotted on it.
% msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
msk = validmsk&anysucsmsk;
msk_col = {V1msk&msk, V4msk&msk, ITmsk&msk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_ws_indiv_Area_anysucs", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski);hold on
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:),'LineWidth',1.5})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:),'LineWidth',1.5})
for Expi = find(msk)' % plot individual experiments' optim traj in the mask. 
    plot(block_trajs{Expi},movmean(score_C_m_trajs{Expi},3),'Color',[Corder(2,:),0.4],'lineWidth',0.5)
    plot(block_trajs{Expi},movmean(score_G_m_trajs{Expi},3),'Color',[Corder(1,:),0.4],'lineWidth',0.5)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");

%%
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_indiv_Area", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
for Expi = find(msk)
    shadedErrorBar(block_trajs{Expi},score_C_m_trajs{Expi},score_C_s_trajs{Expi},'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1},'patchSaturation',0.08)
    shadedErrorBar(block_trajs{Expi},score_G_m_trajs{Expi},score_G_s_trajs{Expi},'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1},'patchSaturation',0.08)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Collect of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary statistics plot: Load summary statistics for trajectories and plot them.
%  Replicate Analysis from Table and Stats
figdir = "E:\OneDrive - Washington University in St. Louis\Evol_ReducDim\summary";
ExpType = "RDEvol";
Animal = "Both";
load(fullfile(figdir, Animal+"_RDEvol_summaryStats.mat"), "RDEvol_Stats")
RDEvolTab = readtable(fullfile(figdir, Animal+"_RDEvol_trajcmp.csv"));
%% Masks for testing 
validmsk = ones(numel(RDStats), 1, 'logical');
for iTr=1:numel(RDStats)
if ~all(RDStats(iTr).evol.optim_names == ["ZOHA Sphere lr euclid", "ZOHA Sphere lr euclid ReducDim"])
validmsk(iTr) = false;
end
end
anysucsmsk = any(RDEvolTab.t_p_succ<0.01,2);
allsucsmsk = all(RDEvolTab.t_p_succ<0.01,2);

Alfamsk = (RDEvolTab.Animal=="Alfa");
Betomsk = (RDEvolTab.Animal=="Beto");
V1msk = (RDEvolTab.pref_chan<=48 & RDEvolTab.pref_chan>=33);
V4msk = (RDEvolTab.pref_chan>48);
ITmsk = (RDEvolTab.pref_chan<33);
%% Test on individual Session and Collect T stats on population
%% Test on the aggregated mean value at population level with a T test. 
diary(fullfile(figdir,"progression_summary.log"))
msk = validmsk&anysucsmsk;%validmsk
fprintf("Inclusion criterion: Use the right Optimizer pair and any of the two threads succeed.\n")
testProgression(RDEvolTab, "Dpr_int_norm", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "middle_cmp_dpr", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "last23_cmp_dpr", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "last23_m_ratio", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "middle_m_ratio", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "traj_int_ratio", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");

diary off

%% Statistics Separate by Area. 
h = stripe_plot(RDEvolTab, "last23_m_ratio", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                    "Both Monk All Exp", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%
h = stripe_plot(RDEvolTab, "last23_cmp_dpr", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                    "Both Monk All Exp", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%
h = stripe_plot(RDEvolTab, "Dpr_int_norm", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                    "Both Monk All Exp", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%%

h = stripe_minor_plot(RDEvolTab, "last23_cmp_dpr", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
    {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp", "area_anim_sep", {[1,2],[2,3],[1,3]}, 'marker','MarkerEdgeAlpha',0.9);                

%% Statistics Separate by Area and Animal. 
msk = validmsk&anysucsmsk;%validmsk
h = stripe_minor_plot(RDEvolTab, "last23_cmp_dpr", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], ...
                    {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp (Any success)", "area_anim_sep_anysucs", {[1,2],[2,3],[1,3]}, 'marker','MarkerEdgeAlpha',0.9);                

h = stripe_minor_plot(RDEvolTab, "Dpr_int_norm", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], ...
                   {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp (Any success)", "area_anim_sep_anysucs", {[1,2],[2,3],[1,3]}, 'marker', 'MarkerEdgeAlpha',0.9);

h = stripe_minor_plot(RDEvolTab, "traj_int_ratio", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], ...
                   {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp (Any success)", "area_anim_sep_anysucs", {[1,2],[2,3],[1,3]}, 'marker', 'MarkerEdgeAlpha',0.9);

h = stripe_minor_plot(RDEvolTab, "last23_m_ratio", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1","V4","IT"], ...
                    {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp (Any success)", "area_anim_sep_anysucs", {[1,2],[2,3],[1,3]}, 'marker','MarkerEdgeAlpha',0.9);                

%%
% score2cmp = score_vec_col(:,2:3);
% S = score_cmp_stats(score2cmp, "init23");
% S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23", S);
% Dprime_integral(score_vec_col)

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
% Comparing cell array of firing rates and write stats 
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

% function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% % sort an array of scores according to the block array labels. compute the
% % mean and std for each block. 
% % really useful function to summarize multiple evolution trajectories into
% % a mean one. 
% blockvec = min(blockarr):max(blockarr);
% score_m = [];score_s = [];
% for blocki = min(blockarr):max(blockarr)
%     score_m(blocki) = mean(scorearr(blockarr==blocki));
%     score_s(blocki) = sem(scorearr(blockarr==blocki));
% end
% end
% % 
% function saveallform(figdir,fignm,h,sfxlist)
% % Save a (current) figure with all suffices in a figdir. 
% if nargin <=3, h=gcf; end
% if nargin <=4, sfxlist = ["fig","pdf","png"]; end
% for sfx = sfxlist
% if strcmp(sfx, "fig")
%    savefig(h,fullfile(figdir,fignm+"."+sfx))
% else
%    saveas(h,fullfile(figdir,fignm+"."+sfx))
% end
% end
% end