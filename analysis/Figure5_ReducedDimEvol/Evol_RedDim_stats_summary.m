%% Manifold paper: Reduced Dimension Evolution vs Full space Evolution
%% Summary statistics plot: 
%      Load summary statistics for trajectories and plot them.
%  Replicate Analysis from Table and Stats
global figdir
figdir = "O:\Evol_ReducDim\summary";
ExpType = "RDEvol";
Animal = "Both"; Set_Path;
load(fullfile(figdir, Animal+"_RDEvol_summaryStats.mat"), "RDEvol_Stats")
RDEvolTab = readtable(fullfile(figdir, Animal+"_RDEvol_trajcmp.csv"));

%% Create Masks for testing 
validmsk = ones(numel(RDStats), 1, 'logical');
for iTr=1:numel(RDStats)
if ~all(RDStats(iTr).evol.optim_names == ["ZOHA Sphere lr euclid", "ZOHA Sphere lr euclid ReducDim"])
validmsk(iTr) = false;
end
end
% any thread succeeds vs both thread succeeds mask.
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
