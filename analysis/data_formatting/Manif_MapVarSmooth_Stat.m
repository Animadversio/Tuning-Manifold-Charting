%% 
%  Collect the Baseline and Activation trial by trial variability for Manifold
%  Experiments
Animal = "Alfa";Set_Path;
mat_dir = "O:\Mat_Statistics";
figdir = "O:\ImMetricTuning";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "selectivity")&...
            ExpRecord.Expi > 0);
[meta_new,rasters_new,~,Trials_new] = loadExperiments(find(expftr),Animal);
%%
MapVarStats = repmat(struct(), 1, length(Stats));
%%
for Triali = 1:length(Stats)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);Expi=Expi(end); % hack this for beto exp 35
if isnan(Expi) || ~all(contains(ExpRecord.expControlFN(exp_rowi),'selectivity')) ...
        || ~all(contains(ExpRecord.Exp_collection(exp_rowi),'Manifold'))
    % add this filter to process a sequence of Trials_new 
    keyboard
    continue
end
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
% savepath = fullfile(result_dir, sprintf("%s_Exp%02d", Animal, Expi));
% mkdir(savepath)
% if Expi == 16, keyboard; end
MapVarStats(Expi).Animal = Animal;
MapVarStats(Expi).Expi = Expi;
MapVarStats(Expi).imageName = Stats(Expi).imageName;
MapVarStats(Expi).meta = meta;
%%
% pref_chan = Trials.TrialRecord.User.prefChan;
% Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
% unit_name_arr = generate_unit_labels(meta.spikeID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
% unit_in_pref_chan = 0;
pref_chan = Stats(Expi).units.pref_chan;
unit_num_arr = Stats(Expi).units.unit_num_arr;
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                unit_num_arr > 0); % match for unit number
MapVarStats(Expi).units = Stats(Expi).units;

if Animal == "Beto"
    if Expi <= 10, subsp_n = 3; else, subsp_n = 1;end
elseif Animal == "Alfa"
    subsp_n = 1;
end
% subsp_templ = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d','norm_%d_RND1_%d_RND2_%d'};
% sphere_norm = infer_norm_from_imgnm(Trials.imageName);
MapVarStats(Expi).manif.sphere_norm = Stats(Expi).manif.sphere_norm;
MapVarStats(Expi).manif.subsp_n = subsp_n;
for subsp_i = 1:subsp_n
idx_grid = Stats(Expi).manif.idx_grid{subsp_i};
activ_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx), 2)), idx_grid, 'Uni', false);
basel_col = cellfun(@(idx) squeeze(mean(rasters(:, 1:45, idx), 2)), idx_grid, 'Uni', false);
MapVarStats(Expi).manif.act_col{subsp_i} = activ_col;
MapVarStats(Expi).manif.bsl_col{subsp_i} = basel_col;
end
%%
MapVarStats(Expi).ref.didGabor = Stats(Expi).ref.didGabor;
MapVarStats(Expi).ref.didPasu = Stats(Expi).ref.didPasu;
if Stats(Expi).ref.didPasu
pasu_idx_grid = Stats(Expi).ref.pasu_idx_grid;
pasu_activ_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),2)), pasu_idx_grid, 'Uni', false);
pasu_basel_col = cellfun(@(idx) squeeze(mean(rasters(:, 1:45, idx),2)), pasu_idx_grid, 'Uni', false);
MapVarStats(Expi).ref.pasu_act_col = pasu_activ_col;
MapVarStats(Expi).ref.pasu_bsl_col = pasu_basel_col;
end

if Stats(Expi).ref.didGabor
gab_idx_grid = Stats(Expi).ref.gab_idx_grid;
pasu_activ_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),2)), gab_idx_grid, 'Uni', false);
pasu_basel_col = cellfun(@(idx) squeeze(mean(rasters(:, 1:45, idx),2)), gab_idx_grid, 'Uni', false);
MapVarStats(Expi).ref.gab_act_col = activ_col;
MapVarStats(Expi).ref.gab_bsl_col = basel_col;
end
end
%%
save(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')