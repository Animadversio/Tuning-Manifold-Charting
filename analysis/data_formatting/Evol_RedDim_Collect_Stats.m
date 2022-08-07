%% RedDim Collect Stats
Animal = "Alfa"; Set_Path;
result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_ReducDim";
expftr = contains(ExpRecord.Exp_collection, "ReducDimen_Evol") & ...
        contains(ExpRecord.expControlFN,"generate") &...
        ExpRecord.Expi>0;
rowis = find(expftr);
% [meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr),Animal);
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(rowis,Animal,false,true);
%%
RDStats = repmat(struct(), 1, length(meta_new));
%%
% The key analysis is how a unit responds to both optimizer thread
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
if isnan(Expi) || ~contains(ExpRecord.expControlFN{exp_rowi},'generate') ...
        || ~contains(ExpRecord.Exp_collection{exp_rowi},'ReducDimen_Evol')
    keyboard
    continue
end
fprintf("\nProcessing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
% save the most basic stuffs
RDStats(Expi).Animal = Animal;
RDStats(Expi).Expi = Expi;
RDStats(Expi).imageName = Trials.imageName;
RDStats(Expi).meta = meta;
%% Compute the less basi info about units, store them 
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                       unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
assert(all(unit_in_pref_chan==unit_in_pref_chan(1)), "threads should evolve to the same units in RD exp")
assert(all(pref_chan_id==pref_chan_id(1)), "threads should evolve to the same unit in RD exp")
assert(all(pref_chan==pref_chan(1)), "threads should evolve to the same channel in RD exp")

RDStats(Expi).units.pref_chan = pref_chan(1);
RDStats(Expi).units.pref_chan_id = pref_chan_id(1);
RDStats(Expi).units.unit_in_pref_chan = unit_in_pref_chan(1);
RDStats(Expi).units.unit_name_arr = unit_name_arr;
RDStats(Expi).units.unit_num_arr = unit_num_arr;
RDStats(Expi).units.activ_msk = activ_msk;
RDStats(Expi).units.spikeID = meta.spikeID;

record_id = find(meta.spikeID==pref_chan(1));
RDStats(Expi).units.rec_chan_id = record_id;
% Record the optimizer settings
if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
RDStats(Expi).evol.optim_names = Optim_names;
RDStats(Expi).evol.optim_opts = Trials.TrialRecord.User.optim_opts;
RDStats(Expi).evol.thread_num = thread_num;
RDStats(Expi).evol.imgpos = imgpos;
RDStats(Expi).evol.imgsize = imgsize;
RDStats(Expi).evol.unit_in_pref_chan = unit_in_pref_chan; 
RDStats(Expi).evol.pref_chan = pref_chan;

% sort the images
% seperate the thread natural images and generated images 
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 (and maybe thread 3 4)
thread_msks = cell(1, thread_num);
for threadi = 1:thread_num
    msk = contains(imgnm, compose("thread%03d", threadi - 1));
    thread_msks{threadi} = msk; % store masks in a structure for the ease to iterate
end
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm)) % images comes from all these threads
RDStats(Expi).stim.gen_msk = row_gen;
RDStats(Expi).stim.nat_msk = row_nat;
RDStats(Expi).stim.thread_msks = thread_msks;

% split blocks and put images into arrays
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
color_seq = brewermap(block_num, 'spectral');
RDStats(Expi).color_seq = color_seq;
RDStats(Expi).evol.block_arr = block_arr;
RDStats(Expi).evol.block_n = block_num; % a block number for each trial
% Save the indices and the corresponding PSTHs
gen_idx_seq = cell(thread_num, block_num);
nat_idx_seq = cell(thread_num, block_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        gen_idx_seq{threadi, blocki} = find(gen_msk);
        nat_idx_seq{threadi, blocki} = find(nat_msk);
    end
end
% Note the pref_chan_id here could be the same for 2 channels
gen_psth_col = cellfun(@(idx) rasters(record_id, :, idx), gen_idx_seq, ...
                'UniformOutput', false);
nat_psth_col = cellfun(@(idx) rasters(record_id, :, idx), nat_idx_seq, ...
                'UniformOutput', false);

RDStats(Expi).evol.idx_seq = gen_idx_seq;
RDStats(Expi).evol.psth = gen_psth_col;
RDStats(Expi).ref.idx_seq = nat_idx_seq;
RDStats(Expi).ref.psth = nat_psth_col;
%% Sort the reference images by name and find their psth, this part doesn't split thread
nat_imgnm = cellfun(@(thr_msk)unique(imgnm(row_nat & thr_msk)),thread_msks, 'UniformOutput', false); % sorted names of reference images 
nat_imgnm = cat(2,nat_imgnm{:});
nat_imgidx = cellfun(@(nm) find(contains(imgnm, nm)), nat_imgnm, 'UniformOutput', false);
nat_psth = cellfun(@(idx) rasters(record_id, :, idx), nat_imgidx, 'UniformOutput', false);

RDStats(Expi).ref.imgnm = string(nat_imgnm);
RDStats(Expi).ref.idx_arr = nat_imgidx; % Save the idx of pref chan for reference image, by name
RDStats(Expi).ref.psth_arr = nat_psth; % Save the psth of pref chan for reference image, by name

end
%%
%save("D:\Alfa_Evol_stats.mat", 'EStats')
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(compose("D:\\%s_RDEvol_stats.mat", Animal), 'RDStats')
save(fullfile(MatStats_path, compose("%s_RDEvol_stats.mat", Animal)), 'RDStats')