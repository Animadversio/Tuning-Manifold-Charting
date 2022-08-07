%% Evol_collect_Stats (for Manifold Exp mainly)
% Collect statisitcs from evolution experiments into a compressed compilation
% Trealy well designed, Note this code support multi-thread well! 
% Really succint and well written code! 
clearvars -except EStats Stats meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug  ExpSpecTable_Aug_alfa ExpRecord
%%
Animal = "Alfa";Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
% expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
%             contains(ExpRecord.expControlFN, "generate"));
expftr = (contains(ExpRecord.expControlFN,"200429"));
rowis = find(expftr);
% Expi_col = [1,2,3,6,7,10,11,12,19,27];
% Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
% assert(all(ExpRecord.Expi(rowis(Expi_col))==Expi_col')) % assert you are getting what you want. 
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(rowis, Animal);
%%
EStats = repmat(struct(), 1, length(meta_new));
%%
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
%% 
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
if isnan(Expi) || ~contains(ExpRecord.expControlFN{exp_rowi},'generate') ...
        || ~contains(ExpRecord.Exp_collection{exp_rowi},'Manifold')
    keyboard
    continue
end
fprintf("\nProcessing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
% savepath = fullfile(result_dir, sprintf("%s_Exp%02d", Animal, Expi));
% mkdir(savepath)
% if Expi == 16, keyboard; end
EStats(Expi).Animal = Animal;
EStats(Expi).Expi = Expi;
EStats(Expi).imageName = Trials.imageName;
EStats(Expi).meta = meta;

pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num % note here we use the unit_num_arr which exclude the null channels.
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
EStats(Expi).units.pref_chan = pref_chan;
EStats(Expi).units.unit_name_arr = unit_name_arr;
EStats(Expi).units.unit_num_arr = unit_num_arr;
EStats(Expi).units.activ_msk = activ_msk;
EStats(Expi).units.spikeID = meta.spikeID;
EStats(Expi).units.pref_chan_id = pref_chan_id;
EStats(Expi).units.unit_in_pref_chan = unit_in_pref_chan;

if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
EStats(Expi).evol.optim_names = Optim_names;
EStats(Expi).evol.thread_num = thread_num;
EStats(Expi).evol.imgpos = imgpos;
EStats(Expi).evol.imgsize = imgsize;
EStats(Expi).evol.unit_in_pref_chan = unit_in_pref_chan; 
EStats(Expi).evol.pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,1));

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
EStats(Expi).stim.gen_msk = row_gen;
EStats(Expi).stim.nat_msk = row_nat;
EStats(Expi).stim.thread_msks = thread_msks;

% split blocks and put images into arrays
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
color_seq = brewermap(block_num, 'spectral');
EStats(Expi).color_seq = color_seq;
EStats(Expi).evol.block_arr = block_arr;
EStats(Expi).evol.block_n = block_num;
gen_idx_seq = cell(thread_num, block_num); % generated image idx cell as a horizontal array. 
nat_idx_seq = cell(thread_num, block_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        gen_idx_seq{threadi, blocki} = find(gen_msk);
        nat_idx_seq{threadi, blocki} = find(nat_msk);
    end
end
gen_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), gen_idx_seq, 'Uni', 0);
nat_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), nat_idx_seq, 'Uni', 0);
EStats(Expi).evol.idx_seq = gen_idx_seq;
EStats(Expi).evol.psth = gen_psth_col;
EStats(Expi).ref.idx_seq = nat_idx_seq;
EStats(Expi).ref.psth = nat_psth_col;

% Sort the reference images by name and find their psth. 
nat_imgnm = unique(imgnm(row_nat)); % sorted names of reference images 
nat_imgidx = cellfun(@(nm) find(contains(imgnm, nm)), nat_imgnm, 'Uni', 0);
nat_psth = cellfun(@(idx) rasters(pref_chan_id, :, idx), nat_imgidx, 'Uni', 0);

EStats(Expi).ref.imgnm = string(nat_imgnm);
EStats(Expi).ref.idx_arr = nat_imgidx; % Save the idx of pref chan for reference image, by name
EStats(Expi).ref.psth_arr = nat_psth; % Save the psth of pref chan for reference image, by name
end
%%
%save("D:\Alfa_Evol_stats.mat", 'EStats')
MatStats_path = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(compose("D:\\%s_Evol_stats.mat", Animal), 'EStats')
save(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
%% Add the ref image stats to it
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats','-v6')
