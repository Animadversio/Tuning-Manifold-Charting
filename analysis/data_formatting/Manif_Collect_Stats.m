%% Manif_Collect_Stats into mat
%%
Animal = "Beto";Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "selectivity")&...
            ExpRecord.Expi > 0);
rowis = find(expftr);
% Expi_col = [1,2,3,6,7,10,11,12,19,27];
% Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
% assert(all(ExpRecord.Expi(rowis(Expi_col))==Expi_col')) % assert you are getting what you want. 
[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowis,Animal); % ,false, true
%%
Stats = repmat(struct(), 1, length(meta_new));
%% If there is stats saved, load it! 
% mat_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%%
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
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
Stats(Expi).Animal = Animal;
Stats(Expi).Expi = Expi;
Stats(Expi).imageName = Trials.imageName;
Stats(Expi).meta = meta;
%%
pref_chan = Trials.TrialRecord.User.prefChan;
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
% unit_in_pref_chan = 0;
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                unit_num_arr > 0); % match for unit number

Stats(Expi).units.pref_chan = pref_chan;
Stats(Expi).units.unit_name_arr = unit_name_arr;
Stats(Expi).units.unit_num_arr = unit_num_arr;
Stats(Expi).units.activ_msk = activ_msk;
Stats(Expi).units.spikeID = meta.spikeID;
Stats(Expi).units.pref_chan_id = pref_chan_id;

if Animal == "Beto"
    if Expi <= 10, subsp_n = 3; else, subsp_n = 1;end
elseif Animal == "Alfa"
    subsp_n = 1;
end
subsp_templ = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d','norm_%d_RND1_%d_RND2_%d'};
sphere_norm = infer_norm_from_imgnm(Trials.imageName);
Stats(Expi).manif.sphere_norm = sphere_norm;
Stats(Expi).manif.subsp_n = subsp_n;
for subsp_i = 1:subsp_n
name_pattern = subsp_templ{subsp_i};
[idx_grid, ~,~,~] = build_idx_grid(Trials.imageName, name_pattern, sphere_norm);
%%
psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), idx_grid, ...
                    'UniformOutput', false);
% activ_map = cellfun(@(c) mean(c(1,50:200,:), [2,3]), psths_col, 'UniformOutput', true);
Stats(Expi).manif.idx_grid{subsp_i} = idx_grid;
Stats(Expi).manif.psth{subsp_i} = psths_col;
end
%%
Stats(Expi).ref.didGabor = false;
Stats(Expi).ref.didPasu = false;
if sum(contains(Trials.imageName, "pasu")) > 186
Stats(Expi).ref.didPasu = true;
[pasu_idx_grid,~,~,~] = build_Pasu_idx_grid(Trials.imageName);
pasu_psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), pasu_idx_grid, ...
                'UniformOutput', false);
Stats(Expi).ref.pasu_psths = pasu_psths_col;
Stats(Expi).ref.pasu_idx_grid = pasu_idx_grid;
end

if sum(contains(Trials.imageName, "gab")) > 12
Stats(Expi).ref.didGabor = true;
[gab_idx_grid,~,~,~] = build_Gabor_idx_grid(Trials.imageName);
gab_psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), gab_idx_grid, ...
                'UniformOutput', false);
Stats(Expi).ref.gab_psths = gab_psths_col;
Stats(Expi).ref.gab_idx_grid = gab_idx_grid;
end
end
%%
savepath = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(compose("D:\\%s_Manif_stats.mat", Animal), 'Stats')
save(fullfile(savepath, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%save("D:\Alfa_Manif_stats.mat", 'Stats')
%%
savepath = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(fullfile(savepath, compose("%s_Manif_stats.mat", Animal)), 'Stats','-v6')
%%
[pasu_idx_grid,id_grid,p1_grid,p2_grid] = build_Pasu_idx_grid(Trials.imageName);
pasu_psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), pasu_idx_grid, ...
                'UniformOutput', false);
uniti=1;
act_mat = cellfun(@(psth) squeeze(mean(psth(uniti, 51:200, :), [2,3])), pasu_psths_col, ...
                'UniformOutput', true);
act_col = cellfun(@(psth) squeeze(mean(psth(uniti, 51:200, :), [2])), pasu_psths_col, ...
                'UniformOutput', false);
anovan()

%%
[pasu_idx_grid,~,~,~] = build_Pasu_idx_grid(Trials.imageName);
[gab_idx_grid,~,~,~] = build_Gabor_idx_grid(Trials.imageName);
Stats(Expi).ref.pasu_idx_grid = pasu_idx_grid;
Stats(Expi).ref.gab_idx_grid = gab_idx_grid;
%%
savefast(fullfile(mat_dir, Animal+'_Manif_stats.mat'),'Stats')
%%
% psths_mean = cellfun(@(c) mean(c, 3), psths_col, 'UniformOutput', false);
% psths_1 = cellfun(@(c) squeeze(c(1,:)), psths_mean, 'UniformOutput', false);
% %%
% activ_map = cellfun(@(c) mean(c(1,50:200,:), [2,3]), psths_col, 'UniformOutput', true);
% figure
% imagesc(activ_map)
% axis equal tight; 
% %%
% figure
% for i = 1:size(idx_grid, 1)
% for j = 1:size(idx_grid, 2)
% subplottight(11,11,j+(i-1)*11);
% plot(psths_1{i,j})
% end
% end
%%
function [idx_grid, id_mat, p1_mat, p2_mat] = build_idx_grid(imgnms, name_pattern, sphere_norm)
ang_step = 18;
idx_grid = cell(11,11);
id_mat = nan(11, 11);
p1_mat = nan(11, 11);
p2_mat = nan(11, 11);
for i =  -5:5 % i code for PC2 
    for j =  -5:5 % j code for PC3
        cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
        img_idx = find(contains(imgnms, cur_fn));
        if j ~= 5 && j ~= -5
            id = 11 * i + j;
        elseif j == 5 
            id = 5;
        elseif j == -5
            id = -5;
        end
        id_mat(i+6, j+6) = id; 
        p1_mat(i+6, j+6) = i*ang_step;
        p2_mat(i+6, j+6) = j*ang_step;
        idx_grid{i+6, j+6} = img_idx; 
    end
end
end

function [idx_grid, id_mat, p1_mat, p2_mat] = build_Pasu_idx_grid(imgnms)
idx_grid = cell(51, 4);
id_mat = nan(51, 4);
p1_mat = nan(51, 4);
p2_mat = nan(51, 4);
for i = 1:51
    for j = 1:4
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(imgnms, cur_fn));
        if i==1 || i==2 || i==3
            id = 4 * (i-1) + 1;
        else 
            id = 4 * (i-1) + j;
        end
        id_mat(i, j) = id;
        p1_mat(i, j) = i;
        p2_mat(i, j) = 2*j-1;
        idx_grid{i, j} = img_idx;
    end
end
end

function [idx_grid, id_mat, p1_mat, p2_mat] = build_Gabor_idx_grid(imgnms)
ori = 0:30:150; sf = [0.5, 1];
idx_grid = cell(2,6);
id_mat = nan(2, 6);
p1_mat = nan(2, 6);
p2_mat = nan(2, 6);
for i = 1:numel(sf) 
    for j =  1:numel(ori) 
        cur_fn = sprintf('gab_ori_%.1f_%.1f', ori(j), sf(i));
        img_idx = find(contains(imgnms, cur_fn));
        id = 6 * (i - 1) + j;
        id_mat(i, j) = id; 
        idx_grid{i, j} = img_idx; 
    end
end
end


function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern, sphere_norm)
    global  Trials rasters channel ang_step Reps
    score_mat = nan(11,11,Reps); 
    bsl_mat = nan(11,11,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(11,11); 
    id_mat = zeros(11,11); % record the id correspond to i,j
    for i =  -5:5 % i code for PC2 
        for j =  -5:5 % j code for PC3
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
            baseline = squeeze(mean(psths(1, 1:50, :)));
            score_mat(i+6, j+6, 1:length(img_idx)) = scores; % first idx is PC2. 
            bsl_mat(i+6, j+6, 1:length(img_idx)) = baseline;
            if j ~= 5 && j ~= -5
                id = 11 * i + j;
            elseif j == 5 
                id = 5;
            elseif j == -5
                id = -5;
            end
            id_mat(i+6, j+6) = id; 
        end
    end
    [theta_mat, phi_mat] = meshgrid(ang_step*(-5:5), ang_step*(-5:5));
    mean_fr_mat = bsl_mat + score_mat;
    id_vec_nan = reshape(repmat(id_mat, 1,1, Reps), 1, []);
    score_vec_nan = reshape(score_mat, 1, []);
    bsl_vec_nan = reshape(bsl_mat, 1, []);
    mean_fr_vec_nan = bsl_vec_nan + score_vec_nan;
    % Do statistics
    [p,tbl,stats] = anova1(score_vec_nan, id_vec_nan, 'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    %
    [p2,tbl2,stats2] = anovan(score_vec_nan, {reshape(repmat(theta_mat, 1,1,Reps),1,[]), ...
                              reshape(repmat(phi_mat, 1,1,Reps),1,[])}, 'model', 'interaction','display' ,'off');
    stats2.p = p2;
    stats2.F = [tbl2{2:4,6}];
    summary.anova2_p = p2; % p for theta, phi and interaction 
    summary.anova2_F = [tbl2{2:4,6}]; % F for theta, phi and interaction
    %
    [~,P,CI] = ttest(mean_fr_vec_nan, bsl_vec_nan);
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],CI(1), CI(2), P, stats.F, stats.p, ...
            stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));

end