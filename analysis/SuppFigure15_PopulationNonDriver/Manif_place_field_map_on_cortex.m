%% Manifold tuning plotted on cortex for all exps
% Using MapVarStats instead of original data to work 
% Basic structure is kept and renovate this code 
% global  Trials rasters channel sphere_norm ang_step Reps
% storedStruct = load("D:\\Manifold_Exps.mat");
Animal = "Beto"; Set_Path; 
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"), 'EStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')

%% 
result_dir = "O:\\CorticTuneMaps";
mkdir(result_dir)
%% Set up canvas for plotting
figIT = figure(9);clf;hold on
figV1 = figure(10);clf;hold on
figV4 = figure(11);clf;hold on
[ax_arrA,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4);
figITB = figure(12);clf;hold on
figV1B = figure(13);clf;hold on
figV4B = figure(14);clf;hold on
[ax_arrB,tITB,tV1B,tV4B] = Cortex_Channel_Tile_Layout_All(figITB, figV1B, figV4B);
%%
% constant for maximum number of repetitions (as long as it's larger than the maximum, it's fine)
Reps = 11;
sign_filter = true;
ang_step = 18;
for Expi = 1:numel(Stats)
% assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
fprintf("Processing %s Exp %d, %s\n", Animal, Expi, Stats(Expi).meta.comments)

% sphere_norm = norm_arr(Expi);
pref_chan = EStats(Expi).units.pref_chan;%pref_chan_arr(Expi);
unit_name_arr = MapVarStats(Expi).units.unit_name_arr;
unit_num_arr = MapVarStats(Expi).units.unit_num_arr;
activ_msk = MapVarStats(Expi).units.activ_msk;
spikeID = MapVarStats(Expi).units.spikeID;
imgpos = EStats(Expi).evol.imgpos;
imgsize = EStats(Expi).evol.imgsize;
savepath = sprintf("O:\\PC_space_tuning\\%s_Exp%d_chan%02d", Animal, Expi, pref_chan);
mkdir(savepath);

% assert(sum(contains(unit_name_arr, "C"))==0, "Too many units in a channel! FIX THE CODE TO HANDLE ME! ")
% Select the significant channels
si = 1;
act_col = MapVarStats(Expi).manif.act_col{si};
bsl_col = MapVarStats(Expi).manif.bsl_col{si};
act_tsr = cell2mat(reshape(cellfun(@(A)mean(A,2),act_col,'uni',0),[1,11,11]));
bsl_tsr = cell2mat(reshape(cellfun(@(A)mean(A,2),bsl_col,'uni',0),[1,11,11]));
FStats = [];
for iCh = 1:numel(spikeID)
    actmap_col = cellfun(@(A)squeeze(A(iCh,:)),act_col,'uni',0);
    FStat = anova_cells(actmap_col);
    FStats = [FStats, FStat];
end
FTab = struct2table(FStats);
% load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
signif_chan = find(FTab.F_P<1E-2);
% rsp_mat_arr = {};
signif_F = FTab.F(signif_chan);
signif_P = FTab.F_P(signif_chan);
fprintf("=======\nExp %02d Significantly modulated channels:\n", Expi)
% for channel = 1:numel(spikeID)
%     chan_label_str = sprintf("Channel %s  ", unit_name_arr{channel});
%     if Stat_summary{channel,1}.anova_p < 0.01
%         signif_chan(end+1) = channel;
%         signif_F(end+1) = Stat_summary{channel,1}.anova_F;
%         signif_P(end+1) = Stat_summary{channel,1}.anova_p;
%         [score_mat, bsl_mat] = get_score_mat('norm_%d_PC2_%d_PC3_%d');
%         rsp_mat_arr{end+1} = score_mat;
%         fprintf("%s\t", unit_name_arr{channel})
%     else
        
%     end
% end
fprintf('\n')
fprintf(num2str(signif_F,'%.2f\t'))
fprintf("\n%d units in total \n", numel(signif_chan))
%% 
Exp_label = compose("Exp %d Evolv channel %d PC23 space %.1f deg cent [%.1f %.1f]", Expi, pref_chan, ...
    imgsize, imgpos(1), imgpos(2));
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:64, spikeID, activ_msk); % add active mask to get rid of bad zero units. 
if all(chan_idxA==chan_idxB) % each channel has only 1 unit exactly
    plot2unit = false;
else % some channels have more than 1 unit
    plot2unit = true;
end
for ax = [ax_arrA,ax_arrB]
cla(ax{1},'reset');%colorbar(ax{1},'off')
end
title(tIT,strcat(Exp_label, " IT array"))
title(tV1,strcat(Exp_label, " V1V2 array"))
title(tV4,strcat(Exp_label, " V4 array"))
if plot2unit
    title(tITB,strcat(Exp_label, " IT array"))
    title(tV1B,strcat(Exp_label, " V1V2 array"))
    title(tV4B,strcat(Exp_label, " V4 array"))
end
for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan); % the 
    if isnan(channel) % if the channel has no unit. 
        set(ax_arrA{arr_chan},'Visible','off')
        set(ax_arrB{arr_chan},'Visible','off')
        continue
    end
    chan_label_str = compose("Ch %s F%.2f(p=%.1e)", unit_name_arr{channel}, ...
        FTab.F(channel), FTab.F_P(channel));
    % [score_mat, bsl_mat] = get_score_mat('norm_%d_PC2_%d_PC3_%d');
    plot_contour_heatmap(squeeze(act_tsr(channel,:,:)), ax_arrA{arr_chan})
    title([chan_label_str])
    
    if plot2unit
        channel = chan_idxB(arr_chan); % the c
        chan_label_str = compose("Ch %s F%.2f(p=%.1e)", unit_name_arr{channel}, ...
            FTab.F(channel), FTab.F_P(channel));
        % [score_mat, bsl_mat] = get_score_mat('norm_%d_PC2_%d_PC3_%d');
        plot_contour_heatmap(squeeze(act_tsr(channel,:,:)), ax_arrB{arr_chan})
        title([chan_label_str])
    end
    
end
% suptitle(Exp_label)
ITimnm = compose("IT_array_%s_Exp%02d_PC23_placemap", Animal, Expi);
V1imnm = compose("V1_array_%s_Exp%02d_PC23_placemap", Animal, Expi);
V4imnm = compose("V4_array_%s_Exp%02d_PC23_placemap", Animal, Expi);
saveas(figIT, fullfile(savepath, ITimnm+"A.jpg"))
saveas(figIT, fullfile(result_dir, ITimnm+"A.jpg"))
saveas(figIT, fullfile(result_dir, ITimnm+"A.pdf"))
saveas(figV1, fullfile(savepath, V1imnm+"A.jpg"))
saveas(figV1, fullfile(result_dir, V1imnm+"A.jpg"))
saveas(figV1, fullfile(result_dir, V1imnm+"A.pdf"))
saveas(figV4, fullfile(savepath, V4imnm+"A.jpg"))
saveas(figV4, fullfile(result_dir, V4imnm+"A.jpg"))
saveas(figV4, fullfile(result_dir, V4imnm+"A.pdf"))
if plot2unit
saveas(figITB, fullfile(savepath, ITimnm+"B.jpg"))
saveas(figITB, fullfile(result_dir, ITimnm+"B.jpg"))
saveas(figITB, fullfile(result_dir, ITimnm+"B.pdf"))
saveas(figV1B, fullfile(savepath, V1imnm+"B.jpg"))
saveas(figV1B, fullfile(result_dir, V1imnm+"B.jpg"))
saveas(figV1B, fullfile(result_dir, V1imnm+"B.pdf"))
saveas(figV4B, fullfile(savepath, V4imnm+"B.jpg"))
saveas(figV4B, fullfile(result_dir, V4imnm+"B.jpg"))
saveas(figV4B, fullfile(result_dir, V4imnm+"B.pdf"))
end
%% if sign filter then Plot only the significant channels 
if sign_filter
    for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan);
    if isnan(channel) % The channel is not in the data! 
        continue
    end
    if ~ (FTab.F_P(channel) < 0.01) % Note not < is different from > because there are NaN p values
        set(get(ax_arrA{arr_chan},'Children'),'Visible','off');
        colorbar(ax_arrA{arr_chan},'off')
    end
    if plot2unit
    channel = chan_idxB(arr_chan);
    if ~ (FTab.F_P(channel) < 0.01)
        set(get(ax_arrB{arr_chan},'Children'),'Visible','off');
        colorbar(ax_arrB{arr_chan},'off')
    end
    end
    end
    saveas(figIT, fullfile(savepath, ITimnm+"_signifA.jpg"))
    saveas(figIT, fullfile(result_dir, ITimnm+"_signifA.jpg"))
    saveas(figIT, fullfile(result_dir, ITimnm+"_signifA.pdf"))
    saveas(figV1, fullfile(savepath, V1imnm+"_signifA.jpg"))
    saveas(figV1, fullfile(result_dir, V1imnm+"_signifA.jpg"))
    saveas(figV1, fullfile(result_dir, V1imnm+"_signifA.pdf"))
    saveas(figV4, fullfile(savepath, V4imnm+"_signifA.jpg"))
    saveas(figV4, fullfile(result_dir, V4imnm+"_signifA.jpg"))
    saveas(figV4, fullfile(result_dir, V4imnm+"_signifA.pdf"))
    if plot2unit
    saveas(figITB, fullfile(savepath, ITimnm+"_signifB.jpg"))
    saveas(figITB, fullfile(result_dir, ITimnm+"_signifB.jpg"))
    saveas(figITB, fullfile(result_dir, ITimnm+"_signifB.pdf"))
    saveas(figV1B, fullfile(savepath, V1imnm+"_signifB.jpg"))
    saveas(figV1B, fullfile(result_dir, V1imnm+"_signifB.jpg"))
    saveas(figV1B, fullfile(result_dir, V1imnm+"_signifB.pdf"))
    saveas(figV4B, fullfile(savepath, V4imnm+"_signifB.jpg"))
    saveas(figV4B, fullfile(result_dir, V4imnm+"_signifB.jpg"))
    saveas(figV4B, fullfile(result_dir, V4imnm+"_signifB.pdf"))
    end
end

end
%%
function plot_contour_heatmap(score_mat_avg, ax)
    % Given a averaged score matrix, plot a contour heatmap to certain
    % axis! 
    maximum = max(score_mat_avg, [], 'all');
    axes(ax);
    imagesc(-90:18:90, -90:18:90, score_mat_avg) % sum(score_mat,3)./cnt_mat
    colormap('parula')
    hold on 
    % Add contour to the plot for visibility
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.9] * maximum,...
        'LineColor',[1 0 0],'LineWidth',1) % Note single value will cause it to interpret that value as contour line number! 
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.75] * maximum,...
        'LineColor',[0.6350 0.0780 0.1840],'LineWidth',1) 
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.5] * maximum,...
        'LineColor','w','LineWidth',1) 
    %contourcmap('hot')
    %ylabel("PC 2 degree");xlabel("PC 3 degree")
    axis off;shading flat;axis image;colorbar
end
