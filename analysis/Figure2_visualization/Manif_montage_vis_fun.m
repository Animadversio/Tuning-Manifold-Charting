% This plot function is used to get the color framed montage plot for any
% manifold experiments for the use in paper. 
% Plot part of this function is obtained from PC_space_Pasu_tuning_analysis.m

Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";

% for Animal = ["Alfa", "Beto"]
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"),'MapVarStats')

%%
Animal = "Alfa"; Expi=4; 
% channum=61;chan2plot = find(any(MapVarStats(Expi).units.spikeID==[channum],2));
pref_chan = Stats(Expi).units.pref_chan;
chan2plot = find(any(MapVarStats(Expi).units.spikeID==[pref_chan],2));
%%
si=1; 
figroot = "E:\OneDrive - Washington University in St. Louis\PC_space_tuning";
savedir = fullfile(figroot,compose("%s_Exp%d_chan%02d", Animal, Expi, pref_chan));
tic
manif_imgfn = cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'uni',0);
[manif_imgfps, mapper] = map2fullpath(manif_imgfn, Stats(Expi).meta.stimuli);
manif_images = cellfun(@(fp)imread(fp),manif_imgfps,'uni',0);
fprintf("Load images took %.2f sec\n",toc)
manif_imgtile = imtile(manif_images','BorderSize',2);
imwrite(manif_imgtile,fullfile(savedir,"PC23_manif_imgtile.jpg"))
%% chan2plot = chan2plot;%find(any(MapVarStats(Expi).units.spikeID==[pref_chan],2));
for iCh = reshape(chan2plot,1,[])
tic
unitstr = Stats(Expi).units.unit_name_arr{iCh};
chan_label_str = sprintf("%s Exp%d (Driver %d) Channel %s", Animal, Expi, pref_chan, unitstr);

actvec = cell2mat(reshape(cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.act_col{si},'uni',0),1,[]));
actmap = cellfun(@(A)mean(A(iCh,:),'all'),MapVarStats(Expi).manif.act_col{si},'uni',1);
actcol = cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.act_col{si},'uni',0);
bslvec = cell2mat(reshape(cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.bsl_col{si},'uni',0),1,[]));
bslmean = mean(bslvec);
bslstd = std(bslvec);
scoremap = actmap - bslmean;
[~,P,CI,TST] = ttest(actvec, bslvec);
tval = TST.tstat; pval = P;
anovaS = anova_cells(actcol(:));
stat_str = compose(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)'],CI(1), CI(2), P, anovaS.F, anovaS.F_P);
tic
figure(3); clf; set(3,'pos',[ 805         197        1559         781]); % all manifold images montaged
set(0,'CurrentFigure',3); clf;T=tiledlayout(1,2,'tilesp','compact','padd','compact');
ax1 = nexttile(T,1);%subplot(1,2,1);
imagesc(-90:18:90, -90:18:90, scoremap) 
ylabel("PC 2 degree");xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
axis image;colorbar%shading flat;
% Get the colormap from the heatmap in order to visualize the map below. 
frame_img_list = score_frame_image_arr(manif_images, scoremap...
    , caxis(ax1), colormap(ax1), 16);
ax2 = nexttile(T,2); %subplot(1,2,2);
montage(frame_img_list', 'Size', [11, 11]);
fprintf("Fully visualize images took %.2f sec\n",toc)
saveallform(savedir,compose("PC23_tune_chan%s.pdf",unitstr),3,["pdf"])
% saveas(3,fullfile(savedir,compose("PC23_tune_chan%s.pdf",unitstr)))
% set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
% set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
%%
frame_imgtile = imtile(frame_img_list');
imwrite(frame_imgtile,fullfile(savedir,compose("PC23_tune_chan%s_scored_imgs.jpg",unitstr)))
end
winopen(savedir)