function [score_m,score_s,blockvec] = summarize_paired_trajs_tile(traj_col,major_msk,minor_msk)
% msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
% anim_msks = {Alfamsk&validmsk, Betomsk&validmsk};
% label_col = ["V1", "V4", "IT"];
% anim_col = ["Alfa","Beto"];
if nargin<2, norm_mode="max";end
if nargin<3, movmean_N=3;end
if nargin<4, major_msk={ones(size(traj_col,1),1,'logical')}; major_lab=["all"]; end
if nargin<6 || isempty(minor_msk), 
    minor_msk={ones(size(traj_col,1),1,'logical')}; 
end
if nargin<7 || isempty(minor_lab), minor_lab = ["all"]; end 
if nargin<8, figdir=""; end
if nargin<9, fignm=compose("%sNorm_scoreTraj_avg_Area_Anim_movmean",norm_mode); end
% Normalize the trajs 
if strcmp(norm_mode,"max") % normalize each traj to the max in that exp. 
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
elseif strcmp(norm_mode,"rng") % normalize each traj, min to 0, max to 1
normed_traj_col = cellfun(@(traj) (traj-min(traj))/(max(traj)-min(traj)),traj_col,'uni',0);
elseif strcmp(norm_mode,"movmean_rng") % normalize each traj to min max of smooth traj
normed_traj_col = [];
for i = 1:numel(traj_col)
smooth_traj = movmean(traj_col{i},3);
MIN = min(smooth_traj); 
MAX = max(smooth_traj);
normed_traj_col{i} = (traj_col{i}-MIN)/(MAX-MIN);
end
normed_traj_col = reshape(normed_traj_col, size(traj_col));
else
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
end
% block 
block_col = cellfun(@(traj)1:numel(traj),traj_col,'uni',0);
[score_m,score_s,blockvec] = sort_scoreblock(block_col, normed_traj_col)

h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean", Animal);
set(h,'pos',[125   258   935   715])
nrow = numel(minor_msk); ncol = numel(major_msk);
T = tiledlayout(nrow,ncol,"pad",'compact',"tilespac",'compact');
for animi=1:nrow
for mski=1:ncol
nexttile(T,mski+3*(animi-1))
msk = major_msk{mski} & minor_msk{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_col(msk),cat(2,score_C_m_trajs{msk}));
% [score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_col(msk),cat(2,score_G_m_trajs{msk}));
score_C_col_mov_m = movmean(score_C_col_m,3);
% score_G_col_mov_m = movmean(score_G_col_m,3);
score_C_col_mov_s = movmean(score_C_col_s,3);
% score_G_col_mov_s = movmean(score_G_col_s,3);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2,:)})
% shadedErrorBar(block_colvec,score_G_col_mov_m,score_G_col_mov_s,'lineProps',{'Color',Corder(1,:)})
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
end