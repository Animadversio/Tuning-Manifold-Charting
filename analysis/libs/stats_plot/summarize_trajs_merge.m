function [h,score_m,score_s,blockvec] = summarize_trajs_merge(traj_col,norm_mode,movmean_N,major_msk,major_lab,minor_msk,minor_lab,figdir,fignm)
% msk on the experimental session level. 
% Merge the trajectories in each tile. 
% Same logic as `summarize_trajs_tile`
% traj_col: 1d cell array of trajectories
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
if nargin<9, fignm=compose("%sNorm_scoreTraj_avg_Area_Anim_movmean_merge",norm_mode); end
Corder = colororder();
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
else % By default, normalize each traj to the max in that exp. 
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
end
% block cell array
block_col = cellfun(@(traj)1:numel(traj),traj_col,'uni',0);
% [score_m,score_s,blockvec] = sort_scoreblock(block_col, normed_traj_col)
nrow = numel(minor_msk); ncol = numel(major_msk);
h=figure; %fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean", Animal);
set(h,'pos',[125,   258,   310*nrow, 100+310])
T = tiledlayout(1,nrow,"pad",'compact',"tilespac",'compact');
for animi=1:nrow
traj_lab = strings(1,ncol);
for mski=1:ncol
nexttile(T,animi);hold on
msk = major_msk{mski} & minor_msk{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_col(msk),normed_traj_col(msk));
score_C_col_mov_m = movmean(score_C_col_m,movmean_N);
score_C_col_mov_s = movmean(score_C_col_s,movmean_N);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2*mski-1,:)})
traj_lab(mski) = major_lab(mski) + compose(" n=%d",sum(msk));
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory\n%s %s cmp",minor_lab(animi),strjoin(major_lab)))
legend(traj_lab,'Location','best')
end
title(T, compose("Summary of mean evol trajectory (%d block movmean) for each area",movmean_N), 'FontSize',16)
saveallform(figdir,fignm);
for animi=1:nrow
nexttile(T,animi)
XLIMIT = 0;
for mski=1:ncol
    msk = major_msk{mski} & minor_msk{animi};
    XLIMIT = max(XLIMIT, prctile(cellfun(@numel,block_col(msk)),80));
end
xlim([1,XLIMIT]);
end
saveallform(figdir,fignm+"_Xlim");
end