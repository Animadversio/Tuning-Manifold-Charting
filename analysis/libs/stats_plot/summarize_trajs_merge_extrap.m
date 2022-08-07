function [h,score_m,score_s,blockvec] = summarize_trajs_merge_extrap(traj_col,norm_mode,movmean_N,...
    major_msk,major_lab,minor_msk,minor_lab,figdir,fignm,plot_indivi_traj)
% msk on the experimental session level. 
% Merge the trajectories in each tile. 
% Same logic as `summarize_trajs_tile`
% Directly adapted from `summarize_trajs_merge` but do extrapolation to
%   each traces to avoid noise after exp finished.
% Inputs: 
%   traj_col: 1d cell array of trajectories, containing mean rate for each block
%   msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
%   anim_msks = {Alfamsk&validmsk, Betomsk&validmsk};
%   label_col = ["V1", "V4", "IT"];
%   anim_col = ["Alfa","Beto"];
%   norm_mode: How to normalize single trajectory. Choose from ["max", "rng", "movmean_rng"], default to be "max". 
%   movmean_N: Number of blocks to do moving average. default to be 3. 
% 
if nargin<2, norm_mode="max";end
if nargin<3, movmean_N=3;end
if nargin<4, major_msk={ones(size(traj_col,1),1,'logical')}; major_lab=["all"]; end
if nargin<6 || isempty(minor_msk), 
    minor_msk={ones(size(traj_col,1),1,'logical')}; 
end
if nargin<7 || isempty(minor_lab), minor_lab = ["all"]; end 
if nargin<8, figdir=""; end
if nargin<9, fignm=compose("%sNorm_scoreTraj_avg_Area_Anim_movmean_merge",norm_mode); end
if nargin<10, plot_indivi_traj=false; end
Corder = colororder();

% Normalize the trajs 
if strcmp(norm_mode,"max") % normalize each traj to the max in that exp. 
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
elseif strcmp(norm_mode,"movmean_max") % normalize each traj to max of smooth traj
normed_traj_col = [];
for i = 1:numel(traj_col)
smooth_traj = movmean(traj_col{i}, 3);
MAX = max(smooth_traj);
normed_traj_col{i} = smooth_traj / MAX; %smooth_traj / MAX;%traj_col{i} / MAX
end
normed_traj_col = reshape(normed_traj_col, size(traj_col));
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

% extrapolate the trajs using the mean 
gennum = cellfun(@numel,normed_traj_col);
max_gen = max(gennum);
norm_extrp_traj_col = cell(size(normed_traj_col));
for Expi = 1:numel(normed_traj_col)
    finact = mean(normed_traj_col{Expi}(end-1:end)); % average over the last 2 full blocks
    gen2append = max_gen - gennum(Expi); 
    norm_extrp_traj_col{Expi} = [normed_traj_col{Expi}, ones(1,gen2append)*finact];
end

% block cell array
extrp_block_col = cellfun(@(traj)1:numel(traj), norm_extrp_traj_col, 'uni',0);
nrow = numel(minor_msk); ncol = numel(major_msk);
h=figure; %fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean", Animal);
set(h,'pos',[125,   258,   310*nrow, 100+310])
T = tiledlayout(1,nrow,"pad",'compact',"tilespac",'compact');
for animi=1:nrow
traj_lab = strings(1,ncol);
for mski=1:ncol
nexttile(T,animi);hold on
msk = major_msk{mski} & minor_msk{animi};
if plot_indivi_traj
    trajids = reshape(find(msk),1,[]);
    for i = trajids
        ind_cutoff = gennum(i);
        plot(extrp_block_col{i}(1:ind_cutoff), norm_extrp_traj_col{i}(1:ind_cutoff),'Color',[Corder(2*mski-1,:),0.4],'LineWidth',0.4,'HandleVis','off')
        plot(extrp_block_col{i}(ind_cutoff:end), norm_extrp_traj_col{i}(ind_cutoff:end),':','Color',[Corder(2*mski-1,:),0.4],'LineWidth',0.4,'HandleVis','off')
    end
end
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(extrp_block_col(msk), norm_extrp_traj_col(msk));%normed_traj_col
score_C_col_mov_m = movmean(score_C_col_m, movmean_N);
score_C_col_mov_s = movmean(score_C_col_s, movmean_N);
cutoff = maxk(gennum(msk),2); 
cutoff = cutoff(2);
shadedErrorBar(block_colvec(1:cutoff),score_C_col_mov_m(1:cutoff),score_C_col_mov_s(1:cutoff),'lineProps',{'Color',Corder(2*mski-1,:),'LineWidth',1.5})
shadedErrorBar(block_colvec(cutoff:end),score_C_col_mov_m(cutoff:end),score_C_col_mov_s(cutoff:end),'lineProps',{':','Color',Corder(2*mski-1,:),'LineWidth',1.5,'HandleVis','off'},...
                            'patchSaturation',0.05)
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
    XLIMIT = max(XLIMIT, prctile(gennum(msk),80));
end
xlim([1,XLIMIT]);
end
saveallform(figdir,fignm+"_Xlim");
end