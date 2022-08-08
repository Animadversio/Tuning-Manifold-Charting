function [figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   mskcol, labcol, mskrow, labrow, thread_labels, plot_indiv, Xlim)
% Visualize paired Evolution experiments (CMA vs GA; Full vs RedDim; BigGAN vs FC6)
%   Form a tile grid, separate the tiles by 0,1,2 variables one in row, one in col.
% Plotting based on output of `optim_traj_process` function (newer API)
% 
% Input:
%  score_m_traj_extrap_col : 
%  block_traj_extrap_col : 
%  extrap_mask_col : 
%  mskcol : cell array of masks to plot across columns. e.g. `{V1msk,V4msk,ITmsk}`. 
%           If `isempty(mskcol)`, then there is one column, no separating variable.
%  labcol : string array, labels for each masks. e.g. `["V1","V4","IT"]`
%  mskrow : cell array of masks to plot across rows. e.g. `{Alfamsk,Betomsk}`
%           If `isempty(mskrow)`, then there is one row, no separating variable.
%  labrow : string array, labels for each masks. e.g. `["Alfa","Beto"]`
%  thread_labels : cell arry of string, or string array. 
%                usually len = 2. Example: 
%                       `thread_labels = ["4096d", "50d"];`
%  plot_indiv : Bool. if true, also plot individual traces; if false, just the summary average traj. 
%  Xlim : Bool. if true, use the `XLIM_Cutoff_Prctl` percentile of trajectory length as cutoff. 
%               if false, just plot all the trajs, so xlim will be the max of traj length. 
% 
% Example:
%   [score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
%      extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
%      optim_traj_process(block_traces, score_traces, ["Full", "50D"], "max12", 55);
%   [figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
%      {V1msk&sucsmsk_end,V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V1","V4","IT"], {Alfamsk,Betomsk}, ["Alfa","Beto"], ["Full", "50D"], true);
%   title(T, "Reduced Dimension Evolution Trajectory Comparison (End gen success P<0.01)")
%   saveallform([figdir,outdir],"Both_MaxNorm_extrapSmth_windiv_NOYLIM",figh,["png","pdf"])
% 
if nargin<10, Xlim = true;end
expN = numel(block_traj_extrap_col);
if isempty(mskcol), mskcol={ones(expN,1,'logical')}; labcol="";
else, mskcol = cellfun(@(msk)reshape(msk,[],1),mskcol,'uni',0); end
if isempty(mskrow), mskrow={ones(expN,1,'logical')}; labrow="";
else, mskrow = cellfun(@(msk)reshape(msk,[],1),mskrow,'uni',0); end
XLIM_Cutoff_Prctl = 75;
Corder = colororder(); % by default the first condition is orange, 2nd is blue.
blockN_arr = cellfun(@(M)sum(~M), extrap_mask_col); 
ncol = numel(mskcol);
nrow = numel(mskrow);
figh = figure('position', [200,100,300*ncol,300*nrow]);
T = tiledlayout(nrow,ncol,'TileSp','compact','Pad','compact');
for ci = 1:ncol
   for ri = 1:nrow
   ax = nexttile(T,ci+ncol*(ri-1));hold on
   msk = mskcol{ci} & mskrow{ri};
   [score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(msk),score_m_traj_extrap_col(msk,1));
   [score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(msk),score_m_traj_extrap_col(msk,2));
   % [score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_extrap_col{msk}),cat(2,score_m_traj_extrap_col{msk,1}));
   % [score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_extrap_col{msk}),cat(2,score_m_traj_extrap_col{msk,1}));
   shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:),'LineWidth',1.5})
   shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:),'LineWidth',1.5})
   title(ax,compose("%s %s (N=%d)",labrow(ri),labcol(ci),sum(msk)))
    if plot_indiv
        for iTr=1:expN
          if ~msk(iTr),continue; end % only plot those within the mask, jump
          blkmsk = extrap_mask_col{iTr};
          datablocks = block_traj_extrap_col{iTr}(~blkmsk); % block numbers with real data  
          extrapblocks = block_traj_extrap_col{iTr}(blkmsk); % block numbers without. (extrapolation)
          % if plot individual traces it perform 3 blocks mov mean smoothing by default. 
          score_trajs = {movmean(score_m_traj_extrap_col{iTr,1}(~blkmsk),3),...
                           movmean(score_m_traj_extrap_col{iTr,2}(~blkmsk),3)};
          % show the real data with solid line
          plot(datablocks, score_trajs{1},'Color',[Corder(2,:),0.2],'LineWidth',0.5)
          plot(datablocks, score_trajs{2},'Color',[Corder(1,:),0.2],'LineWidth',0.5)
          % show the constant extrapolation part with blkmsk with dashed line
          plot([datablocks(end),extrapblocks], [score_trajs{1}(end),score_m_traj_extrap_col{iTr,1}(blkmsk)],'Color',[Corder(2,:),0.2],'LineWidth',0.5,'LineStyle','-.')
          plot([datablocks(end),extrapblocks], [score_trajs{2}(end),score_m_traj_extrap_col{iTr,2}(blkmsk)],'Color',[Corder(1,:),0.2],'LineWidth',0.5,'LineStyle','-.')
      end
    end
    if Xlim
      Xmax = prctile(blockN_arr(msk),XLIM_Cutoff_Prctl);
      xlim([0,Xmax])
    end
   end
end
nexttile(T,1);legend(thread_labels,'location','best')
title(T,"Trajectory Compare")
end