function [score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
	extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
	optim_traj_process(block_traces, score_traces, thread_labels, norm_scheme, blockN_extrap)
% Stand alone preprocessing function that modularize the analsis done to: 
%   CMA vs GA; Full vs RedDim; BigGAN vs FC6 evolution comparison. 
%  
%  Comparing pairs of optimization trajectories, the key inputs are 
%  
% Input
%   block_traces: N by 2 cell array. Full block number vector per thread
%   score_traces: N by 2 cell array. Full Score vector per thread 
%   thread_labels: string array. label for thread. e.g. ["CMA", "GA"]; ["Full", "RedDim"]
%   norm_scheme: string. Selecting the normalization scheme for the
%      trajectories. Can be "max1","max2","max12"
% 
% Return 
%   score_m_traj_col
%   block_traj_col
%   score_m_traj_extrap_col
%   block_traj_extrap_col
%   extrap_mask_col
%   sucsmsk_end
%   sucsmsk_max
%   tval_end_arr, pval_end_arr, 
%   tval_max_arr, pval_max_arr
if nargin <= 3, norm_scheme = "max1"; end
if nargin <= 3, blockN_extrap = 50; end 
%% Normalize score by the max in one thread. 
expN = size(score_traces,1);

score_m_traj_col = cell(expN,2);
block_traj_col = cell(expN,1);
score_m_traj_extrap_col = cell(expN,2);
block_traj_extrap_col = cell(expN,1);
extrap_mask_col = cell(expN,1);
for iTr=1:expN
    if (isempty(block_traces{iTr,1}) || isempty(block_traces{iTr,2})) % invalid trial
        score_m_traj_col(iTr,:) = {[],[]};
        block_traj_col(iTr,:) = {[]};
        score_m_traj_extrap_col(iTr,:) = {[],[]};
        block_traj_extrap_col(iTr,:) = {[]};
        extrap_mask_col(iTr,:) = {[]};
    else
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{iTr,1},score_traces{iTr,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{iTr,2},score_traces{iTr,2});
    normmin = 0;
    if strcmp(norm_scheme,"max1")
    normmax = max(score_C_m);
    elseif strcmp(norm_scheme,"max2")
    normmax = max(score_G_m);
    elseif strcmp(norm_scheme,"max12")
    normmax = max(cat(2,score_C_m,score_G_m));
    elseif strcmp(norm_scheme,"init1")
    normmax = mean(score_C_m(1:2));
    elseif strcmp(norm_scheme,"init2")
    normmax = mean(score_G_m(1:2));
    elseif strcmp(norm_scheme,"init12")
    normmax = mean([score_C_m(1:2),score_G_m(1:2)]);
    else
        
    end
    scaling = abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_m_traj_col{iTr,1} = score_C_m_norm;
    score_m_traj_col{iTr,2} = score_G_m_norm;
    block_traj_col{iTr} = blockvec;
    % Flat extrapolation to the same gen
    extrap_val_C = mean(score_C_m_norm(end-1:end));
    extrap_val_G = mean(score_G_m_norm(end-1:end));
    extrap_N = blockN_extrap - numel(blockvec);
    if numel(blockvec) < blockN_extrap
        score_m_traj_extrap_col{iTr,1} = [score_C_m_norm,extrap_val_C * ones(1,extrap_N)];
        score_m_traj_extrap_col{iTr,2} = [score_G_m_norm,extrap_val_G * ones(1,extrap_N)];
        block_traj_extrap_col{iTr} = 1:blockN_extrap;
        extrap_mask_col{iTr} = (1:blockN_extrap) > numel(blockvec);
    else
        score_m_traj_extrap_col{iTr,1} = score_C_m_norm;
        score_m_traj_extrap_col{iTr,2} = score_G_m_norm;
        block_traj_extrap_col{iTr} = blockvec;
        extrap_mask_col{iTr} = ~ones(1,numel(blockvec),'logical');
    end
    end
end

%% Compare the first few and last few / max gen activation. 
tval_end_arr = []; pval_end_arr = [];
tval_max_arr = []; pval_max_arr = [];
for iTr = 1:size(score_traces, 1)
if (isempty(block_traces{iTr,1}) || isempty(block_traces{iTr,2})) % invalid trial
tval_end_arr(iTr,:) = nan; pval_end_arr(iTr,:) = nan;
tval_max_arr(iTr,:) = nan; pval_max_arr(iTr,:) = nan;
else
for ithread=1:2
fprintf(thread_labels(ithread)+" ")
blockN = max(block_traces{iTr,ithread});
[tval,pval] = ttest2_print(score_traces{iTr,ithread}(any(block_traces{iTr,ithread}==[1,2],2)),...
    score_traces{iTr,ithread}(any(block_traces{iTr,ithread}==[blockN-1,blockN],2)),"init12","last12");
tval_end_arr(iTr,ithread) = tval;
pval_end_arr(iTr,ithread) = pval;
[score_m,score_s,blockvec] = sort_scoreblock(block_traces{iTr,ithread},...
                score_traces{iTr,ithread});
[~,maxN] = max(score_m);
[tval,pval] = ttest2_print(score_traces{iTr,ithread}(any(block_traces{iTr,ithread}==[1,2],2)),...
    score_traces{iTr,ithread}(any(block_traces{iTr,ithread}==[maxN-1,maxN],2)),"init12","max12");
tval_max_arr(iTr,ithread) = tval;
pval_max_arr(iTr,ithread) = pval;
end
end
end
% masks of experiments that any one of it succeed
pThresh = 0.01;
sucsmsk_end = any((tval_end_arr<0)&(pval_end_arr<pThresh), 2);% <0 enforce that the evolution is going up, not down.
sucsmsk_max = any((tval_max_arr<0)&(pval_max_arr<pThresh), 2);
end