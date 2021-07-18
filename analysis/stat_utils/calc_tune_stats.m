function reportStats = calc_tune_stats(psth_cells, wdw_vect)
% similar function to calc_tuning_stats, but with simpler input format: just a cell
% array of **psth** and it can compute T and F for you. 
% Input format:
%   psth_cells: cell array of psth array. each cell is a condition. Eahc
%       psth array is of shape (N, T, R) N neuron, T time, R trials. 
%       These axes order is the same as in `rasters`
%       It assumes the N==1, and use the first index to get PSTH. 
%   wdw_vect: Default to be empty. If it's not empty,  `wdw_vect` can be an 
%       array of time windows to compute the F stats for each time window
% 
if nargin == 1
    wdw_vect = [];
    MovF = false; % calculate F for a single time windows 50:200
else
    MovF = true; % calculate F for moving time windows
end
% PSTH Chacteristics F to baseline, t of baseline 
reportStats = struct();
groupsize = cellfun(@(psth) size(psth,3), psth_cells);
indices = reshape(1:numel(psth_cells),size(psth_cells));
idx_vect = arrayfun(@(L, idx) idx*ones(L,1), groupsize, indices, 'uni', false);

act_wdw = 51:200; bsl_wdw=1:50; % default value for activation and baseline window
score_vect = cellfun(@(psth)squeeze(mean(psth(1,act_wdw,:),2)),psth_cells,'uni',false); % non-baseline subtracted activation.
basel_vect = cellfun(@(psth)squeeze(mean(psth(1,bsl_wdw,:),2)),psth_cells,'uni',false);
idx_vect = cell2mat(reshape(idx_vect,[],1));
score_vect = cell2mat(reshape(score_vect,[],1));
basel_vect = cell2mat(reshape(basel_vect,[],1));
[P,ANOVATAB,STATS] = anova1(score_vect,idx_vect,'off');
reportStats.F = ANOVATAB{2,5};
reportStats.F_P = P;
[P,ANOVATAB,STATS] = anova1(basel_vect,idx_vect,'off');
reportStats.F_bsl = ANOVATAB{2,5};
reportStats.F_P_bsl = P;
[H,P,CI,STATS] = ttest2(score_vect,basel_vect);
reportStats.T = STATS.tstat;
reportStats.t_P = P;
if MovF
F_wdw = [];
F_P_wdw = [];
for fi = 1:size(wdw_vect,1)
    movscore_vect = cellfun(@(psth)squeeze(mean(psth(1,wdw_vect(fi,1):wdw_vect(fi,2),:),2)),psth_cells,'uni',false);
    movscore_vect = cell2mat(reshape(movscore_vect,[],1));
    [P,ANOVATAB,STATS] = anova1(movscore_vect,idx_vect,'off');
    if isempty(ANOVATAB{2,5}) % in case anova1 returns nothing as F.
        F_wdw(end+1) = nan;
        F_P_wdw(end+1) = nan;
    else
        F_wdw(end+1) = ANOVATAB{2,5};
        F_P_wdw(end+1) = P;
    end
end
reportStats.F_wdw = F_wdw;
reportStats.F_P_wdw = F_P_wdw;
end
end