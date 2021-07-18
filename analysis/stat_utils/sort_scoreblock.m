function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% sort an array of scores according to the block array labels. compute the
% mean and std for each block. 
% If `blockarr` and `scorearr` are cell array it will be concatenated to be an array
% else, they will be array or vector. 
% 
% Note: really useful function to summarize multiple evolution trajectories into
% a mean one. 
% Originally developped in the Evol_Optimizer_CMAGA_cmp.m
if iscell(blockarr) && iscell(scorearr)
blockarr = cellfun(@(S)reshape(S,[],1),blockarr,'uni',0);
scorearr = cellfun(@(S)reshape(S,[],1),scorearr,'uni',0);
blockarr = cat(1, blockarr{:});
scorearr = cat(1, scorearr{:});
end
assert (numel(blockarr)==numel(scorearr))
blockvec = min(blockarr):max(blockarr);
score_m = [];score_s = [];
for blocki = min(blockarr):max(blockarr)
    score_m(blocki) = mean(scorearr(blockarr==blocki));
    score_s(blocki) = sem(scorearr(blockarr==blocki));
end
end