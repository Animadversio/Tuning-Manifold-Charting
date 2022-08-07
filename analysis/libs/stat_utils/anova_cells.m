function [reportStats] = anova_cells(score_cell)
reportStats	= struct();
groupsize = cellfun(@(S) numel(S), score_cell); % how many elements under that condition
indices = reshape(1:numel(score_cell),size(score_cell)); % index for that condition
idx_cell = arrayfun(@(L, idx) idx*ones(L,1), groupsize, indices, 'uni', 0); % create index cell array.
score_cell = cellfun(@(S)reshape(S,[],1), score_cell,'uni',0);
assert(all(size(idx_cell) == size(score_cell)))
idx_vect = cell2mat(reshape(idx_cell,[],1)); 
score_vect = cell2mat(reshape(score_cell,[],1)); 
[P,ANOVATAB,STATS] = anova1(score_vect,idx_vect,'off');
reportStats.F = ANOVATAB{2,5};
reportStats.F_P = P;
reportStats.STATS = STATS;
reportStats.ANOVATAB = ANOVATAB;
end