function [tval,pval,sumstr,mean_arr,sem_arr] = ttest_print(group1,g1label)
% Signature:
%   ttest_print(group1,g1label)
%    
% Examples:
%    ttest_print( val_arr );
if nargin==1, g1label="Group 1"; end
g1_mean = nanmean(group1);
g1_sem = sem(group1);
mean_arr = [g1_mean];
sem_arr = [g1_sem];
[~,P,CI,TST] = ttest(group1);
tval = TST.tstat; pval = P;
sumstr = sprintf("%s (%.3f+-%.3f) - 0: t=%.3f(df=%d), P=%.1e, CI=[%.2f,%.2f]",g1label,g1_mean,sem_arr,tval,TST.df,pval,CI(1),CI(2));
fprintf("%s\n",sumstr)
end