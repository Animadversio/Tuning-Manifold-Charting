function [tval,pval,sumstr,mean_arr,sem_arr] = ttest2corr_print(group1,group2,g1label,g2label)
if nargin==2, g1label="Group 1";g2label="Group 2";end
% g1_mean = nanmean(group1);
% g2_mean = nanmean(group2);
g1_mean = tanh(nanmean(fisherR2Z(group1)));
g2_mean = tanh(nanmean(fisherR2Z(group2)));
mean_arr = [g1_mean,g2_mean];
sem_arr = [sem(group1),sem(group2)];
[~,P,CI,TST] = ttest2(fisherR2Z(group1),fisherR2Z(group2));
tval = TST.tstat;
pval = P;
sumstr = sprintf("%s (%.3f+-%.3f n=%d) - %s (%.3f+-%.3f n=%d): t=%.3f(df=%d), P=%.1e, CI=[%.1f,%.1f]\n",g1label,g1_mean,sem_arr(1),numel(group1),g2label,g2_mean,sem_arr(2),numel(group2),tval,TST.df,P,CI(1),CI(2));
fprintf(sumstr)
end