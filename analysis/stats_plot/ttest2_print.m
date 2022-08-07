function [tval,pval,sumstr,mean_arr,sem_arr] = ttest2_print(group1,group2,g1label,g2label,paired)
if nargin==2, paired=false;g1label="Group 1";g2label="Group 2";end
if nargin==3, paired=g1label;g1label="Group 1";g2label="Group 2";end
if nargin==4, paired=false;end
g1_mean = nanmean(group1);
g2_mean = nanmean(group2);
g1_sem = sem(group1);
g2_sem = sem(group2);
mean_arr = [g1_mean,g2_mean];
sem_arr = [sem(group1),sem(group2)];
if ~ paired
[~,P,CI,TST] = ttest2(group1,group2);
else
[~,P,CI,TST] = ttest(group1,group2);
end
tval = TST.tstat;
pval = P;
sumstr = sprintf("%s (%.3f+-%.3f n=%d) - %s (%.3f+-%.3f n=%d): t=%.3f(df=%d), P=%.1e, CI=[%.1f,%.1f]\n",...
                g1label,g1_mean,g1_sem,numel(group1),...
                g2label,g2_mean,g2_sem,numel(group2),...
                tval,TST.df,P,CI(1),CI(2));
fprintf(sumstr)
end