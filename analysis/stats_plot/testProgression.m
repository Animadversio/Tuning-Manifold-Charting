function title_str=testProgression(StatsTab_sum, statname, masks, labels, sepvarnm, titstr)
% A high level API for testing the Progression of one variable (statname) in table
% w.r.t. another categorical variable (sepvar). the sepvar is given in a
% series of masks e.g. {V1msk, V4msk, ITmsk}
% Currently it does Spearman correlation, ANOVA and linear regression model.
% and then generate and return the summary string. 
% 
% Example:
%  testProgression(RDEvolTab, "Dpr_int_norm", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
%      "Both Monk All Exp");
%  testProgression(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
%      "Both Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
title_str = compose("Test progression of %s ~ %s, for %s\n",statname, sepvarnm, titstr);
valuevec = [];
predvec = [];
for i = 1:numel(masks)
    value = StatsTab_sum.(statname)(masks{i});
    valuevec = [valuevec; value];
    predvec = [predvec; ones(numel(value),1)*(i-1)];
    title_str = title_str+compose("%s %.3f+-%.3f (n=%d)\t",labels(i),mean(value),sem(value),numel(value));
end
[cval, pval] = corr(predvec, valuevec, 'Type', 'Spearman');
corrstr = compose("Spearman Correlation between %s ~ %s %.3f (%.1e)\n",sepvarnm,statname,cval,pval);
lm = fitlm(predvec, valuevec);
[F_p,F_tbl] = anova1(valuevec,predvec,'off');
Fval = F_tbl{2,5}; F_df = F_tbl{4,3};
anovastr = compose("\nANOVA F=%.3f p=%.1e(df=%d,%d)\n",Fval,F_p,F_tbl{2,3},F_tbl{3,3});
lmstr = compose("Linear Regres %s = %.3f + %s * %.3f \n Intercept %.3f+-%.3f, Slope %.3f+-%.3f\n Slope!=0: T=%.1f P=%.1e\n Fit Rsquare=%.3f\n",...
                statname, lm.Coefficients.Estimate(1), sepvarnm, lm.Coefficients.Estimate(2),...
                lm.Coefficients.Estimate(1), lm.Coefficients.SE(1), ...
                lm.Coefficients.Estimate(2), lm.Coefficients.SE(2), ...
                lm.Coefficients.tStat(2), lm.Coefficients.pValue(2), ...
                lm.Rsquared.Ordinary);
title_str = title_str + anovastr + corrstr + lmstr;
% disp(lm)
disp(title_str);
end