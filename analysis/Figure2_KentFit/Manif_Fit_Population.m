%% Manifold paper. Kent function fitting for every single tuning map. 
%% Newer version to fit Kent function for all the channels in all Manifold Experiments
%  Result from this will be analyzed in PopVis and summary, 

Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_Fitting\popstats";
Animal="Alfa"; 
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
%%
addpath D:\Github\Fit_Spherical_Tuning
addpath e:\Github_Projects\Fit_Spherical_Tuning
%% Fitting Kent function without a baseline parameter. 
% pool = parpool(6);
T00 = tic;
parfor Expi = 1:numel(Stats) % parallelize at session level. 
% Basic information to collect into the Stats
nCh = numel(MapVarStats(Expi).units.spikeID);
Animal_tab = array2table(repmat(Animal,nCh,1),'VariableNames',{'Animal'});
Expi_tab = array2table(repmat(Expi,nCh,1),'VariableNames',{'Expi'});
unitstr_tab = array2table(MapVarStats(Expi).units.unit_name_arr,'VariableNames',{'unitstr'});
unit_tab = array2table(MapVarStats(Expi).units.unit_num_arr,'VariableNames',{'unitnum'});
chan_tab = array2table(MapVarStats(Expi).units.spikeID,'VariableNames',{'chan'});
prefchan_tab = array2table(repmat(MapVarStats(Expi).units.pref_chan,nCh,1),'VariableNames',{'prefchan'});
tic
% compute the activation maps for each channel. as a cell array
for si = 1:MapVarStats(Expi).manif.subsp_n
space_tab = array2table(repmat(si,nCh,1),'VariableNames',{'space'});
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
actmap_mean = arrayfun(@(iCh)cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats(Expi).manif.act_col{si},'uni',1),...
    1:nCh,'uni',0);
toc
% Apply Kent fitting and ANOVA to each cell i.e. each channel's tuning map
KentStats_exp = cellfun(@(actmap)fit_Kent_stats(actmap),actmap_mean);
anovaStats_exp = cellfun(@anova_cells,actmap_col);
toc
anovaStats_tab = struct2table(anovaStats_exp);
KentStats_tab = struct2table(KentStats_exp);
StatsTab = [Animal_tab, Expi_tab, unitstr_tab, unit_tab, chan_tab, prefchan_tab, space_tab, anovaStats_tab, KentStats_tab];
writetable(StatsTab,fullfile(tabdir,compose("%s_Exp%d_sp%d.csv",Animal,Expi,si))); % save results for each exp. 
% StatsTab_sum = [StatsTab_sum;StatsTab];
% disp(size(StatsTab_sum))
end
end
toc(T00);
% compile all exp from an animal into a big csv
StatsTab_sum = [];
for Expi = 1:numel(Stats)
for si = 1:MapVarStats(Expi).manif.subsp_n
curtab = readtable(fullfile(tabdir,compose("%s_Exp%d_sp%d.csv",Animal,Expi,si)));
StatsTab_sum = [StatsTab_sum;curtab];
end
end
writetable(StatsTab_sum,fullfile(tabdir,compose("%s_Exp_all_KentStat_pole.csv",Animal)));


%% Fitting Kent function with a baseline parameter. 
%  This part differs from the upper part only by that one toggle setting. 
% pool = parpool(6);
T00 = tic;
parfor Expi = 1:numel(Stats)
% Basic information to collect into the Stats
nCh = numel(MapVarStats(Expi).units.spikeID);
Animal_tab = array2table(repmat(Animal,nCh,1),'VariableNames',{'Animal'});
Expi_tab = array2table(repmat(Expi,nCh,1),'VariableNames',{'Expi'});
unitstr_tab = array2table(MapVarStats(Expi).units.unit_name_arr,'VariableNames',{'unitstr'});
unit_tab = array2table(MapVarStats(Expi).units.unit_num_arr,'VariableNames',{'unitnum'});
chan_tab = array2table(MapVarStats(Expi).units.spikeID,'VariableNames',{'chan'});
prefchan_tab = array2table(repmat(MapVarStats(Expi).units.pref_chan,nCh,1),'VariableNames',{'prefchan'});
tic
for si = 1:MapVarStats(Expi).manif.subsp_n
space_tab = array2table(repmat(si,nCh,1),'VariableNames',{'space'});
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
actmap_mean = arrayfun(@(iCh)cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats(Expi).manif.act_col{si},'uni',1),...
    1:nCh,'uni',0);
toc
KentStats_exp = cellfun(@(actmap)fit_Kent_stats(actmap,true),actmap_mean);
anovaStats_exp = cellfun(@anova_cells,actmap_col);
toc
anovaStats_tab = struct2table(anovaStats_exp);
KentStats_tab = struct2table(KentStats_exp);
StatsTab = [Animal_tab, Expi_tab, unitstr_tab, unit_tab, chan_tab, prefchan_tab, space_tab, anovaStats_tab, KentStats_tab];
writetable(StatsTab,fullfile(tabdir,compose("%s_Exp%d_sp%d_bsl.csv",Animal,Expi,si)));
% StatsTab_sum = [StatsTab_sum;StatsTab];
% disp(size(StatsTab_sum))
end
end
%% 
% compile all exp from an animal into a big csv
toc(T00);
StatsTab_sum = [];
for Expi = 1:numel(Stats)
for si = 1:MapVarStats(Expi).manif.subsp_n
curtab = readtable(fullfile(tabdir,compose("%s_Exp%d_sp%d_bsl.csv",Animal,Expi,si)));
StatsTab_sum = [StatsTab_sum;curtab];
end
end
writetable(StatsTab_sum,fullfile(tabdir,compose("%s_Exp_all_KentStat_bsl_pole.csv",Animal)));



function S = fit_Kent_stats(funcval,fitbsl)
if nargin==1, fitbsl = false; end
% Take the mean value for north pole and south pole. 
funcval(:,1) = mean(funcval(:,1));
funcval(:,end) = mean(funcval(:,end));
if fitbsl
[Parameter, gof] = fit_Kent_w_baseline(funcval);
else
[Parameter, gof] = fit_Kent(funcval);
end
S.sse = gof.sse;
S.R2 = gof.rsquare;
S.dfe = gof.dfe;
S.adjR2 = gof.adjrsquare;
S.rmse = gof.rmse;
for i = 1:numel(gof.coefname)
S.(gof.coefname(i)) = gof.coef(i);
S.(gof.coefname(i)+"_CI") = gof.confint(:,i);
end
end



