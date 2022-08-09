%% Manif_Tune_Smooth_Popul
%  Code to compute Smoothness statistics for all the units for all exps. 
Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_MapSmooth\popstats";
% load(fullfile(Matdir, "Beto_ManifPopDynamics.mat"))
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
%%
% StatsTab_sum = [];
pool = parpool(5); % distribute the computation since it can take quite a while
%%
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
toc
smthStats_exp = cellfun(@(actmap)smooth_stats(actmap),actmap_col);
anovaStats_exp = cellfun(@anova_cells,actmap_col);
toc
anovaStats_tab = struct2table(anovaStats_exp);
smthStats_tab = struct2table(smthStats_exp);
StatsTab = [Animal_tab, Expi_tab, unitstr_tab, unit_tab, chan_tab, prefchan_tab, space_tab, anovaStats_tab, smthStats_tab];
writetable(StatsTab,fullfile(tabdir,compose("%s_Exp%d_sp%d.csv",Animal,Expi,si)));
% StatsTab_sum = [StatsTab_sum;StatsTab];
% disp(size(StatsTab_sum))
end
end
toc(T00);
%% Merge the stats to a full table.
StatsTab_sum = [];
for Expi = 1:numel(Stats)
for si = 1:MapVarStats(Expi).manif.subsp_n
curtab = readtable(fullfile(tabdir,compose("%s_Exp%d_sp%d.csv",Animal,Expi,si)));
StatsTab_sum = [StatsTab_sum;curtab];
end
end
writetable(StatsTab_sum,fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));

%% Load and Summarize the stats (see the summary for more)
StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
%%
drivermsk = (StatsTab_sum.chan==StatsTab_sum.prefchan);
tunemsk = (StatsTab_sum.F_P<1E-6);
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk));
%%
figure;hold on 
histogram(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk))
histogram(StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk))
title(compose("Comparison of Gap of Dirichlet Energy for Driver\nand Non-Driver channels in all Experiments"))
legend(["Driver", "Non-Drivers"])
%% Utility functions to compute smoothness stats
% SmthStat = smooth_stats(act_col,"sphere");
function SmthStat = smooth_stats(act_col, varargin)
if nargin==1, stat_opt = "sphere"; 
elseif nargin>=2, stat_opt = varargin{1}; end
if strcmp(stat_opt,"sphere"), energy_fun = @sph_map_energy;
elseif strcmp(stat_opt,"euclid"), energy_fun = @map_energy;
else, energy_fun = @sph_map_energy;end
actmap = cellfun(@mean, act_col);
[D2E, TVE, D1E, D2map, GPmap] = energy_fun(actmap);
btrp_Scol = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col); % Resample trial responses
[D2E_btrp, TVE_btrp, D1E_btrp, D2map_btrp, GPmap_btrp] = energy_fun(actmap_btrp);
btrp_Scol = [btrp_Scol;[D2E_btrp, TVE_btrp, D1E_btrp]];
end
shfl_Scol = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_shfl, TVE_shfl, D1E_shfl, D2map_shfl, GPmap_shfl] = energy_fun(actmap_shfl);
shfl_Scol = [shfl_Scol;[D2E_shfl, TVE_shfl, D1E_shfl]];
end
[~,P_D2E,~,tSt_D2E] = ttest2(btrp_Scol(:,1), shfl_Scol(:,1));
Dpr_D2E = computeCohen_d(btrp_Scol(:,1), shfl_Scol(:,1));
[~,P_TVE,~,tSt_TVE] = ttest2(btrp_Scol(:,2), shfl_Scol(:,2));
Dpr_TVE = computeCohen_d(btrp_Scol(:,2), shfl_Scol(:,2));
[~,P_D1E,~,tSt_D1E] = ttest2(btrp_Scol(:,3), shfl_Scol(:,3));
Dpr_D1E = computeCohen_d(btrp_Scol(:,3), shfl_Scol(:,3));
SmthStat.D2E = D2E; 
SmthStat.TVE = TVE; 
SmthStat.D1E = D1E; 
SmthStat.D2map = D2map; 
SmthStat.GPmap = GPmap; 
SmthStat.D2E_btrp_mean = mean(btrp_Scol(:,1));
SmthStat.D2E_shfl_mean = mean(shfl_Scol(:,1));
SmthStat.D2E_btrp_std = std(btrp_Scol(:,1));
SmthStat.D2E_shfl_std = std(shfl_Scol(:,1));
SmthStat.TVE_btrp_mean = mean(btrp_Scol(:,2));
SmthStat.TVE_shfl_mean = mean(shfl_Scol(:,2));
SmthStat.TVE_btrp_std = std(btrp_Scol(:,2));
SmthStat.TVE_shfl_std = std(shfl_Scol(:,2));
SmthStat.D1E_btrp_mean = mean(btrp_Scol(:,3));
SmthStat.D1E_shfl_mean = mean(shfl_Scol(:,3));
SmthStat.D1E_btrp_std = std(btrp_Scol(:,3));
SmthStat.D1E_shfl_std = std(shfl_Scol(:,3));
SmthStat.D2E_t = tSt_D2E.tstat;
SmthStat.D2E_P = P_D2E;
SmthStat.D2E_Dpr = Dpr_D2E;
SmthStat.TVE_t = tSt_TVE.tstat;
SmthStat.TVE_P = P_TVE;
SmthStat.TVE_Dpr = Dpr_TVE;
SmthStat.D1E_t = tSt_D1E.tstat;
SmthStat.D1E_P = P_D1E;
SmthStat.D1E_Dpr = Dpr_D1E;
end

function [laplsEng,TVEng,dirEng,laplsMap,gradPowMap] = map_energy(actmap)
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
% laplsMap = conv2(actmap, kerd2);
laplsMap = conv2(padarray(actmap,[1,1],'replicate'), kerd2,'valid');
[actGx,actGy] = gradient(actmap);
gradPowMap = actGx.^2 + actGy.^2;
dirEng = sum(gradPowMap,'all');
TVEng = sum(sqrt(gradPowMap),'all');
laplsEng = sum(abs(laplsMap),'all');
end
function [laplsEng_sp, TVEng_sp, dirEng_sp, laplsMap_sp, gradNormMap_sp, detJac_map] = sph_map_energy(actmap)
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); 
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
% detJac = abs(cosd(PHI));
Ginv11 = 1./cosd(PHI).^2;Ginv11(:,1)=0;Ginv11(:,end)=0; % These components explode, set to 0. 
Ginv22 = ones(size(PHI));
% Sample points for the gradients are far away.
detJac_map = cosd([mean(PHI(:,[1,2]),2),PHI(:,2:end-1),mean(PHI(:,[end-1,end]),2)]);
actmap_sp = actmap;
actmap_sp(:,1)=mean(actmap(:,1));
actmap_sp(:,end)=mean(actmap(:,end));
laplsMap_sp = conv2(padarray(actmap_sp,[1,1],'replicate'), kerd2,'valid');%del2(actmap_sp);%
[actGx_sp,actGy_sp] = gradient(actmap_sp);
gradNormMap_sp = actGx_sp.^2 .* Ginv22 + actGy_sp.^2 .*Ginv11;
dirEng_sp = sum(gradNormMap_sp.*detJac_map,'all');
TVEng_sp = sum(sqrt(gradNormMap_sp).*detJac_map,'all');
laplsEng_sp = sum(abs(laplsMap_sp).*detJac_map,'all');
end