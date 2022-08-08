%% Compute the nonparametric statistics for tuning maps of all exp all channels pop. 
%  Similar content to Manif_NonParametrics_Pop.m 
mat_dir = "O:\Mat_Statistics";
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];

%% Compute and curate non parametric stats about single channel tuning map
tic
S_col = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
for Expi = 1:numel(Stats)
for spi = 1:MapVarStats(Expi).manif.subsp_n
exptab = poptab(poptab.Animal==Animal & poptab.Expi==Expi & poptab.space==spi,:);
assert(size(exptab,1)==numel(MapVarStats(Expi).units.spikeID))
for iCh = 1:numel(MapVarStats(Expi).units.spikeID)
S = struct();
S.Animal = Animal;
S.Expi = Expi;
S.chan = MapVarStats(Expi).units.spikeID(iCh);
S.unitnum = MapVarStats(Expi).units.unit_num_arr(iCh);
S.prefchan = MapVarStats(Expi).units.pref_chan;
S.unitstr = MapVarStats(Expi).units.unit_name_arr(iCh);
S.area = area_map(S.chan);
S.iCh = iCh;
S.F = exptab.F(iCh);
S.F_P = exptab.F_P(iCh);
S.space = spi;
% for nm = ["Animal", "Expi", "unitstr", "unitnum", "chan", "prefchan", "space","F", "F_P"] % copy stats.
% S.(nm) = Tab.(nm);
% end
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2), MapVarStats(Expi).manif.act_col{spi}, 'uni',1);
bslmat = cell2mat(reshape(MapVarStats(Expi).manif.bsl_col{spi},1,[]));
bslmean = mean(bslmat(iCh,:));
% bslstd = std(bslmat(iCh,:));
% Major workhorse here 
S = Manif_Map_NonParam_stat_fun(actmap_mean, bslmean, S);
S_col = [S_col;S];
end
end
toc
end
end % 74 sec for computation
%%
nonpardir = "O:\Manif_NonParam\summary";
NonParTab = struct2table(S_col);
writetable(NonParTab,fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"))
% Check `fullfile(nonpardir,"Both"+"_Pop_NonParamStat.csv"))`


function S = Manif_Map_NonParam_stat_fun(actmap_mean, bslmean, S)
if nargin==2, S = struct(); end
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
XX = cosd(theta_grid).* cosd(phi_grid);
YY = sind(theta_grid) .* cosd(phi_grid);
ZZ = sind(phi_grid);
%  Integration Weight Matrix: 
%  Integrate the area of the patch that each node governs
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9, -90) /180 *pi;
theta2_grid = min(theta_grid + 9,  90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);

[maxAct, maxIdx] = max(actmap_mean,[],'all','linear');
[minAct, ~] = min(actmap_mean,[],'all','linear');
S.maxAct = maxAct;
S.minAct = minAct;
S.Act_range = maxAct - minAct;
S.bslAct = bslmean;
thresh = 0.9 * maxAct; 
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_90 = CoMtheta;
S.CoMphi_90 = CoMphi;
S.CoMrho_90 = CoMrho;
S.Mtheta_90 = Mtheta;
S.Mphi_90 = Mphi;
S.Mrho_90 = Mrho;
% Core of the loop, use the Spherical Integration to get the activation maps. 
integ = SphIntegration(actmap_mean,0,theta_grid,phi_grid,Wgrid);
S.AUS = integ; % Area under tuning map
S.normAUS = integ / maxAct; % Normalized area under tuning map.
integ = SphIntegration(actmap_mean,bslmean,theta_grid,phi_grid,Wgrid);
S.AUS_bsl = integ;
S.normAUS_bsl = integ / (maxAct - bslmean);
end

function [CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap, thresh, XX, YY, ZZ, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
msk = actmap > thresh;
Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
Actvec = actmap(msk);
CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
% [TH,PHI,R] = cart2sph(CoMvec(),CoMvec(),CoMvec())
[CoMtheta, CoMphi, CoMrho] = cart2sph(CoMvec(1),CoMvec(2),CoMvec(3));
[Mtheta, Mphi, Mrho] = cart2sph(Mvec(1),Mvec(2),Mvec(3));
end

function [Integ,meanInteg] = SphIntegration(actmap, baseline, theta_grid, phi_grid, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
mask = (actmap - baseline) > 0;
Integ = sum(mask .* (actmap - baseline) .* Weight,'all');
meanInteg = sum(mask .* (actmap - baseline) .* Weight,'all') / sum(Weight,'all');
% Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
% Actvec = actmap(msk);
% CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
% Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
end