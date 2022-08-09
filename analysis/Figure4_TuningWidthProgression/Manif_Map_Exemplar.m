%% Manifold Paper Figure 4. Example Manifold tuning maps within different masks. 
Animal="Both";Set_Path;
addpath D:\Github\Fit_Spherical_Tuning
addpath e:\Github_Projects\Fit_Spherical_Tuning
% tabdir = "O:\Manif_Fitting\Kent_summary";
% alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
% betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
% preftab = [];
% Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
% preftab = [preftab; Animal_tab, alfatab];
% Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
% preftab = [preftab; Animal_tab, betotab];

%% Merge the stats for 2 monkeys 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"), 'Stats')
load(fullfile(mat_dir, Animal+"_Evol_stats.mat"), 'EStats')
MapVar_both.(Animal) = MapVarStats;
Stats_both.(Animal) = Stats;
EStats_both.(Animal) = EStats;
end
%% Load the population stats
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat.csv"));
poptab = [alfatab_pop;betotab_pop];
%% Create the masks and find the drivers 
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = (poptab.chan<=48 & poptab.chan>=33);
V4msk = (poptab.chan>48);
ITmsk = (poptab.chan<33);
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_both.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
prefchmsk = poptab.chan==poptab.prefchan;
Fsigmsk = poptab.F_P<1E-3;

%% Main part
sumdir = "O:\Manif_Fitting\summary";
%% Plot random samples from several masks
h = figure(1);
pubmsk = Fsigmsk&drivermsk&poptab.R2>0.6;
msks = {pubmsk&V1msk&Alfamsk, pubmsk&V4msk&Alfamsk, pubmsk&ITmsk&Alfamsk;
        pubmsk&V1msk&Betomsk, pubmsk&V4msk&Betomsk, pubmsk&ITmsk&Betomsk};
plot_samples(h, msks, MapVar_both, poptab)
RND = randi(1000,1,1);
%%
saveallform(sumdir, compose("TuneMapExample%03d.png",RND), h)


%% Explicit no function version
Fsigmsk = poptab.F_P<1E-3;
msks = {Fsigmsk&drivermsk&V1msk, Fsigmsk&drivermsk&V4msk, Fsigmsk&drivermsk&ITmsk};
idxlist = [];
for msk = msks
idx = randsample(find(msk{1}), 1);
idxlist = [idxlist, idx];
end
figure(1);
T=tiledlayout(size(msks,1),size(msks,2),'TileSpacing','compact','padding','compact');
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Animal = Tab.Animal{1}; Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan; unitstr = Tab.unitstr{1};
iCh = find((MapVar_both.(Animal)(Expi).units.spikeID==ci & MapVar_both.(Animal)(Expi).units.unit_num_arr==ui));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2),MapVar_both.(Animal)(Expi).manif.act_col{spi},'uni',1);

ax = nexttile(T,i);
imagesc(-90:18:90,-90:18:90,actmap_mean);axis image
% sphere_plot(ax,theta_grid,phi_grid,actmap_mean);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
     Animal, Expi, unitstr, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end

function plot_samples(h, msks, MapVar_both, poptab)
% Plot manifold experiment on as imagesc in a tiled layout. 
% Tiles will be arranged in the same size as msks. 
% The Manifold experiment will be randomly sampled from the exp/chan under each mask. 
% 
idxlist = [];
for mi = 1:size(msks,1)
for mj = 1:size(msks,2)
idx = randsample(find(msks{mi,mj}), 1);
idxlist = [idxlist, idx];
end
end
figure(h);
T=tiledlayout(size(msks,1),size(msks,2),'TileSpacing','compact','padding','compact');
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Animal = Tab.Animal{1}; Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan; unitstr = Tab.unitstr{1};
iCh = find((MapVar_both.(Animal)(Expi).units.spikeID==ci & MapVar_both.(Animal)(Expi).units.unit_num_arr==ui));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2),MapVar_both.(Animal)(Expi).manif.act_col{spi},'uni',1);

ax = nexttile(T,i);
imagesc(-90:18:90,-90:18:90,actmap_mean);axis image;colorbar() % simple 2D plot 
% sphere_plot(ax,theta_grid,phi_grid,actmap_mean);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
	 Animal, Expi, unitstr, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end
end


