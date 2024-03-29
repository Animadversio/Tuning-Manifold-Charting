%% Manifold Paper Figure 4. Example Manifold tuning maps within different masks. 
Animal="Both";Set_Path;
addpath analysis\libs\Kent_func

%% Merge the stats for 2 monkeys 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
MapVar_both.(Animal) = MapVarStats;
end
%% Load the population stats
poptabdir = "O:\Manif_Fitting\popstats";
poptab = readtable(fullfile(poptabdir,"Both_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
%% Create the masks and find the drivers 
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = poptab.area == "V1";
V4msk = poptab.area == "V4";
ITmsk = poptab.area == "IT";
drivermsk = poptab.is_driver;
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
for msk = msks % random sample channel x exp to plot 
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


