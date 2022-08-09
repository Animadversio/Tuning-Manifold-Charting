%% Manifold paper Supplementary figure 06/07
%% Plot montage of numerous tuning maps (of driving units) eps. grouped by area or animal
Animal="Both";Set_Path;
ang_step = 18;
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
poptabdir = "O:\Manif_Fitting\popstats";
sumdir = "O:\Manif_Fitting\summary";
mtgdir = "O:\Manif_Fitting\montage";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');%Old default behavior
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
poptab = [alfatab_pop;betotab_pop];
% load in stats
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats')
MapVarStats_all.(Animal) = MapVarStats;
Stats_all.(Animal) = Stats;
EStats_all.(Animal) = EStats;
end

%% Prep the poptab and masks 
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = (poptab.chan<=48 & poptab.chan>=33);
V4msk = (poptab.chan>48);
ITmsk = (poptab.chan<33);
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
prefchmsk = poptab.chan==poptab.prefchan;
Fsigmsk = poptab.F_P<1E-3;

%% Collect tuning maps of drivers. 
MapsCol = {};
idxlist = find(drivermsk);% & Alfamsk);
selectTab = poptab(idxlist,:);
selectTab.area = arrayfun(@area_map,selectTab.chan);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Anim = Tab.Animal{1}; Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
iCh = find((MapVarStats_all.(Anim)(Expi).units.spikeID==ci & MapVarStats_all.(Anim)(Expi).units.unit_num_arr==ui));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats_all.(Anim)(Expi).manif.act_col{spi},'uni',1); % 51:200 activity
bslvec = cell2mat(reshape(cellfun(@(A)A(iCh,:),MapVarStats_all.(Anim)(Expi).manif.bsl_col{spi},'uni',0),1,[]));
bslmean = mean(bslvec);
scoremap = actmap_mean - bslmean;
MapsCol{i} = scoremap; % baseline subtracted score map. 
end

%% Visualize maps separated by masks
TuningMapMontage(MapsCol, selectTab, "flat", selectTab.area=="V1", [5,-1], "Alfa_V1")
saveallform(mtgdir, "MapsMtg_Alfa_V1")
TuningMapMontage(MapsCol, selectTab, "flat", selectTab.area=="V4", [5,-1], "Alfa_V4")
saveallform(mtgdir, "MapsMtg_Alfa_V4")
TuningMapMontage(MapsCol, selectTab, "flat", selectTab.area=="IT", [5,-1], "Alfa_IT")
saveallform(mtgdir, "MapsMtg_Alfa_IT")

%% Plot raw maps flat as heatmap 
valmsk = selectTab.space==1&selectTab.F_P<0.001;
msk_col = {selectTab.area=="V1"&valmsk,...
           selectTab.area=="V4"&valmsk,...
           selectTab.area=="IT"&valmsk};
label_col = {"Both_Fsig_V1_bsl","Both_Fsig_V4_bsl","Both_Fsig_IT_bsl"};
for i = 1:numel(msk_col)
msk = msk_col{i}; label = label_col{i};
TuningMapMontage(MapsCol, selectTab, "flat", msk, [5,-1], label)
saveallform(mtgdir, "MapsMtg_"+label)
end

%% Plot Kent fit on the sphere and montage. 
valmsk = selectTab.space==1&selectTab.F_P<0.001;
msk_col = {selectTab.area=="V1"&valmsk,...
           selectTab.area=="V4"&valmsk,...
           selectTab.area=="IT"&valmsk};
label_col = {"Both_Fsig_V1_bsl_fit","Both_Fsig_V4_bsl_fit","Both_Fsig_IT_bsl_fit"};
for i = 1:numel(msk_col)
msk = msk_col{i}; label = label_col{i};
% TuningMapMontage(MapsCol, selectTab, "fit_flat", msk, [5,-1], label)
% saveallform(mtgdir, "MapsMtg_"+label+"_flat")
TuningMapMontage(MapsCol, selectTab, "fit_sphere", msk, [-1,7], label)
saveallform(mtgdir, "MapsMtg_"+label+"_sphere")
end

function TuningMapMontage(MapsCol, selectTab, mode, mask, tile_size, suptit)
% Unified plotting function for plotting raw or fit maps as a montage. 
% 
%  MapsCol: cell array of raw maps each `MapsCol{mapi}` will be an 11x11 tuning map. 
%           length of it is N
%  selectTab: a length N table, containing the Kent fitting information for each channel. 
%  mode: string. it could be "flat", "sphere", "fit_flat", "fit_sphere"
%        if it contains `fit` then the Kent fitting is plotted. else, the raw data
%        if it contains `flat` then use imagesc, else plot it on 3d sphere surface. 
%  mask:  a length `N` boolean array showing which map to plot. 
%  tile_size: 1x2 array. [row, col], if one of them is `-1` then it will be determined from the other one auto. 
%  suptit: string. supertitle string / label for the mask.
% 
if nargin == 2, mode = "flat"; end
if nargin < 4, mask = ones(numel(MapsCol), 1, 'logical');end
if nargin < 5, tile_size = [5,10]; end 
if nargin < 6, suptit = "";end
ang_step = 18;
phi_arr = -90:18:90; theta_arr = -90:18:90;
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 

assert(numel(mask)==numel(MapsCol) && numel(mask)==size(selectTab,1))

idxlist = find(mask);
if tile_size(1)==-1 && tile_size(2)~=-1  % fixed column number, flexible row number
    tile_size(1) = ceil(numel(idxlist)/tile_size(2));
elseif tile_size(2)==-1  && tile_size(1)~=-1  % fixed row number, flexible column number
    tile_size(2) = ceil(numel(idxlist)/tile_size(1));
end
figure(26);clf;set(26,'pos',[292, 66, 190*tile_size(2), 150*tile_size(1)+45])
T=tiledlayout(tile_size(1),tile_size(2),'TileSpacing','tight','padding','none');
for i = 1:numel(idxlist)
mapi = idxlist(i);
ax = nexttile(T,i);
Tab = selectTab(mapi,:);
if strcmp(mode,"flat")
imagesc(ax,-90:18:90, -90:18:90, MapsCol{mapi}); colorbar(); axis image;axis off
elseif strcmp(mode,"sphere")
sphere_plot(ax, theta_grid, phi_grid, MapsCol{mapi});axis off
elseif contains(mode,"fit")
    fitval = KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, ...
        reshape(theta_grid,[],1)/180*pi, reshape(phi_grid,[],1)/180*pi) + Tab.bsl;
    fitval = reshape(fitval,size(theta_grid));
    if strcmp(mode,"fit_flat")
    imagesc(ax,-90:18:90, -90:18:90, fitval); colorbar(); axis image;axis off
    elseif strcmp(mode,"fit_sphere")
    sphere_plot(ax, theta_grid, phi_grid, fitval);axis off
    end
end
title(compose("%s Exp%02d %s",Tab.Animal{1}(1), Tab.Expi, Tab.unitstr{1}),'FontSize',10)
% title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
% 	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end
title(T,compose("Tuning Map Montage %s",suptit),'interpreter','none')
end

function area = area_map(chan)
if chan<=48 & chan>=33
area = "V1";
elseif chan>48
area = "V4";
elseif chan<33
area = "IT";
end
end
