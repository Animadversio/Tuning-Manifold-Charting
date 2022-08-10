%  Post processing after Fitting Kent function, add relavent info to the csv for better analysis.
%% Transcribe other useful info from EStats to poptab. 
mat_dir = "O:\Mat_Statistics";
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal+"_Evol_stats.mat"), 'EStats')
EStats_both.(Animal) = EStats;
end
%%
poptabdir = "O:\Manif_Fitting\popstats";
% alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat.csv"));
% betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat.csv"));
% poptab = [alfatab_pop;betotab_pop];
alfafittab = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
betofittab = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
poptab = [alfafittab;betofittab];
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_both.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
poptab.is_driver = drivermsk;
%% Validate the is_driver mask correct
assert (sum(poptab.is_driver) == 91 + 20)
assert (sum(poptab.is_driver & poptab.space==1) == 91)
assert (sum(poptab.is_driver & poptab.Animal=="Alfa") == 46)
%% Transcribe other useful info from EStats to poptab. 
poptab.area = arrayfun(@area_map, poptab.chan);
for i = 1:size(poptab,1)
    poptab.imgsize(i) = EStats_both.(poptab.Animal{i})(poptab.Expi(i)).evol.imgsize;
    poptab.imgposX(i) = EStats_both.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(1);
    poptab.imgposY(i) = EStats_both.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(2);
end
%%
writetable(poptab, fullfile(poptabdir,"Both_Exp_all_KentStat_bsl_pole.csv"))
% drivertab = poptab(drivermsk,:);
% writetable(drivertab,fullfile(tabdir,"Both_Exp_Driver_KentStat_bsl_pole.csv"))
