function Manif_Map_show_fun(MapVarStats, Animal, Expi, chan2plot, figh)
if nargin < 5
    figh = 3;
end
figure(figh); clf; set(figh,'pos',[ 805         197        1559         781]); % all manifold images montaged
set(0,'CurrentFigure',figh); clf;T=tiledlayout('flow','tilesp','compact','padd','compact');
% Put Manifold maps into a flow tiled layout
si = 1;
% T=tiledlayout(1,2,'tilesp','compact','padd','compact');
pref_chan = MapVarStats(Expi).units.pref_chan;
tic
for iCh = reshape(chan2plot,1,[])
% unitstr = Stats(Expi).units.unit_name_arr{iCh};
unitstr = MapVarStats(Expi).units.unit_name_arr{iCh};
chan_label_str = sprintf("%s Exp%d (Driver %d) Channel %s", Animal, Expi, pref_chan, unitstr);
actvec = cell2mat(reshape(cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.act_col{si},'uni',0),1,[]));
actmap = cellfun(@(A)mean(A(iCh,:),'all'),MapVarStats(Expi).manif.act_col{si},'uni',1);
actcol = cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.act_col{si},'uni',0);
bslvec = cell2mat(reshape(cellfun(@(A)A(iCh,:),MapVarStats(Expi).manif.bsl_col{si},'uni',0),1,[]));
bslmean = mean(bslvec);
bslstd = std(bslvec);
scoremap = actmap - bslmean;
[~,P,CI,TST] = ttest(actvec, bslvec);
tval = TST.tstat; pval = P;
anovaS = anova_cells(actcol(:));
stat_str = compose(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)'],CI(1), CI(2), P, anovaS.F, anovaS.F_P);

ax1 = nexttile(T);%1);%subplot(1,2,1);
imagesc(-90:18:90, -90:18:90, scoremap) 
ylabel("PC 2 degree");xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
axis image;colorbar%shading flat;
toc
end
end