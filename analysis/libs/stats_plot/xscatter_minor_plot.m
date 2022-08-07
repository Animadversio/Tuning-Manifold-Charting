function h = xscatter_minor_plot(tab, statname1, statname2, masks, labels, ...
        minormsks, minorlabels, titstr, savestr, sumdir, varargin)
% Scatter 2 variables in a table. 
% with special treatments for theta, phi and others. 
marker = 'o*x^v';
Corder = colororder;
clrs = Corder([1,3,5,7],:);
h = figure;clf;hold on;set(h,'pos',[1106         327         560         526]);
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    statcol1 = tab.(statname1)(M);
    statcol2 = tab.(statname2)(M);
    mean_1 = mean(statcol1);
    sem_2 = sem(statcol1);
    mean_2 = mean(statcol2);
    sem_2 = sem(statcol2);
    N_i = sum(M);
    legstr = compose("%s %.3f(%.3f) %.3f(%.3f) (%d)",labels(i)+minorlabels(j),mean_1,sem_2,mean_2,sem_2,N_i);
    scatter(statcol1, statcol2,'Marker',marker(j),'MarkerEdgeColor',clrs(i,:),'DisplayName',legstr,varargin{:})
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
title_str = compose("Scatter of %s - %s for\n %s channels %s",...
    statname1,statname2,join(labels),titstr);
% FStat = anova_cells(statscol);
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
% if any(strcmp(statname1,["phi","theta"])),
% for i = 1:numel(masks)
% title_str = title_str +compose("\n%s",labels(i));%+
% [~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}));
% title_str = title_str + compose(" %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname1,CI(1),CI(2),TST.tstat,P);
% [~,P,CI,TST] = ttest(StatsTab_sum.(statname2)(masks{i}));
% title_str = title_str + compose(" %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname2,CI(1),CI(2),TST.tstat,P);
% end
% end
% if strcmp(statname1,"kappa") && strcmp(statname2,"beta") % if kappa-beta then draw y=x/2 line
%     YLIM = ylim(); XLIM=xlim();
%     LB = max(2*YLIM(1),XLIM(1));
%     UB = min(2*YLIM(2),XLIM(2));
%     plot([LB,UB],[LB,UB]/2,'k-.')
%     for i = 1:numel(masks)
%     title_str = title_str +compose("\n%s",labels(i));
% %     [~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i}) - 2);
%     [~,P,CI,TST] = ttest(log(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i})) - log(2));
%     title_str = title_str + compose("%s / %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname1,statname2,2*exp(CI(1)),2*exp(CI(2)),TST.tstat,P);
%     end
% end
% Ns = cellfun(@sum,masks); % Num of points within each musk
xlabel(statname1,'interpreter','none');ylabel(statname2,'interpreter','none'); %
title(title_str,'interpreter','none')
legend('Location','best')
% legend(compose("%s (%d)",labels',Ns')) % Num of points marked on label
if islogical(savestr) && (~savestr), return
else
saveallform(sumdir, compose("%s_%s_xscat_%s", statname1, statname2, savestr), h)
end
end