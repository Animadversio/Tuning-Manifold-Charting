function h=violin_plot_masks(tab, varnm, msks, labels, minormsks, minorlabels, varargin)
for mi=1:numel(msks) % substitute empty mask
    if isempty(msks{mi}), msks{mi}=ones(size(tab,1),1,'logical');end
end
for mj=1:numel(minormsks)
    if isempty(minormsks{mj}), minormsks{mj}=ones(size(tab,1),1,'logical');end
end
Yvec = [];
labvec = [];
statstr = "";
lab_seq = string([]);
for mi=1:numel(msks)
    for mj=1:numel(minormsks)
    msk = msks{mi}&minormsks{mj};
    num = sum(msk);
    lab = compose("%s-%s",labels(mi),minorlabels(mj));
    ydata = reshape(tab.(varnm)(msk),[],1);
    statstr = statstr+compose("%s-%s %.3f+-%.3f N=%d\n",labels(mi),minorlabels(mj),mean(ydata),sem(ydata),num);
    Yvec = cat(1, Yvec, ydata);
    labvec = cat(1, labvec, repmat(lab,num,1));
    lab_seq(end+1) = lab;
    end
end
figure;hold on;set(gcf,'pos',[1000         357         430         625])
h=violinplot(Yvec,labvec,varargin{:},'GroupOrder',cellstr(lab_seq));
% legend('Location','best')
ylabel(varnm,'interpreter','none')
annotation('textbox',[.2 .1 .3 .3],'String',statstr,'FitBoxToText','on')

% Nrow = numel(var_cols{1});
% Ncol = numel(var_cols);
% varmat = [];
% for i = 1:Ncol
%     assert(numel(var_cols{i})==Nrow)
%     varmat(:,i) = var_cols{i};
% end
% xjit = 0.08 * randn(Nrow,1);
% mskofs = 0.8 / numel(msks);
% offsets = ([1:numel(msks)] - 1 - (numel(msks)-1)/2) * mskofs;
% figure;hold on;set(gcf,'pos',[1000         357         430         625])
% for i=1:Ncol
% for mi=1:numel(msks)
%     offset = offsets(mi);
%     scatter(i+offset+xjit(msks{mi}), varmat(msks{mi},i),'DisplayName',compose("%s-%s",varnms(i),labels(mi)),varargin{:});
% end
% end
% for mi=1:numel(msks)
% plot(offsets(mi)+[1:Ncol]'+xjit(msks{mi})',varmat(msks{mi},:)','color',[0,0,0,0.1],'HandleVisibility','off')
% end
% xticks(1:Ncol); xticklabels(varnms)
% title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))
end