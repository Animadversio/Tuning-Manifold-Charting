function paired_stripe_plot(var_cols, varnms, msks, labels,varargin)
Nrow = numel(var_cols{1});
Ncol = numel(var_cols);
varmat = [];
for i = 1:Ncol
    assert(numel(var_cols{i})==Nrow)
    varmat(:,i) = var_cols{i};
end
if isempty(msks), msks={ones(Nrow,1,'logical')}; labels="all";
else
for mi=1:numel(msks) % substitute empty mask
    if isempty(msks{mi}), msks{mi}=ones(Nrow,1,'logical');end
end
end
xjit = 0.08 * randn(Nrow,1);
mskofs = 0.8 / numel(msks);
offsets = ([1:numel(msks)] - 1 - (numel(msks)-1)/2) * mskofs;
figure;hold on;set(gcf,'pos',[1000         357         430         625])
for i=1:Ncol
for mi=1:numel(msks)
    offset = offsets(mi);
    varmean = mean(varmat(msks{mi},i));
    varsem = sem(varmat(msks{mi},i));
    varN = numel(varmat(msks{mi},i));
    scatter(i+offset+xjit(msks{mi}), varmat(msks{mi},i),...
        'DisplayName',compose("%s-%s %.1f+-%.1f(N=%d)",varnms(i),labels(mi),varmean,varsem,varN),varargin{:});
end
end
for mi=1:numel(msks)
plot(offsets(mi)+[1:Ncol]'+xjit(msks{mi})',varmat(msks{mi},:)','color',[0,0,0,0.1],'HandleVisibility','off')
end
xticks(1:Ncol); xticklabels(varnms)
legend('Location','best')
% title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))
end