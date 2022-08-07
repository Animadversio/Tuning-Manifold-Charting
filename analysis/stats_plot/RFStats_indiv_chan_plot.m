function pred_S = RFStats_indiv_chan_plot(S, chan2plot, Xq, Yq)
if nargin < 2
    chan2plot = 1:numel(S.unit.unit_num_arr);
end
if nargin < 3, Xq = -8:0.2:8; end
if nargin < 4, Yq = -8:0.2:8; end
P.plotEachChan = true;
pixperdeg = 40;

ephysFN = S.meta.ephysFN;
uniqsize_deg_all = unique(S.stim.size_deg);
rsp_vec = S.psth.rsp_vec;
bsl_mean = S.psth.bsl_mean;
xy_all = S.stim.xy_all;
size_all = S.stim.size_all;
unit_num_arr = S.unit.unit_num_arr;
pred_S = repmat(struct(),1,numel(S.unit.unit_num_arr));
if P.plotEachChan
% figure;plot(psthbest_mean')
for iCh = chan2plot
X = [xy_all,size_all/pixperdeg];
y = squeeze(rsp_vec(iCh,:))-bsl_mean(iCh);
gprMdl = fitrgp(X,y);
% Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
colN = numel(uniqsize_deg_all);
figure(2);clf;T=tiledlayout(3,colN,'tilespac','compact');
set(2,'pos',[400   60   320*colN   912])
for iSz = 1:numel(uniqsize_deg_all)
sz = uniqsize_deg_all(iSz);
nexttile(T,iSz)
% imagesc(S.stim.xpos{iSz},S.stim.ypos{iSz},S.psth.score_mean{iSz}(:,:,iCh))
imagesc(S.stim.xpos{iSz},S.stim.ypos{iSz},S.psth.act_mean{iSz}(:,:,iCh) - S.psth.bsl_mean(iCh))
axis image;set(gca,'YDir','normal');xlim([min(Xq),max(Xq)]);ylim([min(Yq),max(Yq)])
title(compose("size %.1f deg",sz))
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*sz]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,colN+iSz)
imagesc(Xq,Yq,pred_rfmat);
title(compose("size %.1f deg",sz))
axis image;set(gca,'YDir','normal')
nexttile(T,colN*2+iSz)
shadedErrorBar(1:200,S.psth.psth_best_mean{iSz}(iCh,:),S.psth.psth_best_sem{iSz}(iCh,:))
title(compose("size %.1f deg",sz))
end
title(T,compose("%s Unit %s",ephysFN,S.unit.unit_name_arr(iCh)))
if colN > 1
axs = get(T,'Children');
AlignAxisLimits({axs(1:colN:3*colN)});
AlignAxisCLimits(axs(2:colN:3*colN));
AlignAxisCLimits(axs(3:colN:3*colN));
end
pause
end
end
end