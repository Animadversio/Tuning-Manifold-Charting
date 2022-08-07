function figh = RF_contour_plot(RFStat, style, chanmask)
% Plot RF of population of neurons as one contour plot with area filled or not. 
% chanmask: mask for the channels to be plotted in a plot.
% style: a string "edge" "fill"
if nargin == 1, style="edge"; end
if nargin < 3, 
    chanmask = ones(size(RFStat.unit.chan_num_arr),'logical');
end
area_colormap = containers.Map({'V1','V4','IT'},{[1,0,0],[0,1,0],[0,0,1]});

uniqsize_deg = unique(RFStat.stim.size_deg)';
msk = chanmask;
Fmsk = RFStat.stats.F_P<0.001;
Tmsk = RFStat.stats.T_P<0.001;
% activ_chans = find(Fmsk & Tmsk & msk)';
activ_chans = find(Tmsk & msk)';%Fmsk & 
% activ_chans = find(Fmsk & Tmsk & msk)';
% activ_chans = find(msk)';
Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
Nsize = numel(uniqsize_deg);
figh=figure();clf;T=tiledlayout(1,Nsize,'tilespac','tight','padding','compact');
set(figh,'pos',[200   193   306*Nsize  320   ])
for iCh = activ_chans%targmsk_chans%activ_chans%
chan_num = RFStat.unit.chan_num_arr(iCh);
area = area_map(chan_num);
clr = area_colormap(area);
pixperdeg = 40;
X = [RFStat.stim.xy_all,RFStat.stim.size_all/pixperdeg];
y = RFStat.psth.rsp_vec(iCh,:)-RFStat.psth.bsl_mean(iCh);
gprMdl = fitrgp(X,y);
for iSz = 1:Nsize
sz = uniqsize_deg(iSz);
% imagesc(S.stim.xpos{iSz},S.stim.ypos{iSz},S.psth.score_mean{iSz}(:,:,iCh))
% axis image;set(gca,'YDir','normal');xlim([min(Xq),max(Xq)]);ylim([min(Yq),max(Yq)])
% title(compose("size %.1f deg",sz))
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*sz]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,iSz);hold on
if style == "edge"
contour(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'linecolor',[clr],'linewidth',0.5);
elseif style == "fill"
[C,hcf] = contourf(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'facecolor',[clr]);
allH = allchild(hcf);
valueToHide = 0.5*max(pred_rfmat,[],'all');
patchValues = cell2mat(get(allH,'UserData'));
patchesToHide = patchValues == valueToHide;
set(allH(patchesToHide),'FaceColor',[clr],'FaceAlpha',0.05);
end
end
% pause7
end
for iSz = 1:Nsize
sz = uniqsize_deg(iSz);
nexttile(T,iSz);hold on;title(compose("imsize %.1f deg",sz));axis equal
end
title(T,RFStat.meta.ephysFN)
end