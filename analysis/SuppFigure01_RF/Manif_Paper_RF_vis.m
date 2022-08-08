
area_colormap = containers.Map({'V1','V4','IT'},{[1,0,0],[0,1,0],[0,0,1]});
RFCol.Alfa = load(fullfile(matdir,"Alfa"+"_Manif_RFstats.mat"));
RFCol.Beto = load(fullfile(matdir,"Beto"+"_Manif_RFstats.mat"));
%%
% Animal="Alfa"; Mapi = 42; 
Animal="Beto"; Mapi = 15; 
RFStats = RFCol.(Animal).RFStats;
uniqpos = RFStats(Mapi).stim.uniqpos;
xpos = RFStats(Mapi).stim.xpos;
ypos = RFStats(Mapi).stim.ypos;
uniqsize_deg = unique(RFStats(Mapi).stim.imgsize_deg)';
fprintf("%s\nx position from %.1f to %.1f (N=%d)\n y position from %.1f to %.1f (N=%d)\nimage size %.1f\n",...
        RFStats(Mapi).meta.ephysFN,xpos(1),xpos(end),numel(xpos),ypos(1),ypos(end),numel(ypos),uniqsize_deg)
Nsize = RFStats(Mapi).stim.nSize; 
%%
x_vect = uniqpos(:,1); % reshape(posgrid(:,:,1),1,[]);
y_vect = uniqpos(:,2); % reshape(posgrid(:,:,2),1,[]);
nPos = size(uniqpos,1);
posgrid = reshape([uniqpos, (1:nPos)'],[length(xpos),length(ypos),3]);
posvec = uniqpos;
sizevec = ones(size(uniqpos,1),1)*uniqsize_deg;

Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
figure(2);clf;T=tiledlayout(1,Nsize,'tilespac','tight','padding','compact');
set(2,'pos',[1000   193   306*Nsize  320   ])
activ_msk = (RFStats(Mapi).stats.anovaP<0.001) & (RFStats(Mapi).stats.ttestP<0.01);
activ_list = find(activ_msk)';

for iCh = activ_list % 1:numel(RFStats(Mapi).meta.spikeID)
chan_num = RFStats(Mapi).meta.spikeID(iCh);
area = area_map(chan_num);
clr = area_colormap(area);
actmap = cellfun(@(P)mean(P(iCh,51:200),'all'),RFStats(Mapi).psth.psth_mean);
bslmap = cellfun(@(P)mean(P(iCh,1:50),'all'),RFStats(Mapi).psth.psth_mean);

% X = [posvec, sizevec];
X = [posvec];
y = reshape(actmap - mean(bslmap,'all'),[],1);
gprMdl = fitrgp(X,y);
% pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*uniqsize_deg]);
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1)]);%,ones(numel(XX),1)*uniqsize_deg]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,1);hold on
% contour(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'linecolor',[clr],'linewidth',0.5);
[C,hcf] = contourf(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'facecolor',[clr]);
if mean(pred_rfmat>0.5*max(pred_rfmat,[],'all'))>0.6,keyboard;end
allH = allchild(hcf);
valueToHide = 0.5*max(pred_rfmat,[],'all');
patchValues = cell2mat(get(allH,'UserData'));
patchesToHide = patchValues == valueToHide;
set(allH(patchesToHide),'FaceColor',[clr],'FaceAlpha',0.05);
end
%%
saveallform("O:\Manif_RF",Animal+"_RF_Demo",2)
%%
Mapis = [41,42,43];
Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
figure(2);clf;T=tiledlayout(1,Nsize,'tilespac','tight','padding','compact');
set(2,'pos',[1000   193   306*Nsize  320   ])

posvecs = [];
sizevecs = [];
for Mapi = Mapis
uniqpos = RFStats(Mapi).stim.uniqpos;
xpos = RFStats(Mapi).stim.xpos;
ypos = RFStats(Mapi).stim.ypos;
uniqsize_deg = unique(RFStats(Mapi).stim.imgsize_deg)';
posvecs = [posvecs; uniqpos];
sizevecs = [sizevecs; ones(size(uniqpos,1),1)*uniqsize_deg];
end
X = [posvecs,sizevecs];
%%
activ_msk = (RFStats(Mapi).stats.anovaP<0.001) & (RFStats(Mapi).stats.ttestP<0.01);
activ_list = find(activ_msk)';
for iCh = activ_list
chan_num = RFStats(Mapi).meta.spikeID(iCh);
area = area_map(chan_num);
clr = area_colormap(area);

bslvec = [];
actvec = [];
for Mapi = Mapis
actmap = cellfun(@(P)mean(P(iCh,51:200),'all'),RFStats(Mapi).psth.psth_mean);
bslmap = cellfun(@(P)mean(P(iCh,1:50),'all'),RFStats(Mapi).psth.psth_mean);
bslvec = [bslvec; reshape(bslmap,[],1)];
actvec = [actvec; reshape(actmap,[],1)];
end
y = actvec - mean(bslvec,'all');

gprMdl = fitrgp(X,y);
% pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*uniqsize_deg]);
sizelist = [1,3,6];
for iSz = 1:3
sz = sizelist(iSz);
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*sz]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,iSz);hold on
% contour(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'linecolor',[clr],'linewidth',0.5);
[C,hcf] = contourf(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'facecolor',[clr]);
if mean(pred_rfmat>0.5*max(pred_rfmat,[],'all'))>0.6,keyboard;end
allH = allchild(hcf);
valueToHide = 0.5*max(pred_rfmat,[],'all');
patchValues = cell2mat(get(allH,'UserData'));
patchesToHide = patchValues == valueToHide;
set(allH(patchesToHide),'FaceColor',[clr],'FaceAlpha',0.05);
end
end
% get x,y,size
% formulate response
% fit model 
% predict 
% plot contour
%%
saveallform("O:\Manif_RF",Animal+"_RF_Demo_cmb",2)

%% Newer version of API
Animal="Beto";Set_Path;
%%
ExpRecord(find(contains(ExpRecord.expControlFN,"rf")),:)
%%
% expidx = find(contains(ExpRecord.ephysFN,'Beto-23042020-001'));
ephysFN = ["Beto-20042020-001","Beto-23042020-001","Beto-27042020-001"];
expidx = find(contains(ExpRecord.ephysFN,ephysFN));%
[meta_new,rasters_new,lfps_new,Trials_new] = loadExperiments(expidx,"Beto");
S_col = RF_Calc_Stats_fun(meta_new, rasters_new, Trials_new);
%%
for Si = 1:numel(S_col)
h=RF_contour_plot(S_col(Si));
% saveallform("O:\RFstats\RF_contours",S_col(Si).meta.ephysFN+"_RF_merge",h);
end
%%
ephysFN = "Beto-05122019-001";
expidx = find(contains(ExpRecord.ephysFN,ephysFN));%
[meta_new,rasters_new,lfps_new,Trials_new] = loadExperiments(expidx,"Beto");
S_col = RF_Calc_Stats_fun(meta_new, rasters_new, Trials_new);
%%
S_col = RF_Calc_Stats_fun({meta}, {rasters}, {Trials});
h=RF_contour_plot(S_col(Si));
saveallform("O:\RFstats\RF_contours",S_col(Si).meta.ephysFN+"_RF_merge",h);
%%
h=RF_contour_plot(S_col(Si));
saveallform("O:\RFstats\RF_contours",S_col(Si).meta.ephysFN+"_RF_merge",h);
%%
h=RF_contour_plot(S_col(Si),'fill');
saveallform("O:\RFstats\RF_contours",S_col(Si).meta.ephysFN+"_RF_merge_filled",h);