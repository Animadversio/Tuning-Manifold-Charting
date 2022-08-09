%% Manifold Bayes Decoding
%  Final version of bayes decoding 
Set_Path;
mat_dir = "O:\Mat_Statistics";
decodedir = "O:\Manif_PopDecode";
mkdir(decodedir);
%% Set up coordinates and geometry for Manifold
global L2mat angmat cosmat
[PHI, THETA] = meshgrid(-90:18:90.1,-90:18:90.1);
XX = cosd(PHI).*cosd(THETA);
YY = cosd(PHI).*sind(THETA);
ZZ = sind(PHI);
targ_coord = [reshape(XX,[],1),reshape(YY,[],1),reshape(ZZ,[],1)]; % 121 by 3
L2mat = squareform(pdist(targ_coord,'eucl'));
cosmat = squareform(pdist(targ_coord,'cosine'));
angmat = acos(1 - cosmat);

%% Main loop, going through every experiment, perform decoding on it.
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"), 'MapVarStats')
%%
diary(fullfile(decodedir,Animal+"_NaiveBayes_fit.log"))
tic
for Expi = 1:numel(MapVarStats)
si = 1;
prefchanlab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
fprintf("\nProcess Exp %d prefchan %s\n",Expi,prefchanlab)
nCh = numel(MapVarStats(Expi).units.spikeID);
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
% creat masks
FStats = cellfun(@anova_cells,actmap_col); 
Fmsk = struct2table(FStats).F_P < 0.001; 
n_features = sum(Fmsk);
spikeID = MapVarStats(Expi).units.spikeID;
V1msk = spikeID <=48 & spikeID>=33;
V4msk = spikeID <=64 & spikeID>=49;
ITmsk = spikeID <=32 & spikeID>=1;
fprintf("Informative Chan num %d\n",n_features)
actmean = cellfun(@(A)mean(A,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actmean = cell2mat(reshape(actmean,1,[]));
actstd = cellfun(@(A)std(A,1,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actstd = cell2mat(reshape(actstd,1,[]));
actmeanmat = actmean';
actstdmat = actstd';
actmean_sel = actmeanmat(:,Fmsk);
actstd_sel = actstdmat(:,Fmsk);
msklabels = ["All","Fsignif","non Fsignif","V1", "V4", "IT","V1 Fsig", "V4 Fsig", "IT Fsig"];
[DecCol,Dectabs,Dec_summary] = BayesDecwMasks(actmeanmat,actstdmat,targ_coord,...
    {[],Fmsk,~Fmsk,V1msk,V4msk,ITmsk,V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk},...
    ["All","Fsignif","non Fsignif","V1", "V4", "IT","V1 Fsig", "V4 Fsig", "IT Fsig"]);
fprintf("----Shuffled target decoding:------\n")
targ_shfl = targ_coord(randperm(121),:);
[DecCol_shfl,Dectabs_shfl,Dec_summary_shfl] = BayesDecwMasks(actmeanmat,actstdmat,targ_shfl,...
    {[],Fmsk,~Fmsk,V1msk,V4msk,ITmsk,V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk},...
    ["All","Fsignif","non Fsignif","V1", "V4", "IT","V1 Fsig", "V4 Fsig", "IT Fsig"]);
% Single trial response decoding accuracy...???
[sg_actvec, sg_targidx] = get_singletrial_data(MapVarStats(Expi).manif.act_col{si});
sg_targ_coord = targ_coord(sg_targidx, :);
[DecCol_sg,Dectabs_sg,Dec_summary_sg] = BayesDecwMasks(actmeanmat,actstdmat,targ_coord,...
    {[],Fmsk,~Fmsk,V1msk,V4msk,ITmsk,V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk},...
    ["All","Fsignif","non Fsignif","V1", "V4", "IT","V1 Fsig", "V4 Fsig", "IT Fsig"],...
    sg_actvec,sg_targidx);
save(fullfile(decodedir,compose(Animal+"_Manif_Exp%02d_NaiveBayesdecode.mat",Expi)), ...
    "msklabels", "DecCol", "Dectabs", "Dec_summary", "DecCol_shfl", "Dectabs_shfl", "Dec_summary_shfl",  ...
    "DecCol_sg", "Dectabs_sg", "Dec_summary_sg", "FStats", "Fmsk", "actmeanmat", "actstdmat", "targ_coord", "sg_actvec", "sg_targidx")
toc
end
diary off;
end


%% Collect exp into a structure array. BayDecStats
BayDecStats = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
for Expi = 1:numel(Stats)
S = load(fullfile(decodedir,compose(Animal+"_Manif_Exp%02d_NaiveBayesdecode.mat",Expi)), ...
    "msklabels", "DecCol", "Dectabs", "Dec_summary", "DecCol_shfl", "Dectabs_shfl", "Dec_summary_shfl", ...
    "DecCol_sg", "Dectabs_sg", "Dec_summary_sg", "FStats", "Fmsk", "actmeanmat", "actstdmat", "targ_coord", "sg_actvec", "sg_targidx");
S.Animal = Animal;
S.Expi = Expi;
S.imgsize = EStats(Expi).evol.imgsize;
S.imgpos = EStats(Expi).evol.imgpos;
S.prefchan = EStats(Expi).evol.pref_chan;
S.prefunit = EStats(Expi).evol.unit_in_pref_chan;
S.chan2dec = arrayfun(@(D)sum(D.mask),S.Dec_summary);
BayDecStats = [BayDecStats;S];
end
end


%% Summary Decoding result
sumdir = fullfile(decodedir,'summary');
diary(fullfile(sumdir,"Shuffle_summary.log"))
fprintf("Decoding in noise free scenario:\n")
% F signif channel 
msklabels = BayDecStats(1).msklabels; 
for mski = 1:numel(msklabels)
chanN_arr = arrayfun(@(B)B.chan2dec(mski),BayDecStats);
Eang_arr = arrayfun(@(B)B.Dec_summary(mski).EangD_m,BayDecStats);
Eang_arr_shfl = arrayfun(@(B)B.Dec_summary_shfl(mski).EangD_m,BayDecStats);
[~,P,CI,TST] = ttest(Eang_arr,Eang_arr_shfl);
fprintf("%s (N=%.1f+-%.1f) ",msklabels(mski),mean(chanN_arr),std(chanN_arr))
fprintf("Angle Error Exp (%.3f+-%.3f) - Image Shuffled(%.3f+-%.3f): P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",...
    mean(Eang_arr),sem(Eang_arr),mean(Eang_arr_shfl),sem(Eang_arr_shfl),P,TST.tstat,TST.df,CI(1),CI(2));
end
%
msklabels = BayDecStats(1).msklabels; 
for mski = 1:numel(msklabels)
chanN_arr = arrayfun(@(B)B.chan2dec(mski),BayDecStats);
Eang_arr = arrayfun(@(B)B.Dec_summary(mski).EL2D_m,BayDecStats);
Eang_arr_shfl = arrayfun(@(B)B.Dec_summary_shfl(mski).EL2D_m,BayDecStats);
[~,P,CI,TST] = ttest(Eang_arr,Eang_arr_shfl);
fprintf("%s (N=%.1f+-%.1f) ",msklabels(mski),mean(chanN_arr),std(chanN_arr))
fprintf("L2 Error Exp (%.3f+-%.3f) - Image Shuffled(%.3f+-%.3f): P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",...
    mean(Eang_arr),sem(Eang_arr),mean(Eang_arr_shfl),sem(Eang_arr_shfl),P,TST.tstat,TST.df,CI(1),CI(2));
end
diary off


%% Summary single trial decoding result.
global sumdir
sumdir = fullfile(decodedir,'summary');
% diary(fullfile(sumdir,"singletr_Shuffle_summary.log"))
fprintf("Decoding in single trial noisy scenario:\n")
msklabels = BayDecStats(1).msklabels; 
DecTab = table();
for mski = 1:numel(msklabels)% F signif channel
chanN_arr = arrayfun(@(B)B.chan2dec(mski),BayDecStats);
Eang_arr = arrayfun(@(B)B.Dec_summary_sg(mski).EangD_m,BayDecStats);
EL2_arr = arrayfun(@(B)B.Dec_summary_sg(mski).EL2D_m,BayDecStats);
StrctCol = arrayfun(@(B)B.DecCol_sg{mski}, BayDecStats, 'Uni', 0);
idxmap = randperm(121); 
ErrSshfl = cellfun(@(S)error_stats_shfl(S,idxmap),StrctCol); % compute error after image index shuffling
Eang_arr_shfl = struct2table(ErrSshfl).EangD_m;
EL2_arr_shfl = struct2table(ErrSshfl).EL2D_m;
DecTab.(msklabels(mski)+"chanN") = chanN_arr;
DecTab.(msklabels(mski)+"Eang") = Eang_arr;
DecTab.(msklabels(mski)+"EL2") = EL2_arr;
DecTab.(msklabels(mski)+"Eang_shfl") = Eang_arr_shfl;
DecTab.(msklabels(mski)+"EL2_shfl") = EL2_arr_shfl;
[~,P,CI,TST] = ttest(Eang_arr,Eang_arr_shfl);
fprintf("%s (N=%.1f+-%.1f)",msklabels(mski),mean(chanN_arr),std(chanN_arr))
fprintf("Angle Error Exp (%.3f+-%.3f) - Image Shuffled(%.3f+-%.3f): P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",...
    mean(Eang_arr),sem(Eang_arr),mean(Eang_arr_shfl),sem(Eang_arr_shfl),P,TST.tstat,TST.df,CI(1),CI(2));

[~,P,CI,TST] = ttest(EL2_arr,EL2_arr_shfl);
fprintf("%s (N=%.1f+-%.1f)",msklabels(mski),mean(chanN_arr),std(chanN_arr))
fprintf("L2 Error Exp (%.3f+-%.3f) - Image Shuffled(%.3f+-%.3f): P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",...
    mean(EL2_arr),sem(EL2_arr),mean(EL2_arr_shfl),sem(EL2_arr_shfl),P,TST.tstat,TST.df,CI(1),CI(2));
end
%%
writetable(DecTab,fullfile(sumdir,"Both_DecodingStatsTab.csv"))
% diary off

%% Shuffling image index and measure the distance
S = DecCol_sg{1};
[STS] = error_stats_shfl(S);
[STS_shfl] = error_stats_shfl(S,randperm(121));
%% Use some stats above to do a comparison plot. 
%% Plot visualization
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
V1msk = (ExpTab_cmb.pref_chan <=48) & (ExpTab_cmb.pref_chan >= 33);
V4msk =  ExpTab_cmb.pref_chan > 48;
ITmsk =  ExpTab_cmb.pref_chan < 33;
Alfamsk = ExpTab_cmb.Animal=="Alfa";
Betomsk = ExpTab_cmb.Animal=="Beto";
sucs_msk = (ExpTab_cmb.t_p_initend<1E-3)&(ExpTab_cmb.t_p_initmax<1E-3);
%%
h=paired_stripe_plot({Eang_arr,Eang_arr_shfl}, ["Exp", "Shuffled"], "ExpAngleDist", ...
    "Expected Angle Distance cmp\nAll F significant units", "SingTrFsig_Shflcmp");
h=paired_stripe_plot({EL2_arr,EL2_arr_shfl}, ["Exp", "Shuffled"], "ExpL2Dist", ...
    "Expected L2 Distance cmp\nAll F significant units", "SingTrFsig_Shflcmp");
%%
h=paired_stripe_plot({Eang_arr,Eang_arr_shfl}, ["Exp", "Shuffled"], "ExpAngleDist", ...
    "Expected Angle Distance cmp\nAll F significant units", "SingTrFsig_Shflcmp_areasep", {V1msk,V4msk,ITmsk}, ["V1", "V4", "IT"]);
h=paired_stripe_plot({EL2_arr,EL2_arr_shfl}, ["Exp", "Shuffled"], "ExpL2Dist", ...
    "Expected L2 Distance cmp\nAll F significant units", "SingTrFsig_Shflcmp_areasep", {V1msk,V4msk,ITmsk}, ["V1", "V4", "IT"]);
%%
h=paired_stripe_plot({Eang_arr,Eang_arr_shfl}, ["Exp", "Shuffled"], "ExpAngleDist", ...
    "Expected Angle Distance cmp All F significant units", "SingTrFsig_Shflcmp_animsep", {Alfamsk, Betomsk}, ["Alfa","Beto"]);
h=paired_stripe_plot({EL2_arr,EL2_arr_shfl}, ["Exp", "Shuffled"], "ExpL2Dist", ...
    "Expected L2 Distance cmp All F significant units", "SingTrFsig_Shflcmp_animsep", {Alfamsk, Betomsk}, ["Alfa","Beto"]);
%%
h=paired_stripe_plot({Eang_arr,Eang_arr_shfl}, ["Exp", "Shuffled"], "ExpAngleDist", ...
    "Expected Angle Distance cmp All F significant units", "SingTrFsig_Shflcmp_sucssep", {sucs_msk, ~sucs_msk}, ["Success","Not Success"]);
% It seems evolution successfulness do not affect the success of decoding
% out manifold image!!!Though non successful evolution has larger error
% than successful. 
%%

function h=paired_stripe_plot(vars, groupnms, varnm, titstr, savestr, msks, labels)%, msks
if nargin<=5
msks = {[]};
labels=["All"];
end
global sumdir
h=figure; hold on ;set(h,'pos',[1000,360,450,620]);
% msk = msks{i};
% nExp = sum(msk);
Cord = colororder();
ncols = numel(vars); 
npnts = numel(vars{1});
xjit = 0.1*randn(npnts,1);
msk_all = ones(npnts,1,'logical');
varsmat = cat(2,vars{:});

    for mi = 1:numel(msks)
    msk = msks{mi}; 
    if isempty(msk), msk = msk_all; end
    for i = 1:ncols
    scatter(i+xjit(msk)',vars{i}(msk),'MarkerEdgeColor',Cord(mi,:), 'Displ', labels(mi));
    end
    [~,P,CI,TST] = ttest(vars{1}(msk),vars{2}(msk));
    titstr = titstr+compose("\n %s %s (%.3f+-%.3f) - %s(%.3f+-%.3f): P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]",...
        labels(mi),groupnms(1),mean(vars{1}(msk)),sem(vars{1}(msk)),groupnms(2),mean(vars{2}(msk)),sem(vars{2}(msk)),P,TST.tstat,TST.df,CI(1),CI(2));
    end

plot([1:ncols]'+xjit',varsmat(:,:)','color',[0,0,0,0.1],'HandleVis','off')
xticks(1:ncols);xticklabels(groupnms); %ylim(YLIM)
ylabel(varnm)
legend('Location','Best')
title(titstr)
savefn = compose("%s_%s",savestr,varnm);
savefig(h,fullfile(sumdir, savefn+".fig")) 
saveas(h,fullfile(sumdir, savefn+".png"))
saveas(h,fullfile(sumdir, savefn+".pdf"))

end

function [STS] = error_stats_shfl(S,idxmap)
% idxmap: a remapping of image index, such that the ground truth 
global L2mat angmat cosmat
if nargin==1, idxmap = 1:121;end
gtidxarr = arrayfun(@(A)A.imgi,S);
probmat = cell2mat(arrayfun(@(A)A.prob,S','uni',0));
[~,invidxmap] = sort(idxmap);
EangD_vec = sum(probmat(invidxmap,:).*angmat(:,idxmap(gtidxarr)),1);
EL2D_vec = sum(probmat(invidxmap,:).*L2mat(:,idxmap(gtidxarr)),1);
STS.EangD_m = mean(EangD_vec); 
STS.EangD_s =sem(EangD_vec);
STS.EL2D_m = mean(EL2D_vec); 
STS.EL2D_s =sem(EL2D_vec);
end

function [popmat, taridx_vec] = get_singletrial_data(act_col)
idx_arr = reshape(1:numel(act_col),11,11);
idx_arr(:,1)=1;idx_arr(:,11)=121;
taridx_col = arrayfun(@(i)repmat(idx_arr(i),1,size(act_col{i},2)),1:numel(act_col),'uni',0);
taridx_vec = cell2mat(taridx_col)';
popmat = cell2mat(reshape(act_col,1,[]))';
end

function ent = distr_entropy(prob)
    prob = prob ./ sum(prob);
    ent = -nansum(prob.*log(prob));
end

function [S,tab] = BayesDecManif(Xmean,Xstd,Y,X2dec,gtidx2dec,verbose)
    if nargin <=5, verbose = false;end
    if nargin ==3, 
        uniqidx = [1,12:110,121];
        X2dec = Xmean(uniqidx,:);
        gtidx2dec = uniqidx;
    end
    % Using Gaussian model of single neurons variability
    idmap = [ones(1,11),12:110,121*ones(1,11)]; 
    % too small std then that channel cannot signal that position.
    invalidmsk = Xstd<1E-5; 
    ndec = size(X2dec,1);
    S = repmat(struct(),ndec,1);
    for i = 1:ndec
        imgi = gtidx2dec(i);
        vec2decode = X2dec(i,:);
        gt_code = Y(imgi,:);
        % using naive bayes, compute posterior for each channel x each image
        naive_posterior = 0.5*((Xmean - vec2decode) ./ Xstd).^2 + log(Xstd)+ 0.5 * log(2*pi);
        naive_posterior(invalidmsk) = nan; 
        posterior = nansum(naive_posterior,2);
        prob = softmax(-posterior);
        [maxprob,maxidx] = max(prob);
        decidx = idmap(maxidx);
        L2dist = sqrt(sum((Y - gt_code).^2,2));
        angdist = acos((Y*gt_code') / norm(gt_code));
        EL2D = sum(L2dist .* prob) / sum(prob);
        EangD = sum(angdist .* prob) / sum(prob);
        prob_mrg = [sum(prob(1:11));prob(12:110);sum(prob(111:121))];
        ent = distr_entropy(prob_mrg);
        if verbose
        fprintf("L2 dist %.3f Angular dist %.3f Output entropy %.3f maxprob %.3f\n",EL2D,EangD,ent,maxprob)
        end
        S(i).EL2D = EL2D;
        S(i).EangD = EangD;
        S(i).ent = ent;
        S(i).maxprob = maxprob;
        S(i).decidx = decidx;
        S(i).imgi = imgi;
        S(i).correct = (imgi==decidx);
        S(i).prob = prob;
    end
    tab = struct2table(S);
    fprintf("Expected L2 dist %.3f (%.3f)\t",mean(tab.EL2D), sem(tab.EL2D))
    fprintf("Expected Angle dist %.3f (%.3f)\n",mean(tab.EangD), sem(tab.EangD))
    fprintf("Entropy %.3f (%.3f)\t",mean(tab.ent), sem(tab.ent))
    fprintf("Max probability %.3f (%.3f)\n",mean(tab.maxprob), sem(tab.maxprob))
    fprintf("Correct rate %.3f (%.3f)\n",mean(tab.correct), sem(tab.correct))
end

function [Scol,tabs,summary] = BayesDecwMasks(Xmean_all,Xstd_all,Y,msks,labels,X2dec,Yidx2dec)
if nargin <=5, sgtrial = false; 
else, sgtrial = true; fprintf("=======Single Trial Decoding========\n");end
Scol = {};
tabs = {};
summary = [];
for mi = 1:numel(msks)
    msk = msks{mi};
    if isempty(msk), msk = ones(size(Xmean_all, 2),1,'logical'); end
    fprintf("=======%s (Nfeat=%d)======\n",labels(mi),sum(msk))
    Xmean = Xmean_all(:,msk);
    Xstd = Xstd_all(:,msk);
    if ~sgtrial,
    [S,tab] = BayesDecManif(Xmean,Xstd,Y);
    else
    [S,tab] = BayesDecManif(Xmean,Xstd,Y,X2dec(:,msk),Yidx2dec);
    end
    Scol{end+1} = S;
    tabs{end+1} = tab;
    SUM = struct();
    SUM.mask = msks{mi};
    SUM.nfeat = sum(msk);
    SUM.EL2D_m = mean(tab.EL2D);
    SUM.EL2D_s = sem(tab.EL2D);
    SUM.EangD_m = mean(tab.EangD);
    SUM.EangD_s = sem(tab.EangD);
    SUM.ent_m = mean(tab.ent);
    SUM.ent_s = sem(tab.ent);
    SUM.maxprob_m = mean(tab.maxprob);
    SUM.maxprob_s = sem(tab.maxprob);
    SUM.correct_m = mean(tab.correct);
    SUM.correct_s = sem(tab.correct);
    summary = [summary,SUM];
end
end

function newObj = copyObject(obj)
% https://undocumentedmatlab.com/articles/general-use-object-copy
    try
        % R2010b or newer - directly in memory (faster)
        objByteArray = getByteStreamFromArray(obj);
        newObj = getArrayFromByteStream(objByteArray);
    end
end
