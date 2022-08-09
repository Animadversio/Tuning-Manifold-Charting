%% Manifold Bayes Decoding
%  Final version of bayes decoding 
%  Decoding example using one experiment.
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

%%
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"), 'MapVarStats')
%%
Expi = 3;si = 1;
% Compute Fstats and create masks
spikeID = MapVarStats(Expi).units.spikeID;
nCh = numel(MapVarStats(Expi).units.spikeID);
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
FStats = cellfun(@anova_cells,actmap_col); 
Fmsk = struct2table(FStats).F_P < 0.001; 
V1msk = spikeID <=48 & spikeID>=33;
V4msk = spikeID <=64 & spikeID>=49;
ITmsk = spikeID <=32 & spikeID>=1;
%%
actmean = cellfun(@(A)mean(A,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actmean = cell2mat(reshape(actmean,1,[]));
actstd = cellfun(@(A)std(A,1,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actstd = cell2mat(reshape(actstd,1,[]));
actmeanmat = actmean';
actstdmat = actstd';
actmean_sel = actmeanmat(:,Fmsk);
actstd_sel = actstdmat(:,Fmsk);
% [zact_sel, chmean, chstd] = zscore(actmean_sel,1);
%% Bayes decode.
idlist = [1,12:110,121];
imgi = 120;
vec2decode = actmean_sel(imgi,:);
gt_code = targ_coord(imgi,:);
% too small std then that channel cannot signal that position.
invalidmsk = actstd_sel<1E-5; 
% using naive bayes, compute posterior for each channel x each image
naive_posterior = 0.5*((actmean_sel - vec2decode) ./ actstd_sel).^2 + log(actstd_sel)+ 0.5 * log(2*pi);
naive_posterior(invalidmsk) = nan; 
posterior = nansum(naive_posterior,2);
prob = softmax(-posterior);
[maxprob,maxidx] = max(prob);
L2dist = sqrt(sum((targ_coord - gt_code).^2,2));
angdist = acos((targ_coord*gt_code') / norm(gt_code));
EL2D = sum(L2dist .* prob) / sum(prob);
EangD = sum(angdist .* prob) / sum(prob);
prob_mrg = [sum(prob(1:11));prob(12:110);sum(prob(111:121))];
ent = distr_entropy(prob_mrg);
fprintf("L2 dist %.3f Angular dist %.3f Output entropy %.3f\n",EL2D,EangD,ent)
%% Summary for each exp
[S,tab] = BayesDecManif(actmean_sel,actstd_sel,targ_coord);
%%
[Scol,tabs,summary] = BayesDecwMasks(actmeanmat,actstdmat,targ_coord,{V1msk,V4msk,ITmsk}, ["V1", "V4", "IT"]);
% [Scol,tabs,summary] = BayesDecwMasks(actmeanmat,actstdmat,targ_coord,{V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk}, ["V1 Fsig", "V4 Fsig", "IT Fsig"]);
% [Scol,tabs,summary] = BayesDecwMasks(actmeanmat,actstdmat,targ_coord,{Fmsk,~Fmsk},["Fsignif","non Fsignif"]);

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

