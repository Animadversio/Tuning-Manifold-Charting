%% Evol_ref_imgs_Dissim

% Animal = "Beto";
Animal = "Beto";Set_Path;
mat_dir = "O:\Mat_Statistics";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'),'Stats')
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'),'EStats')
load(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
load(fullfile(mat_dir,Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat")
%% Build a pool of ref images to search for
refroot = "N:\Stimuli\2019-Selectivity\2019-Selectivity-Big-Set-01";
refroot1 = "N:\Stimuli\2019-Selectivity\2019-Selectivity-Big-Set-01\forTesting";
refroot2 = "N:\Stimuli\2020-manifold-references";
% refimgfdrs = string(ls("N:\Stimuli\2019-Selectivity\chan*"));
refimgfdrs = arrayfun(@(F)"N:\Stimuli\2019-Selectivity\"+string(F.name),dir("N:\Stimuli\2019-Selectivity\chan*"));
search_dirs = [refimgfdrs; refroot; refroot1; refroot2];
allimgs = [];
for fdr = search_dirs'
files = dir((fdr+"\*"));
allrefimg = arrayfun(@(F)fdr+"\"+string(F.name),files);
allimgs = [allimgs;allrefimg];
end
%% PreProc: Find all the reference images and put them in place! 
for Expi = 1:numel(EStats)
fprintf("Process Experiment %d\n",Expi)
parts = split(EStats(Expi).meta.stimuli,"\");
evolpardir = fullfile(parts{1:end-1});
evorefs = arrayfun(@(F) evolpardir+"\"+string(F.name),dir(evolpardir+"\*"));
refimg_paths = [];
alreadythere = true;
for i = 1:numel(EStats(Expi).ref.imgnm)
    msk = contains(evorefs, EStats(Expi).ref.imgnm(i)+".");
    if sum(msk)==1
        refimg_paths = [refimg_paths, evorefs(msk)];
    elseif sum(msk) > 1
        keyboard
    else
        alreadythere = false;
    end
end
if ~ alreadythere
% EStats(Expi).ref.imgnm;
% allrefimg = string(ls(refroot+"\*"));
suf_idx = cell2mat(strfind(EStats(Expi).ref.imgnm, "_thread000_nat"));
imgnm_orig = arrayfun(@(nm,idx)extractBetween(nm,1,idx-1),EStats(Expi).ref.imgnm,suf_idx);
refimg_paths = [];
for i = 1:numel(imgnm_orig)
    % Search for the images with 
    imnm_part = strcat(imgnm_orig(i),'.');%EStats(Expi).ref.imgnm{i}(1:suf_idx(i)-1);
    evorefmsk = contains(evorefs, imnm_part);
    if sum(evorefmsk) > 0 % If we can find it in the parent folder
        imgpath = evorefs(evorefmsk);
        if numel(imgpath)==0, 
            keyboard; 
        elseif numel(imgpath)>1, 
            disp(imgpath); % Multiple image names contains the search str. 
            nmlen = cellfun(@(S)length(S),imgpath);
            [minlen, minidx] = min(nmlen); % choose the shortest name 
            disp(minidx)
            %  keyboard; 
            imgpath = imgpath(minidx);
        end
    else %
        fullpath_msk = contains(allimgs,imnm_part);
        imgpath = allimgs(fullpath_msk);
        if numel(imgpath)==0, 
            keyboard; 
        elseif numel(imgpath)>1, 
            disp(imgpath); % Multiple image names contains the search str. 
            nmlen = cellfun(@(S)length(S),imgpath);
            [minlen, minidx] = min(nmlen); % choose the shortest name 
            disp(minidx)
            % keyboard; 
            imgpath = imgpath(minidx);
        end
    end
    refimg_paths = [refimg_paths; imgpath];
    suffix = split(imgpath,".");
    tarname = fullfile(evolpardir,strcat(EStats(Expi).ref.imgnm(i),".",suffix{end}));
    disp(tarname+"\n");
    copyfile(imgpath, tarname)
end
end
assert(numel(refimg_paths)==numel(EStats(Expi).ref.imgnm)) % Make sure the number of found path is the same as original ones. 
EStats(Expi).ref.impaths = refimg_paths;
end
%%
save(fullfile(mat_dir, Animal+'_Evol_stats.mat'),'EStats')
%% Process: Compute the Dissimilarity
EvoRefImDistStat = repmat(struct(),1,numel(EStats));
D = torchImDist("squeeze");
D=D.load_net();
%%
for Expi = 1:numel(EStats)
parts = split(EStats(Expi).meta.stimuli,"\");
evolpardir = fullfile(parts{1:end-1});
evorefs = arrayfun(@(F) evolpardir+"\"+string(F.name),dir(evolpardir+"\*"));
refimgs = {};
for i = 1:numel(EStats(Expi).ref.imgnm)
    msk = contains(evorefs, EStats(Expi).ref.imgnm(i)+".");
    img = imread(evorefs(msk));
    [H,W,C] = size(img);
	if C==1, img = repmat(img,1,1,3); end
    refimgs{i} = imresize(img, [256,256]);
end
refimgtsr = cat(4,refimgs{:});
tic
D = D.select_metric("squeeze");
distMat = D.distmat_B(refimgtsr,50); % batch process is faster here
EvoRefImDistStat(Expi).squ = distMat;
D.metric = "L2";
EvoRefImDistStat(Expi).L2 = D.distmat(refimgtsr);
D.metric = "FC6_L2";
EvoRefImDistStat(Expi).FC6 = D.distmat(refimgtsr);
D.metric = "FC6_corr";
EvoRefImDistStat(Expi).FC6_corr = D.distmat(refimgtsr);
D = D.select_metric("SSIM"); 
EvoRefImDistStat(Expi).SSIM = D.distmat(refimgtsr);
toc %% SSIM can take quite a while 
% keyboard;
end
%%
save(fullfile(mat_dir,Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat")
%% Summary Part
%% Compute Radial tuning curves, corr and AUC and save to stats 
Animal = "Alfa";Set_Path;
mat_dir = "O:\Mat_Statistics";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'),'Stats') % Manifold Stats 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'),'EStats') % Evolution Stats
load(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats') % Collection of Radial Tuning Curves in all Spaces
load(fullfile(mat_dir,Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat") % Image Distance between reference images. 
%% Adding the `evoref` field to `RadTuneStats` collecting the stats in this space. 
metric_list = ["squ","SSIM","L2","FC6","FC6_corr"]; % all the metrics we want to use to measure distance 
label_list = ["LPIPS (SqueezeNet)", "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
tic
for Expi=1:numel(Stats)
% Manifold images. 
fprintf("Processing %s Exp %d\n", Animal,Expi)
RadTuneStats(Expi).Animal = Animal; 
RadTuneStats(Expi).Expi = Expi;

score_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
[sortScore,sortId] = sort(score_vec,'Descend');
[maxScore,maxId] = max(score_vec);
RadTuneStats(Expi).evoref.maxScore = maxScore;
for mi = 1:numel(metric_list) % Loop through different distance metrics
	metname = metric_list(mi);
	[gpr_fit, ~, AOC, gprMdl] = GPRfitting(EvoRefImDistStat(Expi).(metname)(:,maxId), score_vec); % Gaussian Process Smoothing or Fitting
	[R2, slope, intercept] = Linearfitting(EvoRefImDistStat(Expi).(metname)(:,maxId), score_vec);
	% titstr{mi} = { compose("Manif: pear %.3f spear %.3f",
	RadTuneStats(Expi).evoref.AOC.(metname) = AOC; 
	RadTuneStats(Expi).evoref.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).evoref.gpr_fit.(metname) = gpr_fit;
	RadTuneStats(Expi).evoref.linR2.(metname) = R2;
	RadTuneStats(Expi).evoref.lincc.(metname) = [slope, intercept];
	RadTuneStats(Expi).evoref.corr.(metname) = corr(EvoRefImDistStat(Expi).(metname)(:,maxId),score_vec);
	RadTuneStats(Expi).evoref.corr_sp.(metname) = corr(EvoRefImDistStat(Expi).(metname)(:,maxId),score_vec,'Type','Spearman');
    fprintf("%d %s corr %.3f ",Expi,metname,RadTuneStats(Expi).evoref.corr.(metname))
end
toc
end
%% Save the Radial Tuning Stats back to the mat file. 
save(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
%  For comparisons of interplotation and AUC integration methods see
%  `vis_cmp_smooth_method.m` 


function [gpr_fit, xlinsp, AOC, gprMdl] = GPRfitting(X, Y)
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AOC = trapz(xlinsp, gpr_fit);  % integrate under the gaussian process fitting curve. 
% Plot curve on a fig
% plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end


function [smth_fit, xlinsp, AOC, smthcrv] = BinSmoothing(X, Y, interpmethod)
if nargin==2, varargin={}; interpmethod = 'linear'; 
elseif nargin==3, varargin={}; end
smthcrv = smooth(X,Y);
xlinsp = linspace(0,max(X),100);
[uniqX, iX, iUniq] = unique(X); % get rid of redundancy in X if there is any.
smth_fit = interp1(X(iX),smthcrv(iX),xlinsp,interpmethod,'extrap'); %'spline'
AOC = trapz(xlinsp, smth_fit);  % integrate under the gaussian process fitting curve. 
% % Plot curve on a fig
% scatter(X, smthcrv, 'DisplayName', 'Data Smooth')
% plot(xlinsp, smth_fit,'Disp',strcat('smooth+interp ',interpmethod), varargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end

function [r,m,b] = Linearfitting(X, Y, varargin) % Util function to add a linear reg line to a scatter
if nargin==2, varargin={};end
[r,m,b] = regression(X, Y);
% Plot line on a fig
% xmin = min(X); xmax = max(X);
% p = plot([xmin,xmax],[xmin,xmax].*m+b,varargin{:});
% set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end