%% Collect Best Image Representation. Build `ReprStats`
%  Feb.2nd,2021. Seems these images are a little bit "biases" maybe I should use better examples by averaging them.
%  Save the best image from the manifold space and the best block averaged
%  image from the evolution.
%%
Animal="Beto"; 
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%% 
ReprStats = repmat(struct(),1,length(EStats));
%% Best Evol Images (single trial averaged)
ExpType = "Evol";
for Expi = 1:length(EStats)
    [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, struct());
    [maxscore, maxidx] = max(score_vect(:,end));
    ReprStats(Expi).(ExpType).BestImg = imread(imgfullnm_vect(maxidx));
end
%% Best Manif Images (Multiple averaged)
ExpType = "Manif";
for Expi = 1:length(Stats)
    [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, struct());
    [maxscore, maxidx] = max(score_vect(:,end));
    ReprStats(Expi).(ExpType).BestImg = imread(imgfullnm_vect(maxidx));
end
%% Best Evol Images (block averaged)
G = FC6Generator();
%%
ExpType = "Evol";
for Expi = 1:length(EStats)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
ui=1;
blockscores = cellfun(@(psth)mean(psth(ui,51:200,:),[2,3]),EStats(Expi).evol.psth);
block_std = cellfun(@(psth)std(mean(psth(ui,51:200,:),[2]),1,3),EStats(Expi).evol.psth);
block_sem = cellfun(@(psth)sem(mean(psth(ui,51:200,:),[2]),3),EStats(Expi).evol.psth);
[maxblkscore, maxblkidx] = max(blockscores);
[codes_block,~,~] = load_codes_all(EStats(Expi).meta.stimuli,1,maxblkidx-1);
code_mean = mean(codes_block,1);
avgimg = G.visualize(code_mean);
ReprStats(Expi).(ExpType).BestBlockAvgImg = avgimg;
end
%%
save(fullfile(MatStats_path, compose("%s_ImageRepr.mat", Animal)),'ReprStats')
function [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags)
si=1;ui=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
end
imgN = length(imgnm_vect);
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{2};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix),imgnm_vect);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
if ExpType == "Manif"
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false); % note there is trial averaging here. 
psth_all = reshape(cell2mat(psth_all),imgN,[]);
elseif ExpType == "Evol"
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
end
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end