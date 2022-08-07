%% Create the summary plots for Compare evolution scores for Alfa/Beto.
%  See if the Evolution works
%%
Animal = "Alfa"; Set_Path;
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, Animal+"_RDEvol_stats.mat"), 'RDStats')
%% summary stats
ExpType = "RDEvol";
for Expi = 1:length(RDStats)
fprintf("Processing Reduced Dimension Evolution Exp %d\n",Expi)
thread_n = RDStats(Expi).evol.thread_num;
% the following snippet is to get rid of 0 unit (null channel)
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == 1) & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
ui = find(prefchan_id==chid);
Window = 50:200;
% PSTH for evolution
evol_stim_fr = cellfun(@(psth)mean(psth,3),RDStats(Expi).evol.psth,'Uni',false);
evol_stim_fr = cell2mat(reshape(evol_stim_fr',1,1,[],thread_n));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),RDStats(Expi).evol.psth,'Uni',false);
evol_stim_sem = cell2mat(reshape(evol_stim_sem',1,1,[],thread_n));
% Scalor score for evolution
score_avg = cellfun(@(psth)mean(psth(:,51:200,:),'all') - mean(psth(:,1:50,:),'all'),RDStats(Expi).evol.psth);
score_sem = cellfun(@(psth)std(squeeze(mean(psth(:,51:200,:),[1,2])))...
    /sqrt(size(psth,3)),RDStats(Expi).evol.psth); 
% Generation number of 50% progress
perct = .5; gen50 = sum(score_avg < ( perct * max(score_avg,[],2) + (1-perct) * score_avg(:,1)),2);
end
%% Get the prefer channel and unit numbering array
prefchan_arr = arrayfun(@(R)R.evol.pref_chan(1),RDStats);
unitnum_arr = zeros(1,length(RDStats));
for Expi = 1:length(RDStats)
% the following snippet is to get rid of 0 unit (null channel)
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == 1) & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
unitnum_arr(Expi) = find(prefchan_id==chid);
end

%% MAIN LOOP: Collect the Paired Exp Statistics into a structure / table to plot
savedir = "E:\OneDrive - Washington University in St. Louis\Evol_ReducDim\summary";
RDThreadStats = repmat(struct(),1); % compress the RDStats even more, Statistics for indiv thread
RDEvolCmpStats = repmat(struct(),1); % compare 2 threads if they are both evolving the same unit using Full with different optimizer.
ipair = 1; % entry counter for `RDEvolCmpStats` 
for Expi = 1:length(RDStats)
% the following snippet is to get rid of 0 unit (null channel)
ui = unitnum_arr(Expi);
score_avg = cellfun(@(psth)mean(psth(ui,51:200,:),'all') - mean(psth(ui,1:50,:),'all'), RDStats(Expi).evol.psth);
blockn = RDStats(Expi).evol.block_n;
threadn = RDStats(Expi).evol.thread_num;

unit_psth = cellfun(@(psth)psth(ui,:,:), RDStats(Expi).evol.psth, 'Uni', false);
window = [51:200];bslwdw = [1:50];
score_avg = cellfun(@(psth)mean(psth(:,window,:),'all') - mean(psth(:,bslwdw,:),'all'), unit_psth);
% Generation number of 50%, 65%, 80% progress
perct = .50; RDThreadStats(Expi).gen50 = sum(score_avg < ( perct * max(score_avg,[],2) + (1-perct) * score_avg(:,1)),2)';
perct = .65; RDThreadStats(Expi).gen65 = sum(score_avg < ( perct * max(score_avg,[],2) + (1-perct) * score_avg(:,1)),2)';
perct = .80; RDThreadStats(Expi).gen80 = sum(score_avg < ( perct * max(score_avg,[],2) + (1-perct) * score_avg(:,1)),2)';
for thr_i = 1:threadn % collect activation data for each thread.
block_mean_score = cellfun(@(psth)mean(psth(1,window,:),'all'),unit_psth(thr_i,1:blockn-1));
[~,peakBlock]=max(block_mean_score,[],2);
if peakBlock == blockn-1, peakBlock = blockn-2; end % avoid the last generation.
peakpsths = cell2mat(reshape(unit_psth(thr_i,peakBlock:peakBlock+1),1,1,[]));
initpsths = cell2mat(reshape(unit_psth(thr_i,2:3),1,1,[]));
initacts = squeeze(mean(initpsths(1,window,:),2));
peakacts = squeeze(mean(peakpsths(1,window,:),2));
[H,P,CI,STATS] = ttest2(peakacts, initacts); % Peak compare to the initial
RDThreadStats(Expi).textacts{thr_i} = squeeze(mean(unit_psth{thr_i,1}(1, window, :), [1,2]));
RDThreadStats(Expi).initacts{thr_i} = initacts;
RDThreadStats(Expi).peakacts{thr_i} = peakacts;
RDThreadStats(Expi).textscore(thr_i) = mean(unit_psth{thr_i,1}(1, window, :),'all') - mean(unit_psth{thr_i,1}(1, 1:50, :),'all');
RDThreadStats(Expi).initscore(thr_i) = mean(initpsths(1, window, :),'all') - mean(initpsths(1, 1:50, :),'all');
RDThreadStats(Expi).peakscore(thr_i) = mean(peakpsths(1, window, :),'all') - mean(peakpsths(1, 1:50, :),'all');
RDThreadStats(Expi).tstat(thr_i) = STATS.tstat; % Peak compare to the initial t stats
RDThreadStats(Expi).DAOA(thr_i) = (mean(peakacts) - mean(initacts)) / mean(initacts) ;
RDThreadStats(Expi).DAOA2(thr_i) = (mean(peakacts) - mean(unit_psth{thr_i,1}(1, window, :),'all')) / mean(unit_psth{thr_i,1}(1, window, :),'all') ;
% RDEvolCmpStats(Expi).DAOA(thr_i) = (RDEvolCmpStats(Expi).peakscore(thr_i) - RDEvolCmpStats(Expi).initscore(thr_i)) / RDEvolCmpStats(Expi).initscore(thr_i);
% RDEvolCmpStats(Expi).DAOA2(thr_i) = (RDEvolCmpStats(Expi).peakscore(thr_i) - RDEvolCmpStats(Expi).textscore(thr_i)) / RDEvolCmpStats(Expi).initscore(thr_i);
end
if threadn == 1
    thr_i = 2;
    RDThreadStats(Expi).initscore(thr_i) = nan;
    RDThreadStats(Expi).peakscore(thr_i) = nan;
    RDThreadStats(Expi).tstat(thr_i) = nan;
    RDThreadStats(Expi).DAOA(thr_i) = nan;
    RDThreadStats(Expi).DAOA2(thr_i) = nan;
end
% Do comparison of activation and collect comparison stats for threads
% Collect the mean and sem for the compared 2 distribution, for visualization purpose.
if ~all(RDStats(Expi).evol.optim_names == ["ZOHA Sphere lr euclid", "ZOHA Sphere lr euclid ReducDim"])%(RDStats(Expi).evol.optim_names(1) ~= "ZOHA Sphere lr euclid") || (RDStats(Expi).evol.optim_names(2) ~= "ZOHA Sphere lr euclid ReducDim")
    fprintf("Exp %d optimizers not standard\n",Expi)
    disp(RDStats(Expi).evol.optim_names)
    keyboard
    % assume both threads here are doing Reduced Dimen Evol,
    % And the last Exp contains a full space evolution. 
    assert(all(RDStats(Expi).evol.optim_names == "ZOHA Sphere lr euclid ReducDim"))
    assert(all(RDStats(Expi-1).evol.pref_chan == RDStats(Expi).evol.pref_chan))
    assert(all(RDStats(Expi-1).evol.optim_names(1) == "ZOHA Sphere lr euclid"))
    for thr_i = 1:threadn
        RDEvolCmpStats(ipair).pref_chan = RDStats(Expi-1).evol.pref_chan(1);
        [H,P,CI,STATS] = ttest2(RDThreadStats(Expi-1).peakacts{1}, RDThreadStats(Expi).peakacts{thr_i});
        RDEvolCmpStats(ipair).peak_cmp_t = STATS.tstat;
        RDEvolCmpStats(ipair).peak_cmp_p = P;
        RDEvolCmpStats(ipair).peak_mean = [mean(RDThreadStats(Expi-1).peakacts{1}), mean(RDThreadStats(Expi).peakacts{thr_i})];
        RDEvolCmpStats(ipair).peak_sem = [sem(RDThreadStats(Expi-1).peakacts{1}), sem(RDThreadStats(Expi).peakacts{thr_i})];
        [H,P,CI,STATS] = ttest2(RDThreadStats(Expi-1).initacts{1}, RDThreadStats(Expi).initacts{thr_i});
        RDEvolCmpStats(ipair).init_cmp_t = STATS.tstat;
        RDEvolCmpStats(ipair).init_cmp_p = P;
        RDEvolCmpStats(ipair).init_mean = [mean(RDThreadStats(Expi-1).initacts{1}), mean(RDThreadStats(Expi).initacts{thr_i})];
        RDEvolCmpStats(ipair).init_sem = [sem(RDThreadStats(Expi-1).initacts{1}), sem(RDThreadStats(Expi).initacts{thr_i})];
        [H,P,CI,STATS] = ttest2(RDThreadStats(Expi-1).textacts{1}, RDThreadStats(Expi).textacts{thr_i});
        RDEvolCmpStats(ipair).text_cmp_t = STATS.tstat;
        RDEvolCmpStats(ipair).text_cmp_p = P;
        RDEvolCmpStats(ipair).text_mean = [mean(RDThreadStats(Expi-1).textacts{1}), mean(RDThreadStats(Expi).textacts{thr_i})];
        RDEvolCmpStats(ipair).text_sem = [sem(RDThreadStats(Expi-1).textacts{1}), sem(RDThreadStats(Expi).textacts{thr_i})];
        RDEvolCmpStats(ipair).gen50 = [RDThreadStats(Expi-1).gen50(1), RDThreadStats(Expi).gen50(thr_i)];
        RDEvolCmpStats(ipair).gen65 = [RDThreadStats(Expi-1).gen65(1), RDThreadStats(Expi).gen65(thr_i)];
        RDEvolCmpStats(ipair).gen80 = [RDThreadStats(Expi-1).gen80(1), RDThreadStats(Expi).gen80(thr_i)];
        ipair = ipair + 1;
    end
else
    RDEvolCmpStats(ipair).pref_chan = RDStats(Expi).evol.pref_chan(1);
    [H,P,CI,STATS] = ttest2(RDThreadStats(Expi).peakacts{:});
    RDEvolCmpStats(ipair).peak_cmp_t = STATS.tstat;
    RDEvolCmpStats(ipair).peak_cmp_p = P;
    RDEvolCmpStats(ipair).peak_mean = cellfun(@mean, RDThreadStats(Expi).peakacts);
    RDEvolCmpStats(ipair).peak_sem = cellfun(@sem, RDThreadStats(Expi).peakacts);
    [H,P,CI,STATS] = ttest2(RDThreadStats(Expi).initacts{:});
    RDEvolCmpStats(ipair).init_cmp_t = STATS.tstat;
    RDEvolCmpStats(ipair).init_cmp_p = P;
    RDEvolCmpStats(ipair).init_mean = cellfun(@mean, RDThreadStats(Expi).initacts);
    RDEvolCmpStats(ipair).init_sem = cellfun(@sem, RDThreadStats(Expi).initacts);
    [H,P,CI,STATS] = ttest2(RDThreadStats(Expi).textacts{:});
    RDEvolCmpStats(ipair).text_cmp_t = STATS.tstat;
    RDEvolCmpStats(ipair).text_cmp_p = P;
    RDEvolCmpStats(ipair).text_mean = cellfun(@mean, RDThreadStats(Expi).textacts);
    RDEvolCmpStats(ipair).text_sem = cellfun(@sem, RDThreadStats(Expi).textacts);
    RDEvolCmpStats(ipair).gen50 = RDThreadStats(Expi).gen50;
    RDEvolCmpStats(ipair).gen65 = RDThreadStats(Expi).gen65;
    RDEvolCmpStats(ipair).gen80 = RDThreadStats(Expi).gen80;
    ipair = ipair + 1;
end
end
RDEvolCmpTable = struct2table(RDEvolCmpStats);
%%
writetable(RDEvolCmpTable, fullfile(savedir,compose("%s_RDEvolcompare.csv",Animal)))

%% Statistics for activations each group, recorded in RDEvolCmpStats 
initacts_m = RDEvolCmpTable.init_mean;
initacts_sem = RDEvolCmpTable.init_sem;
peakacts_m = RDEvolCmpTable.peak_mean;
peakacts_sem = RDEvolCmpTable.peak_sem;
textacts_m = RDEvolCmpTable.text_mean;
textacts_sem = RDEvolCmpTable.text_sem;
npair = size(RDEvolCmpTable,1);
%% Alfa, create mask from the CmpStats table.
V1msk = RDEvolCmpTable.pref_chan <=48 & RDEvolCmpTable.pref_chan >= 33;
V4msk = RDEvolCmpTable.pref_chan <=64 & RDEvolCmpTable.pref_chan >= 49;
ITmsk = RDEvolCmpTable.pref_chan <=32 & RDEvolCmpTable.pref_chan >= 1;
msk.V1 = V1msk;msk.V4 = V4msk; msk.IT = ITmsk; msk.all = ones(1,length(RDEvolCmpStats),'logical');
%%
figure(1);clf;set(1,'pos',[725         343        1431         467]);
ax = subtightplot(1,1,1,0.05, 0.1, 0.05);hold on
% errorbar([1:36]'+[0, 0.3], textacts_m, textacts_sem,'o')
scatter([1:npair]', initacts_m(:,1), 64,'filled')
scatter([1:npair]'+0.3, initacts_m(:,2), 64,'filled')
scatter([1:npair]', peakacts_m(:,1), 64,'filled')
scatter([1:npair]'+0.3, peakacts_m(:,2), 64,'filled')
errorbar([1:npair]'+[0, 0.3], initacts_m, 2*initacts_sem,'Marker','none','LineStyle','none','color','magenta')
errorbar([1:npair]'+[0, 0.3], peakacts_m, 2*peakacts_sem,'Marker','none','LineStyle','none','color','magenta')
line([1:npair]+[0,0,0.3,0.3]', [initacts_m(:,1)';...
                             peakacts_m(:,1)';...
                             peakacts_m(:,2)';...
                             initacts_m(:,2)'],'color', 'k')
ylabel("Firing Rate");xlabel("Pair id");title(compose("%s activation comparison for all experiments",Animal))
legend(["init rate full", "init rate 50D", "peak rate full", "peak rate 50D"])
saveas(1, fullfile(savedir,compose("%s_ActivDiff_t_allexp.jpg",Animal)))
savefig(1, fullfile(savedir,compose("%s_ActivDiff_t_allexp.fig",Animal)))
%%
figure(2);clf;set(2,'position',[440         473        1233         505])
subtightplot(1,3,1,0.05, 0.1, 0.05);hold on
errorbar([1:sum(msk.V1)]'+[0, 0.3], initacts_m(msk.V1,:), 2*initacts_sem(msk.V1,:),'Marker','o','LineStyle','none','color','blue')
errorbar([1:sum(msk.V1)]'+[0, 0.3], peakacts_m(msk.V1,:), 2*peakacts_sem(msk.V1,:),'Marker','o','LineStyle','none','color','red')
line([1:sum(msk.V1)]+[0,0,0.3,0.3]', [initacts_m(msk.V1,1)';...
		                             peakacts_m(msk.V1,1)';...
		                             peakacts_m(msk.V1,2)';...
		                             initacts_m(msk.V1,2)'],'color', 'k')
ylabel("Firing Rate");xlabel("pref chan");title("V1 Exps");xticks([1:sum(msk.V1)]);xticklabels(RDEvolCmpTable.pref_chan(V1msk))
subtightplot(1,3,2,0.05, 0.1, 0.05);hold on
errorbar([1:sum(msk.V4)]'+[0, 0.3], initacts_m(msk.V4,:), 2*initacts_sem(msk.V4,:),'Marker','o','LineStyle','none','color','blue')
errorbar([1:sum(msk.V4)]'+[0, 0.3], peakacts_m(msk.V4,:), 2*peakacts_sem(msk.V4,:),'Marker','o','LineStyle','none','color','red')
line([1:sum(msk.V4)]+[0,0,0.3,0.3]', [initacts_m(msk.V4,1)';...
		                             peakacts_m(msk.V4,1)';...
		                             peakacts_m(msk.V4,2)';...
		                             initacts_m(msk.V4,2)'],'color', 'k')
xlabel("pref chan");title("V4 Exps");xticks([1:sum(msk.V4)]);xticklabels(RDEvolCmpTable.pref_chan(V4msk))%ylabel("Firing Rate");
subtightplot(1,3,3,0.05, 0.1, 0.05);hold on
errorbar([1:sum(msk.IT)]'+[0, 0.3], initacts_m(msk.IT,:), 2*initacts_sem(msk.IT,:),'Marker','o','LineStyle','none','color','blue')
errorbar([1:sum(msk.IT)]'+[0, 0.3], peakacts_m(msk.IT,:), 2*peakacts_sem(msk.IT,:),'Marker','o','LineStyle','none','color','red')
line([1:sum(msk.IT)]+[0,0,0.3,0.3]', [initacts_m(msk.IT,1)';...
		                             peakacts_m(msk.IT,1)';...
		                             peakacts_m(msk.IT,2)';...
		                             initacts_m(msk.IT,2)'],'color', 'k')
xlabel("pref chan");title("IT Exps");xticks([1:sum(msk.IT)]);xticklabels(RDEvolCmpTable.pref_chan(ITmsk))%ylabel("Firing Rate");
suptitle(compose("%s activation comparison for all experiments",Animal))
legend(["init rate full", "init rate 50D", "peak rate full", "peak rate 50D"],'location','best')
saveas(2, fullfile(savedir,compose("%s_ActivDiff_lines_area.jpg",Animal)))
savefig(2, fullfile(savedir,compose("%s_ActivDiff_lines_area.fig",Animal)))
%% Evolution time
figure(3);clf
subtightplot(1,3,1,0.05, 0.1, 0.05);hold on
plot([1:sum(msk.V1)]'+[0  ], RDEvolCmpTable.gen50(msk.V1,1), 'Marker','o','color','red','LineStyle','none')
plot([1:sum(msk.V1)]'+[0.7], RDEvolCmpTable.gen50(msk.V1,2), 'Marker','o','color','blue','LineStyle','none')
% plot([1:sum(msk.V1)]'+[0, 0.7], RDEvolCmpTable.gen65(msk.V1,:), 'Marker','o','color','magenta','LineStyle','none')
% plot([1:sum(msk.V1)]'+[0, 0.7], RDEvolCmpTable.gen80(msk.V1,:), 'Marker','o','color','k','LineStyle','none')
line([1:sum(msk.V1)]+[0,0.7]', [ RDEvolCmpTable.gen50(msk.V1,1)';...
                                 RDEvolCmpTable.gen50(msk.V1,2)'],'color', 'k')
% line([1:sum(msk.V1)]+[0,0.7]', [ RDEvolCmpTable.gen65(msk.V1,1)';...
%                                  RDEvolCmpTable.gen65(msk.V1,2)'],'color', 'k')
% line([1:sum(msk.V1)]+[0,0.7]', [ RDEvolCmpTable.gen80(msk.V1,1)';...
%                                  RDEvolCmpTable.gen80(msk.V1,2)'],'color', 'k')
ylabel("Convergence Time");xlabel("pref chan");title("V1 Exps");xticks([1:sum(msk.V1)]);xticklabels(RDEvolCmpTable.pref_chan(V1msk))
subtightplot(1,3,2,0.05, 0.1, 0.05);hold on
plot([1:sum(msk.V4)]'+[0  ], RDEvolCmpTable.gen50(msk.V4,1), 'Marker','o','color','red','LineStyle','none')
plot([1:sum(msk.V4)]'+[0.7], RDEvolCmpTable.gen50(msk.V4,2), 'Marker','o','color','blue','LineStyle','none')
% plot([1:sum(msk.V4)]'+[0, 0.7], RDEvolCmpTable.gen65(msk.V4,:), 'Marker','o','color','magenta','LineStyle','none')
% plot([1:sum(msk.V4)]'+[0, 0.7], RDEvolCmpTable.gen80(msk.V4,:), 'Marker','o','color','k','LineStyle','none')
line([1:sum(msk.V4)]+[0,0.7]', [ RDEvolCmpTable.gen50(msk.V4,1)';...
                                 RDEvolCmpTable.gen50(msk.V4,2)'],'color', 'k')
% line([1:sum(msk.V4)]+[0,0.7]', [ RDEvolCmpTable.gen65(msk.V4,1)';...
%                                  RDEvolCmpTable.gen65(msk.V4,2)'],'color', 'k')
% line([1:sum(msk.V4)]+[0,0.7]', [ RDEvolCmpTable.gen80(msk.V4,1)';...
%                                  RDEvolCmpTable.gen80(msk.V4,2)'],'color', 'k')
ylabel("Convergence Time");xlabel("pref chan");title("V4 Exps");xticks([1:sum(msk.V4)]);xticklabels(RDEvolCmpTable.pref_chan(V4msk))
subtightplot(1,3,3,0.05, 0.1, 0.05);hold on
plot([1:sum(msk.IT)]'+[0  ], RDEvolCmpTable.gen50(msk.IT,1), 'Marker','o','color','red','LineStyle','none')
plot([1:sum(msk.IT)]'+[0.7], RDEvolCmpTable.gen50(msk.IT,2), 'Marker','o','color','blue','LineStyle','none')
% plot([1:sum(msk.IT)]'+[0, 0.7], RDEvolCmpTable.gen65(msk.IT,:), 'Marker','o','color','magenta','LineStyle','none')
% plot([1:sum(msk.IT)]'+[0, 0.7], RDEvolCmpTable.gen80(msk.IT,:), 'Marker','o','color','k','LineStyle','none')
line([1:sum(msk.IT)]+[0,0.7]', [ RDEvolCmpTable.gen50(msk.IT,1)';...
                                 RDEvolCmpTable.gen50(msk.IT,2)'],'color', 'k')
% line([1:sum(msk.IT)]+[0,0.7]', [ RDEvolCmpTable.gen65(msk.IT,1)';...
%                                  RDEvolCmpTable.gen65(msk.IT,2)'],'color', 'k')
% line([1:sum(msk.IT)]+[0,0.7]', [ RDEvolCmpTable.gen80(msk.IT,1)';...
%                                  RDEvolCmpTable.gen80(msk.IT,2)'],'color', 'k')
ylabel("Convergence Time");xlabel("pref chan");title("IT Exps");xticks([1:sum(msk.IT)]);xticklabels(RDEvolCmpTable.pref_chan(ITmsk))
suptitle(compose("%s activation comparison for all experiments",Animal))
% legend(["Iter 50% 4096D", "Iter 50% 50D", "Iter 65% 4096D", "Iter 65% 50D", "Iter 80% 4096D", "Iter 80% 50D"],'location','best')
legend(["Iter 50% 4096D", "Iter 50% 50D"],'location','best')
saveas(3, fullfile(savedir,compose("%s_EvolTime_lines_area.jpg",Animal)))
savefig(3, fullfile(savedir,compose("%s_EvolTime_lines_area.fig",Animal)))

%%
figure(1);clf
ax = subtightplot(1,1,1,0.05, 0.1, 0.05);hold on 
scatter(1:36, arrayfun(@(S)S.initscore(1), RDThreadStats),64,'filled')
scatter(0.3+[1:36], arrayfun(@(S)S.initscore(2), RDThreadStats),64,'filled')
scatter(1:36, arrayfun(@(S)S.peakscore(1), RDThreadStats),64,'filled')
scatter(0.3+[1:36], arrayfun(@(S)S.peakscore(2), RDThreadStats),64,'filled')
line([1:36]+[0,0,0.3,0.3]', [arrayfun(@(S)S.initscore(1), RDThreadStats);...
                             arrayfun(@(S)S.peakscore(1), RDThreadStats);...
                             arrayfun(@(S)S.peakscore(2), RDThreadStats);...
                             arrayfun(@(S)S.initscore(2), RDThreadStats)],'color', 'k')
legend(["init score 1", "init score 2", "peak score 1", "peak score 2"])
ylabel("Firing Rate");xlabel("Exp id")

%% Figure: Directly compare the firing rate in 4096 and 50D evolution
text_t_arr = RDEvolCmpTable.text_cmp_t;%arrayfun(@(R)R.text_cmp_t,RDEvolCmpStats);
init_t_arr = RDEvolCmpTable.init_cmp_t;%arrayfun(@(R)R.init_cmp_t,RDEvolCmpStats);
peak_t_arr = RDEvolCmpTable.peak_cmp_t;%arrayfun(@(R)R.peak_cmp_t,RDEvolCmpStats);
xcoord_arr = nan(size(V1msk)); xcoord_arr(V1msk) = 1; xcoord_arr(V4msk) = 2; xcoord_arr(ITmsk) = 3;
figure(4);clf;set(4,'position',[680   315   531   663])
subtightplot(1,1,1, 0.1, 0.08, 0.14);hold on
jitter = 0.1*randn(length(xcoord_arr),1);
scatter(xcoord_arr + jitter, text_t_arr, 50, 'LineWidth', 1.5)
scatter(xcoord_arr + jitter, init_t_arr, 50, 'LineWidth', 1.5)
scatter(xcoord_arr + jitter, peak_t_arr, 50, 'LineWidth', 1.5)
tthresh = tinv([0.01,0.99], 162); % 162 is the dof
line([0.5,0.5;3.5,3.5],[tthresh;tthresh])
legend(["t for 1 block (texture)", "t for 2-3 blocks", "t for peak blocks"],'Location',"Best")
ylabel("t(r_{4096D}, r_{50D})");xlabel("Brain Area");xticks(1:3);xticklabels(["V1","V4","IT"])
title(compose(Animal+" Activation Difference Caused by Reduced Dimension\nRate window [51,200]"))
savedir = "E:\OneDrive - Washington University in St. Louis\Evol_ReducDim\summary";
saveas(4, fullfile(savedir,compose("%s_ActivDiff_t_areas.jpg",Animal)))
savefig(4, fullfile(savedir,compose("%s_ActivDiff_t_areas.fig",Animal)))
% subtightplot(1,3,1);hold on
% subtightplot(1,3,2);hold on
% subtightplot(1,3,3);hold on
%%
DAOA_arr = cell2mat(arrayfun(@(R)R.DAOA',RDThreadStats,'Uni',false));
DAOA2_arr = cell2mat(arrayfun(@(R)R.DAOA2',RDThreadStats,'Uni',false));
figure(7);clf
subtightplot(1,2,1);hold on
scatter(DAOA_arr(1,V1msk),DAOA_arr(2,V1msk),64)
scatter(DAOA_arr(1,V4msk),DAOA_arr(2,V4msk),64)
scatter(DAOA_arr(1,ITmsk),DAOA_arr(2,ITmsk),64)
line([0,4],[0,4],"color","k")
axis equal
legend(["V1","V4","IT"])
title("delta Activation over Activation for reduced dimension evolution")
xlabel("Full space evolution");ylabel("Reduced space evolution")

subtightplot(1,2,2);hold on
scatter(DAOA2_arr(1,V1msk),DAOA2_arr(2,V1msk),64)
scatter(DAOA2_arr(1,V4msk),DAOA2_arr(2,V4msk),64)
scatter(DAOA2_arr(1,ITmsk),DAOA2_arr(2,ITmsk),64)
line([0,4],[0,4],"color","k")
axis equal
legend(["V1","V4","IT"])
title("delta Activation over Activation (2) for reduced dimension evolution")
xlabel("Full space evolution");ylabel("Reduced space evolution")
%%
tstat_arr = cell2mat(arrayfun(@(R)R.tstat',RDThreadStats,'Uni',false));
figure(3);clf;hold on
scatter(tstat_arr(1,V1msk),tstat_arr(2,V1msk),64)
scatter(tstat_arr(1,V4msk),tstat_arr(2,V4msk),64)
scatter(tstat_arr(1,ITmsk),tstat_arr(2,ITmsk),64)
line([0,18],[0,18],"color","k")
axis equal
legend(["V1","V4","IT"])
title("t-stats for peak evolution over reduced evolution")
xlabel("Full space evolution");ylabel("Reduced space evolution")
%
peakS_arr = cell2mat(arrayfun(@(R)R.peakscore',RDThreadStats,'Uni',false));
figure(6);clf;hold on
scatter(peakS_arr(1,V1msk),peakS_arr(2,V1msk),64)
scatter(peakS_arr(1,V4msk),peakS_arr(2,V4msk),64)
scatter(peakS_arr(1,ITmsk),peakS_arr(2,ITmsk),64)
line([0,18],[0,18],"color","k")
axis equal
legend(["V1","V4","IT"])
title("peak activation for reduced dimension evolution")
xlabel("Full space evolution");ylabel("Reduced space evolution")
%%
%% (Beto obsolete version)
% initacts_sem = arrayfun(@(R)cellfun(@sem,R.initacts), RDThreadStats, 'Uni', false); %initacts_sem{6}(1,2)=nan;
% initacts_sem = cell2mat(initacts_sem');
% peakacts_sem = arrayfun(@(R)cellfun(@sem,R.peakacts), RDThreadStats, 'Uni', false); %peakacts_sem{6}(1,2)=nan;
% peakacts_sem = cell2mat(peakacts_sem');
% textacts_sem = arrayfun(@(R)cellfun(@sem,R.textacts), RDThreadStats, 'Uni', false); %textacts_sem{6}(1,2)=nan;
% textacts_sem = cell2mat(textacts_sem');
% initacts_m = arrayfun(@(R)cellfun(@mean,R.initacts), RDThreadStats, 'Uni', false); %initacts_m{6}(1,2)=nan;
% initacts_m = cell2mat(initacts_m');
% peakacts_m = arrayfun(@(R)cellfun(@mean,R.peakacts), RDThreadStats, 'Uni', false); %peakacts_m{6}(1,2)=nan;
% peakacts_m = cell2mat(peakacts_m');
% textacts_m = arrayfun(@(R)cellfun(@mean,R.textacts), RDThreadStats, 'Uni', false); %textacts_m{6}(1,2)=nan;
% textacts_m = cell2mat(textacts_m');
%% Beto obsolete version
% V1msk = prefchan_arr <=48 & prefchan_arr>=33 & (1:36 ~= 4);
% V4msk = prefchan_arr <=64 & prefchan_arr>=49 & (1:36 ~= 4);
% ITmsk = prefchan_arr <=32 & prefchan_arr>=1 & (1:36 ~= 4);
% msk.V1 = V1msk;msk.V4 = V4msk; msk.IT = ITmsk; msk.all = ones(1,length(RDStats),'logical');
