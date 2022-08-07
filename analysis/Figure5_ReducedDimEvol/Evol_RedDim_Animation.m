% Create Animation for Dual evolution in Reduced Dimension Evolution
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, Animal+"_RDEvol_stats.mat"), 'RDStats')
%% 
ExpType = "RDEvol";
for Expi = 24:length(RDStats)
fprintf("Processing Reduced Dimension Evolution Exp %d\n",Expi)
thread_n = RDStats(Expi).evol.thread_num;
% ui = RDStats(Expi).evol.unit_in_pref_chan(1);
% assert(all(ui==RDStats(Expi).evol.unit_in_pref_chan))
% the following snippet is to get rid of 0 unit (null channel)
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == 1) & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
ui = find(prefchan_id==chid);
Window = 50:200;
imgColl = repmat("", RDStats(Expi).evol.block_n-1,thread_n);
scoreColl = zeros(RDStats(Expi).evol.block_n-1,thread_n);
meanscoreColl = zeros(RDStats(Expi).evol.block_n-1,thread_n);
for thr_i = 1:thread_n
for blocki = 1:RDStats(Expi).evol.block_n-1
    gen_scores = squeeze(mean(RDStats(Expi).evol.psth{thr_i,blocki}(ui,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    imgidx = RDStats(Expi).evol.idx_seq{thr_i,blocki}(maxIdx);
    imgfullfn = ls(fullfile(RDStats(Expi).meta.stimuli, RDStats(Expi).imageName(imgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(RDStats(Expi).meta.stimuli, RDStats(Expi).imageName(imgidx)+"*"))
    imgColl(blocki,thr_i) = fullfile(RDStats(Expi).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki,thr_i) = maxScore;
    meanscoreColl(blocki,thr_i) = mean(gen_scores);
end
end
%% Generation averaged psth and sem
evol_stim_fr = cellfun(@(psth)mean(psth,3),RDStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_fr = cell2mat(reshape(evol_stim_fr',1,1,[],thread_n));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),RDStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_sem = cell2mat(reshape(evol_stim_sem',1,1,[],thread_n));
% Scalor score for evolution
score_avg = cellfun(@(psth)mean(psth(:,51:200,:),'all') - mean(psth(:,1:50,:),'all'),RDStats(Expi).evol.psth);
score_sem = cellfun(@(psth)std(squeeze(mean(psth(:,51:200,:),[1,2])))...
    /sqrt(size(psth,3)),RDStats(Expi).evol.psth); % - mean(psth(:,1:50,:),[1,2])
%% Generate Movies
result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
savepath = result_dir; % fullfile(result_dir, compose("%s_Evol_"))
color_seq = RDStats(Expi).color_seq;
v = VideoWriter(fullfile(savepath,compose('%s_RDEvol_Exp%02d_Best_PSTH',Animal,Expi)));
v.FrameRate = 4;
open(v);
h3=figure(4);set(4,'position',[263         148        1341         735]);clf;
ax1_1 = subtightplot(2,4,1,0.07);
% set(ax1_1,"position",[0.07,0.578,0.439,0.3412]);
title(compose("%s",RDStats(Expi).evol.optim_names(1)))
scoreYLIM = [0,max(score_avg(:,1:end-1)+score_sem(:,1:end-1),[],'all')];
ax3_1 = subtightplot(2,4,2,0.07);
shadedErrorBar([],score_avg(1,1:end-1),score_sem(1,1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
xlabel("Generations");ylabel("Response fr (Hz)");axis tight;ylim(scoreYLIM)
title(RDStats(Expi).evol.optim_names(1)+" Evol Traj");
psthYLIM = [0,max(evol_stim_fr(ui,:,1:end-1,:)+evol_stim_sem(ui,:,1:end-1,:),[],'all')];
ax2_1 = subtightplot(2,2,3,0.07,0.07,0.05);
bgpsth_1 = plot(squeeze(evol_stim_fr(ui,:,1:end-1,1)),'Color',[0.7,0.7,0.7]);hold on
sEB_1 = shadedErrorBar([],evol_stim_fr(ui, :, 1, 1),evol_stim_sem(ui, :, 1, 1),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)");ylim(psthYLIM)
title("Evoked PSTH")
if thread_n==2
ax1_2 = subtightplot(2,4,3,0.07);
% set(ax1_2,"position",[0.07,0.578,0.439,0.3412]);
title(compose("%s",RDStats(Expi).evol.optim_names(2)))
ax3_2 = subtightplot(2,4,4,0.07);
shadedErrorBar([],score_avg(2,1:end-1),score_sem(2,1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
xlabel("Generations");ylabel("Response fr (Hz)");axis tight;ylim(scoreYLIM)
title(RDStats(Expi).evol.optim_names(2)+" Evol Traj");
ax2_2 = subtightplot(2,2,4,0.07,0.07,0.05);
bgpsth_2 = plot(squeeze(evol_stim_fr(ui,:,1:end-1,2)),'Color',[0.7,0.7,0.7]);hold on
sEB_2 = shadedErrorBar([],evol_stim_fr(ui, :, 1, 2),evol_stim_sem(ui, :, 1, 2),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)");ylim(psthYLIM)
title("Evoked PSTH")
end
ST = suptitle(compose("%s Reduced Dimen Evol Exp %02d pref chan %s", ...
    Animal, Expi, RDStats(Expi).units.unit_name_arr(RDStats(Expi).units.pref_chan_id)));
%
for blocki = 1:RDStats(Expi).evol.block_n-1
    set(h3,"CurrentAxes",ax1_1)%subplot(2,2,1);
    imshow(imgColl{blocki,1}); 
    title(compose("Gen%d Best Rate %.1f", blocki, scoreColl(blocki,1))) ; 
    if thread_n==2,set(h3,"CurrentAxes",ax1_2)
    imshow(imgColl{blocki,2}); 
    title(compose("Gen%d BestRate %.1f", blocki, scoreColl(blocki,2))) ; 
    end
%     set(h3,"CurrentAxes",ax2_1)
    ax2_1.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, meanscoreColl(blocki,1));
%     set(h3,"CurrentAxes",ax2_2)
    if thread_n==2
    ax2_2.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, meanscoreColl(blocki,2));
    end
    psthcur = evol_stim_fr(ui, :, blocki, 1);
    psthsemcur = evol_stim_sem(ui, :, blocki, 1);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB_1.mainLine.YData = psthcur;
    sEB_1.edge(1).YData = lE;
    sEB_1.edge(2).YData = uE;
    sEB_1.patch.Vertices(:,2) = [lE,fliplr(uE)];
    if thread_n==2
    psthcur = evol_stim_fr(ui, :, blocki, 2);
    psthsemcur = evol_stim_sem(ui, :, blocki, 2);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB_2.mainLine.YData = psthcur;
    sEB_2.edge(1).YData = lE;
    sEB_2.edge(2).YData = uE;
    sEB_2.patch.Vertices(:,2) = [lE,fliplr(uE)];
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h3);
    writeVideo(v,Fs);
end
close(v);
end