%% Manifold paper Supplementary figure compare evolution speed across areas
%  Evol_Converg_Speed. This script does the computation majorly. 
%  the `Evol_Converg_Speed_summary` does the population statistics and visualization. 
%  
Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";
%% Collect score trajectory, do gaussian process regression and put them into a table. 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
tic
EvolTrajStat = repmat(struct(),1,numel(EStats));
for Expi = 1:numel(EStats)
    ui = EStats(Expi).evol.unit_in_pref_chan;
    act_col = cellfun(@(psth)squeeze(mean(psth(ui,51:200,:))), EStats(Expi).evol.psth, 'uni',0);
    bsl_col = cellfun(@(psth)squeeze(mean(psth(ui,1:50,:))), EStats(Expi).evol.psth, 'uni',0);
    gen_col = arrayfun(@(i)i*ones(size(EStats(Expi).evol.psth{i},3),1),1:numel(EStats(Expi).evol.psth),'uni',0);
    gen_vec = cell2mat(gen_col');
    act_vec = cell2mat(reshape(act_col,[],1));
    bsl_vec = cell2mat(reshape(bsl_col,[],1));
    bsl_mean = mean(bsl_vec);
    bsl_std = std(bsl_vec);
    bsl_sem = sem(bsl_vec);
    act_traj_mean = cellfun(@mean, act_col);
    act_traj_std = cellfun(@std, act_col);
    act_traj_sem = cellfun(@sem, act_col);
    % fit gaussian process regression (GPR) model to smooth the evol traj
    gprMdl = fitrgp(gen_vec, act_vec);
    gpr_fit = gprMdl.predict([1:max(gen_vec)-1]');
    init_act = gpr_fit(1);
    fina_act = gpr_fit(end);
    [max_act, maxid] = max(gpr_fit);
    bnd = [1, maxid]; %max(gen_vec)-1];
    step50 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.5 + init_act), bnd);
    step63 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.6321 + init_act), bnd);
    step85 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.85 + init_act), bnd);
    EvolTrajStat(Expi).step50 = step50;
    EvolTrajStat(Expi).step63 = step63;
    EvolTrajStat(Expi).step85 = step85;
    EvolTrajStat(Expi).init_act = init_act;
    EvolTrajStat(Expi).fina_act = fina_act;
    EvolTrajStat(Expi).max_act = max_act;
    EvolTrajStat(Expi).act_vec = act_vec;
    EvolTrajStat(Expi).bsl_vec = bsl_vec;
    EvolTrajStat(Expi).gprMdl = gprMdl;
    EvolTrajStat(Expi).gpr_fit = gpr_fit;
    EvolTrajStat(Expi).act_traj_mean = act_traj_mean;
    EvolTrajStat(Expi).act_traj_sem = act_traj_sem;
    EvolTrajStat(Expi).act_traj_std = act_traj_std;
    [~,T_P,~,TST] = ttest2(act_col{1}, act_col{maxid});
    EvolTrajStat(Expi).t_initmax = TST.tstat;
    EvolTrajStat(Expi).t_p_initmax = T_P;
    [~,T_P,~,TST] = ttest2(act_col{1}, act_col{end-1});
    EvolTrajStat(Expi).t_initend = TST.tstat;
    EvolTrajStat(Expi).t_p_initend = T_P;
    EvolTrajStat(Expi).DAOA_initmax = (max_act - init_act) / init_act;
    EvolTrajStat(Expi).DAOA_initend = (fina_act - init_act) / init_act;
%     figure(4);clf;
%     plot([1:max(gen_vec) - 1]',gpr_fit);hold on 
%     scatter(gen_vec, act_vec);
%     title(Expi)
%     box off
    toc
%     break
end
%%
save(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
%%
end 

%% Transcribe the mat file into a table and save it as csv.
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
StatTab = repmat(struct(),1,numel(EStats));
for Expi = 1:numel(EStats)
	StatTab(Expi).Expi = Expi;
	StatTab(Expi).Animal = Animal;
	StatTab(Expi).pref_chan = EStats(Expi).evol.pref_chan;
	StatTab(Expi).unit_in_pref_chan = EStats(Expi).evol.unit_in_pref_chan;
	StatTab(Expi).imgsize = EStats(Expi).evol.imgsize;
	StatTab(Expi).imgpos = EStats(Expi).evol.imgpos;
    StatTab(Expi).step50 = EvolTrajStat(Expi).step50;
    StatTab(Expi).step63 = EvolTrajStat(Expi).step63;
    StatTab(Expi).step85 = EvolTrajStat(Expi).step85;
    StatTab(Expi).init_act = EvolTrajStat(Expi).init_act;
    StatTab(Expi).fina_act = EvolTrajStat(Expi).fina_act;
    StatTab(Expi).t_initmax = EvolTrajStat(Expi).t_initmax;
	StatTab(Expi).t_p_initmax = EvolTrajStat(Expi).t_p_initmax;
	StatTab(Expi).t_initend = EvolTrajStat(Expi).t_initend;
	StatTab(Expi).t_p_initend = EvolTrajStat(Expi).t_p_initend;
	StatTab(Expi).DAOA_initmax = EvolTrajStat(Expi).DAOA_initmax;
	StatTab(Expi).DAOA_initend = EvolTrajStat(Expi).DAOA_initend;
end
StatTab = struct2table(StatTab);
%%
writetable(StatTab, fullfile(mat_dir, Animal+"_EvolTrajStats.csv"))
end
%% Reload csvs and merge the stats for evolution trajectories. 
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
StatTab = cat(1, tabA, tabB);
writetable(StatTab, fullfile(mat_dir, "Both"+"_EvolTrajStats.csv"));


function [gpr_fit, xlinsp, AOC, gprMdl] = GPRfitting(X, Y)
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AOC = trapz(xlinsp, gpr_fit);  % integrate under the gaussian process fitting curve. 
% Plot curve on a fig
% plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end
