%% Make final population summary figures out of the Kent fittings
% both for prefer driver channels and for non-driver channels
% The plot functions are really WELL written! 
%% Full population analysis
Animal="Both";Set_Path;
global sumdir
sumdir = "O:\Manif_Fitting\summary";
tabdir = "O:\Manif_Fitting\Kent_summary";
mat_dir = "O:\Mat_Statistics";
%% Current version, Get the fitting statistics for all channels with baseline.
poptabdir = "O:\Manif_Fitting\popstats";
poptab = readtable(fullfile(poptabdir,"Both_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
%% Create the masks and find the drivers 
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = poptab.area == "V1";
V4msk = poptab.area == "V4";
ITmsk = poptab.area == "IT";
drivermsk = poptab.is_driver;
prefchmsk = poptab.chan==poptab.prefchan;
Fsigmsk = poptab.F_P<0.001;
R2msk = poptab.R2>0.5;
%%
% msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = hist_plot(poptab, "R2", {validmsk & drivermsk}, ["driver"], ...
            "All Exp (driver)", "bslfit_all", "count");
        %%
h = hist_plot(poptab, "R2", {validmsk & drivermsk, validmsk & drivermsk & poptab.F_P>0.001}, ["driver","driver F>0.001"], ...
            "All Exp (driver)", "bslfit_all_F_P", "count");
h = hist_plot(poptab, "R2", {validmsk & drivermsk & poptab.F_P<0.001}, ["driver F<0.001"], ...
    "All Exp (driver)", "bslfit_all_Fsig", "count");
h = hist_plot(poptab, "R2", {validmsk & drivermsk, validmsk & drivermsk & poptab.F_P>0.001}, ["driver","driver F>0.001"], ...
            "All Exp (driver)", "bslfit_all_F_P", "count");
%%
msk = validmsk & drivermsk & poptab.F_P<0.001;
h = hist_plot(poptab, "R2", {msk&Alfamsk,msk&Betomsk}, ["Alfa","Beto"], ...
    "All Exp (driver F<0.001)", "bslfit_anim_sep_Fsig", "count");
%%
%% Hist comparison of Kappa value
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = hist_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "bslfit_all", "count");
h = hist_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "bslfit_area_sep", "count");
h = hist_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "bslfit_anim_sep", "count");
%%
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = stripe_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "bslfit_all");
h = stripe_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "bslfit_anim_sep");
h = stripe_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "bslfit_area_sep",{[3,1],[2,1],[3,2]});
%%
h = hist_plot(poptab, "R2", {nondrvmsk}, ["non driver"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&Alfamsk, nondrvmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_anim_sep", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&V1msk, nondrvmsk&V4msk, nondrvmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_area_sep", "count",{[3,1],[2,1],[3,2]});
%% Plot the center of tunings
msk = validmsk & drivermsk & poptab.R2 > 0.4 & poptab.space==1;
xscatter_plot(poptab,"phi","theta",{msk},["All"],"All Exp (Driver, R2>0.4, PC23)", "bslfit_PC23all")
xscatter_plot(poptab,"phi","theta",{msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    "All Exp (driver, R2>0.4, PC23)", "bslfit_PC23area_sep")
xscatter_plot(poptab,"phi","theta",{msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"],...
    "All Exp (driver, R2>0.4, PC23)", "bslfit_PC23anim_sep")
%% Plot the kappa and beta (shape parameters)
msk = validmsk & drivermsk & poptab.F_P < 0.001;
% goodfitmsk = poptab.R2 > 0.5;
xscatter_plot(poptab,"kappa","beta",{msk},["All"],"All Exp (Driver P<0.001)", "bslfit_all_Fsig")
xscatter_plot(poptab,"kappa","beta",{msk,msk&poptab.R2>0.5},["All","R>0.5"],"All Exp (Driver P<0.001)", "bslfit_all_F_Fsig")
xscatter_plot(poptab,"kappa","beta",{msk&poptab.R2<0.5,msk&poptab.R2>0.5},["R<0.5","R>0.5"],"All Exp (Driver P<0.001)", "bslfit_all_Fsep")
xscatter_plot(poptab,"kappa","beta",{msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    "All Exp (Driver P<0.001)", "bslfit_area_sep_Fsig")
xscatter_plot(poptab,"kappa","beta",{msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"],...
    "All Exp (Driver P<0.001)", "bslfit_anim_sep_Fsig")

%% Test the ratio between the two statis 
testRatio(poptab,"kappa","beta",{msk},["All"]);
testRatio(poptab,"kappa","beta",{msk,msk&poptab.R2>0.5},["All","R>0.5"]);
testRatio(poptab,"kappa","beta",{msk&poptab.R2<0.5,msk&poptab.R2>0.5},["R<0.5","R>0.5"]);

%% Plot  Kappa distribution
msk = validmsk & drivermsk & poptab.R2 > 0.5;
h = stripe_minor_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (driver, R2>0.5)", "bslfit_area_anim", {[3,1],[2,1],[3,2]});

h = stripe_plot(poptab, "kappa", {msk&V1msk&Alfamsk, msk&V4msk&Alfamsk, msk&ITmsk&Alfamsk}, ["V1","V4","IT"],...
    "All Alfa Exp (driver, R2>0.5)", "bslfit_area_Alfa", {[3,1],[2,1],[3,2]});
h = stripe_plot(poptab, "kappa", {msk&V1msk&Betomsk, msk&V4msk&Betomsk, msk&ITmsk&Betomsk}, ["V1","V4","IT"],...
    "All Beto Exp (driver, R2>0.5)", "bslfit_area_Beto", {[3,1],[2,1],[3,2]});

%% Fitting Validity across region. 
msk = drivermsk;
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],{Alfamsk,Betomsk},["Alfa","Beto"],...
    "All Exp (driver)", "bslfit_anim_area", {[3,1],[2,1],[3,2]});
%
msk = drivermsk & Fsigmsk;
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],{Alfamsk,Betomsk},["Alfa","Beto"],...
    "All Exp (driver, ANOVA p<0.001)", "bslfit_anim_area_Fsig", {[3,1],[2,1],[3,2]});

%% Fitting Validity in 2 monkeys separately.
msk = drivermsk & Fsigmsk & Alfamsk;
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],{msk},["Alfa"],...
    "All Exp (Alfa, driver, ANOVA p<0.001)", "bslfit_Alfa_area_Fsig", {[3,1],[2,1],[3,2]});
msk = drivermsk & Fsigmsk & Betomsk;
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],{msk},["Beto"],...
    "All Exp (Beto, driver, ANOVA p<0.001)", "bslfit_Beto_area_Fsig", {[3,1],[2,1],[3,2]});
%%
msk = drivermsk & Fsigmsk ;
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],{msk},["Both"],...
    "All Exp (Both, driver, ANOVA p<0.001)", "bslfit_Both_area_Fsig", {[3,1],[2,1],[3,2]});
%%
msk = drivermsk & Fsigmsk & (poptab.imgsize>1.1);
h = stripe_minor_plot(poptab, "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], {msk},["Both"],...
    "All Exp (driver, ANOVA p<0.001, imsize>1.0)", "bslfit_Both_nosmall_area_Fsig", {[3,1],[2,1],[3,2]});
%%
msk = drivermsk & Fsigmsk;
h = xscatter_plot(poptab, "imgsize", "R2", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
    "All Exp (driver, ANOVA p<0.001)", "bslfit_area_Fsig")

%% Test if image size or position have an correlation on the kappa or R2.
msk = drivermsk & Fsigmsk;
[ctmp,ptmp] = corr(poptab.imgsize(msk),poptab.kappa(msk))
[ctmp,ptmp] = corr(poptab.imgsize(msk),poptab.R2(msk))
%%
msk = drivermsk & Fsigmsk & Alfamsk;
[ctmp,ptmp] = corr(poptab.imgsize(msk),poptab.R2(msk))
%%
msk = validmsk & drivermsk & poptab.R2 > 0.5 & (poptab.imgsize>1.05) ;
h = stripe_minor_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (driver, R2>0.5, imgsize>1.0)", "bslfit_area_anim_nosmall", {[3,1],[2,1],[3,2]});
%%
msk = validmsk & drivermsk & poptab.R2 > 0.5 & (poptab.imgsize>1.05) ;
h = stripe_minor_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (driver, R2>0.5, imgsize>1.0)", "bslfit_area_anim_nosmall", {[3,1],[2,1],[3,2]});
%%
msk = validmsk & drivermsk & poptab.R2 > 0.5 & (poptab.imgsize>1.05) & Alfamsk;
h = stripe_minor_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    {msk&Alfamsk}, ["Alfa"], "All Exp (driver, R2>0.5, imgsize>1.0)", "bslfit_area_Alfa_nosmall", {[3,1],[2,1],[3,2]});
%%
RTtab = readtable(fullfile(summarydir,"Both"+"_RadialTuningStatsTab_squ.csv"));

%% Test for progression for non-driving channels.
popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0;
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Alfamsk,popmsk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Both_allchan", {[3,1],[2,1],[3,2]});
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Alfamsk}, ["Alfa"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Alfa_allchan", {[3,1],[2,1],[3,2]});
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Betomsk}, ["Beto"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Beto_allchan", {[3,1],[2,1],[3,2]});
%%
popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk;
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Alfamsk,popmsk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Both_nondriver", {[3,1],[2,1],[3,2]});
popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk & Alfamsk;
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Alfamsk}, ["Alfa"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Alfa_nondriver", {[3,1],[2,1],[3,2]});
popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk & Betomsk;
h = stripe_minor_plot(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"],...
    {popmsk&Betomsk}, ["Beto"], "All Exp (R2>0.5, ANOVA P<0.00001)", "bslfit_area_Beto_nondriver", {[3,1],[2,1],[3,2]});

%% Print string for test parameter progression: kappa, beta, R2. across 3 areas in Alfa, Beto, or Both.
diary(fullfile(sumdir,"popul_area_prog.log"))
popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk;
testProgression(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "R2", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "beta", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");

popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk & Alfamsk;
testProgression(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Alfa Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "R2", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Alfa Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "beta", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Alfa Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");

popmsk = poptab.F_P < 1E-5 & poptab.R2 > 0.5 & poptab.kappa > 0 & ~drivermsk & Betomsk;
testProgression(poptab, "kappa", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Beto Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "R2", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Beto Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
testProgression(poptab, "beta", {popmsk&V1msk, popmsk&V4msk, popmsk&ITmsk}, ["V1","V4","IT"], "area", ...
    "Beto Monk All Exp (R2>0.5, ANOVA P<0.00001, non-driver)");
diary off

%% Paired kappa comparison for SU and MU
msk1 = []; msk2 = [];
for idx = find(poptab.unitnum==2)' % collect the pairs
    if poptab.F_P(idx)<1E-5 && poptab.F_P(idx-1)<1E-5 &&...
       poptab.R2(idx)>0.5 && poptab.R2(idx - 1)>0.5 &&...
       poptab.kappa(idx)>0 && poptab.kappa(idx - 1)>0 &&...
       poptab.chan(idx) == poptab.chan(idx - 1) &&...
       poptab.Expi(idx) == poptab.Expi(idx - 1) &&...
       strcmp(poptab.Animal(idx), poptab.Animal(idx - 1))    
    msk1 = [msk1, idx-1];
    msk2 = [msk2, idx];
    end
end
[~,P,CI,TST] = ttest(abs(poptab.kappa(msk1)),abs(poptab.kappa(msk2)))
h = stripe_paired_plot(poptab, "kappa", {msk1, msk2}, ["SU","MU"], "All Chan, both ANOVA P<1E-5, R2>0.5, kappa>0", "SU-MU");
%%
msk1_dr = []; msk2_dr = []; % unit 1 unit 2 pairs that either is driver. 
for idx = find(poptab.unitnum==2)'
    if poptab.F_P(idx)<1E-5 && poptab.F_P(idx-1)<1E-5 &&...
       poptab.R2(idx)>0.5 && poptab.R2(idx - 1)>0.5 &&...
       (drivermsk(idx-1) || drivermsk(idx)) && ...
       poptab.chan(idx) == poptab.chan(idx - 1) &&...
       poptab.Expi(idx) == poptab.Expi(idx - 1) &&...
       strcmp(poptab.Animal(idx), poptab.Animal(idx - 1))    
    msk1_dr = [msk1_dr, idx-1];
    msk2_dr = [msk2_dr, idx];
    end
end
[~,P,CI,TST] = ttest(abs(poptab.kappa(msk1_dr)),abs(poptab.kappa(msk2_dr)))
%
h = stripe_paired_plot(poptab, "kappa", {msk1_dr, msk2_dr}, ["SU","MU"], "Drivers, both, ANOVA P<1E-5, R2>0.5", "SU-MU_driver");

%% Statistical comparison of driver vs Non driver
msk = validmsk & poptab.R2 >0.5;
h = stripe_minor_plot(poptab, "kappa", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5)", "bslfit_drv-non_cmp", {[2,1]});
h = stripe_minor_plot(poptab, "kappa", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5)", "bslfit_drv-non_cmp", {[2,1]});
%%
msk = validmsk & poptab.R2 >0.5;
h = stripe_minor_plot(poptab, "R2", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5)", "bslfit_drv-non_cmp", {[2,1]});
h = stripe_minor_plot(poptab, "R2", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (R2>0.5)", "bslfit_drv-non_cmp", {[2,1]});
%%
msk = validmsk & poptab.F_P < 1E-3;
h = stripe_plot(poptab, "R2", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    "All Exp (ANOVA P<1E-3)", "bslfit_drv-non_Fsign_cmp", {[2,1]});
%%
h = stripe_minor_plot(poptab, "R2", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (ANOVA P<1E-3)", "bslfit_drv-non_Fsign_cmp_anim_sep", {[2,1]});
h = stripe_minor_plot(poptab, "R2", {msk&drivermsk, msk&~drivermsk}, ["Driver","NonDriver"],...
    {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], "All Exp (ANOVA P<1E-3)", "bslfit_drv-non_Fsign_cmp_area_sep", {[2,1]});
%% 
msk = validmsk & poptab.F_P < 1E-3 & poptab.kappa>0;
testRatio(poptab,"kappa","beta",{msk&poptab.R2<0.5,msk&poptab.R2>0.5},["R<0.5","R>0.5"]);
testRatio(poptab,"kappa","beta",{msk&drivermsk,msk&~drivermsk},["Driver","NonDriver"]);
testRatio(poptab,"kappa","beta",{msk&drivermsk&poptab.R2>0.5,msk&~drivermsk&poptab.R2>0.5},["Driver R>0.5","NonDriver R>0.5"]);
%%
msk = validmsk & poptab.F_P < 1E-3&poptab.R2>0.5& poptab.kappa>0;
xscatter_plot(poptab,"kappa","beta",{msk&drivermsk,msk&~drivermsk},["Driver R>0.5","NonDriver R>0.5"],...
    "All Exp (Driver P<0.001 R2>0.5 kappa pos)", "bslfit_all_FsigRsigKpos_drv-non_cmp")
%%
msk = validmsk & poptab.F_P < 1E-3&poptab.R2>0.5;
xscatter_plot(poptab,"kappa","beta",{msk&drivermsk,msk&~drivermsk},["Driver R>0.5","NonDriver R>0.5"],...
    "All Exp (Driver P<0.001 R2>0.5)", "bslfit_all_FsigRsig_drv-non_cmp")


%% OBSOLETE, 
%% Really old version with old Kstats 
alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
preftab = [];
Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, alfatab];
Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, betotab];
%% Making Masks for ploting 
validmsk = ~((alltab.Animal=="Alfa")&(alltab.Expi==10)); % that exp is weird
drivermsk = zeros(size(validmsk)); % Masks of real driver units
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(alltab.Animal(i))(alltab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = driver_unit == alltab.unit(i);
end
Alfamsk = (alltab.Animal=="Alfa");
Betomsk = (alltab.Animal=="Beto");
V1msk = (alltab.chan<=48 & alltab.chan>=33);
V4msk = (alltab.chan>48);
ITmsk = (alltab.chan<33);
%%
figure;hold on
histogram(alltab.R2(validmsk&Alfamsk),20,'FaceAlpha',0.6)
histogram(alltab.R2(validmsk&Betomsk),20,'FaceAlpha',0.6)
%% Show only the real driver units in the channel 
%  Compare the R2 histogram
hist_plot(alltab, "R2", {validmsk&Alfamsk, validmsk&Betomsk},["Alfa","Beto"],...
    "driver valid","anim_sep","count")
hist_plot(alltab, "R2", {validmsk&V1msk, validmsk&V4msk, validmsk&ITmsk},["V1","V4","IT"],...
    "driver valid","area_sep","count")
%%
hist_plot(alltab, "R2", {validmsk& drivermsk& Alfamsk, validmsk& drivermsk& Betomsk},["Alfa","Beto"],...
    "pure driver valid","drv_anim_sep","count")
hist_plot(alltab, "R2", {validmsk& drivermsk& V1msk, validmsk& drivermsk& V4msk, validmsk& drivermsk& ITmsk},["V1","V4","IT"],...
    "pure driver valid","drv_area_sep","count")
hist_plot(alltab, "R2", {validmsk& drivermsk},["All Driver"],...
    "pure driver valid","drv_all","count")

%% All the channels (Popu) OBSOLETE!!!! baseline not taken into account. See above for newer datafile
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat.csv"));
poptab = [alfatab_pop;betotab_pop];
%%
for i = 1:size(poptab,1)
    poptab.imgsize(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgsize;
    poptab.imgposX(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(1);
    poptab.imgposY(i) = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.imgpos(2);
end
%%
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = (poptab.chan<=48 & poptab.chan>=33);
V4msk = (poptab.chan>48);
ITmsk = (poptab.chan<33);
drivermsk = zeros(size(poptab,1),1); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
prefchmsk = poptab.chan==poptab.prefchan;
%% Hist comparison of Kappa value
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = hist_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all", "count");
h = hist_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep", "count");
h = hist_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep", "count");
        %%
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = stripe_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all");
h = stripe_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep");
h = stripe_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep",{[3,1],[2,1],[3,2]});
        %%
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = stripe_plot(poptab, "beta", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all");
h = stripe_plot(poptab, "beta", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep");
h = stripe_plot(poptab, "beta", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep",{[3,1],[2,1],[3,2]});
%%
nondrvmsk = validmsk & ~prefchmsk & poptab.R2 > 0.4 & poptab.F_P < 1E-5;
h = stripe_plot(poptab, "kappa", {nondrvmsk}, ["non driver"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref");
h = stripe_plot(poptab, "kappa", {nondrvmsk&Alfamsk, nondrvmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref_anim_sep");
h = stripe_plot(poptab, "kappa", {nondrvmsk&V1msk, nondrvmsk&V4msk, nondrvmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref_area_sep",{[3,1],[2,1],[3,2]});
%%
nondrvmsk = validmsk & ~prefchmsk & poptab.F_P < 1E-5;
h = hist_plot(poptab, "R2", {nondrvmsk}, ["non driver"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&Alfamsk, nondrvmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_anim_sep", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&V1msk, nondrvmsk&V4msk, nondrvmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_area_sep", "count",{[3,1],[2,1],[3,2]});
%% Hist comparison of Beta value
h = hist_plot(poptab, "beta", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all", "count");
h = hist_plot(poptab, "beta", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep", "count");
h = hist_plot(poptab, "beta", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep", "count");
%% Scatter of the center location of Kent.
set(groot,'defaultAxesTickLabelInterpreter','tex');
msk = validmsk & drivermsk & poptab.R2 > 0.4;
xscatter_plot(poptab,"phi","theta",{msk},["All"],"All Exp (Driver, R2>0.4)", "all")
xscatter_plot(poptab,"phi","theta",{msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    "All Exp (driver, R2>0.4)", "area_sep")
xscatter_plot(poptab,"phi","theta",{msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"],...
    "All Exp (driver, R2>0.4)", "anim_sep")

msk = validmsk & drivermsk & poptab.R2 > 0.4 & poptab.space==1;
xscatter_plot(poptab,"phi","theta",{msk},["All"],"All Exp (Driver, R2>0.4, PC23)", "PC23all")
xscatter_plot(poptab,"phi","theta",{msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    "All Exp (driver, R2>0.4, PC23)", "PC23area_sep")
xscatter_plot(poptab,"phi","theta",{msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"],...
    "All Exp (driver, R2>0.4, PC23)", "PC23anim_sep")
%%
msk = validmsk & drivermsk;
hist_plot(poptab, "R2", {msk},["All Driver"],...
    "pure driver valid","drv_all","count")
hist_plot(poptab, "R2", {msk& Alfamsk, msk& Betomsk},["Alfa","Beto"],...
    "pure driver valid","drv_anim_sep","count")
hist_plot(poptab, "R2", {msk& V1msk, msk& V4msk, msk& ITmsk},["V1","V4","IT"],...
    "pure driver valid","drv_area_sep","count")
%%
stripe_plot(poptab, "R2", {msk& Alfamsk, msk& Betomsk},["Alfa","Beto"],...
    "pure driver valid","drv_anim_sep")
stripe_plot(poptab, "R2", {msk& V1msk, msk& V4msk, msk& ITmsk},["V1","V4","IT"],...
    "pure driver valid","drv_area_sep")


% h = stripe_minor_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
%     {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (driver, R2>0.5, imgsize>1.0)", "bslfit_drv-non_cmp", {[3,1],[2,1],[3,2]});

% validmsk&drivermsk

%% Plain testing function that will pretty print string for use in paper 
function title_str=testRatio(StatsTab_sum, statname1, statname2, masks, labels)
title_str = "";
for i = 1:numel(masks)
    title_str = title_str +compose("\n%s ",labels(i));
    N_i = sum(masks{i});
%     [~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i}) - 2);
    [~,P,CI,TST] = ttest(log(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i})) - log(2));
    title_str = title_str + compose("%s / %s Log: [%.2f,%.2f] (t=%.2f,%.1e)",statname1,statname2,2*exp(CI(1)),2*exp(CI(2)),TST.tstat,P);
    [~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}) - 2*StatsTab_sum.(statname2)(masks{i}));
    title_str = title_str + compose("\t%s - 2* %s Lin: [%.2f,%.2f] (t=%.2f,%.1e,N=%d)",statname1,statname2,CI(1),CI(2),TST.tstat,P,N_i);
end
disp(title_str);
end

function title_str=testProgression(StatsTab_sum, statname, masks, labels, sepvarnm, titstr)
title_str = compose("Test progression of %s ~ %s, for %s\n",statname, sepvarnm, titstr);
valuevec = [];
predvec = [];
for i = 1:numel(masks)
    value = StatsTab_sum.(statname)(masks{i});
    valuevec = [valuevec; value];
    predvec = [predvec; ones(numel(value),1)*(i-1)];
    title_str = title_str+compose("%s %.3f+-%.3f (n=%d)\t",labels(i),mean(value),sem(value),numel(value));
end
lm = fitlm(predvec, valuevec);
[F_p,F_tbl] = anova1(valuevec,predvec,'off');
Fval = F_tbl{2,5}; F_df = F_tbl{4,3};
anovastr = compose("\nANOVA F=%.3f p=%.1e(df=%d,%d)\n",Fval,F_p,F_tbl{2,3},F_tbl{3,3});
lmstr = compose("Linear Regres %s = %.3f + %s * %.3f \n Intercept %.3f+-%.3f, Slope %.3f+-%.3f\n Slope!=0: T=%.1f P=%.1e\n Fit Rsquare=%.3f",...
                statname, lm.Coefficients.Estimate(1), sepvarnm, lm.Coefficients.Estimate(2),...
                lm.Coefficients.Estimate(1), lm.Coefficients.SE(1), ...
                lm.Coefficients.Estimate(2), lm.Coefficients.SE(2), ...
                lm.Coefficients.tStat(2), lm.Coefficients.pValue(2), ...
                lm.Rsquared.Ordinary);
title_str = title_str + anovastr+lmstr;
% disp(lm)
disp(title_str);%+anovastr+lmstr
end

function h = xscatter_corr_plot(StatsTab_sum, statname1, statname2, masks, labels, titstr, savestr)
% Scatter 2 variables in a table. 
global sumdir
h = figure;clf;hold on;set(h,'pos',[1106         327         560         526]);
% statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
title_str = compose("Scatter of %s - %s for\n %s channels %s",...
    statname1,statname2,join(labels),titstr);
for i = 1:numel(masks)
scatter(StatsTab_sum.(statname1)(masks{i}),StatsTab_sum.(statname2)(masks{i}))
[cval,pval] = corr(StatsTab_sum.(statname1)(masks{i}),StatsTab_sum.(statname2)(masks{i}))
title_str = title_str + compose("\n%s: corr %.3f p=%.1e", labels(i), cval, pval)
end
if any(strcmp(statname1,["phi","theta"])), xlim([-pi/2,pi/2]); end
if any(strcmp(statname2,["phi","theta"])), ylim([-pi/2,pi/2]); 
    pbaspect([1 1 1]);daspect([1 1 1]);end % set aspect ratio
% FStat = anova_cells(statscol);
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks); % Num of points within each musk
xlabel(statname1);ylabel(statname2); %
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % Num of points marked on label
saveas(h,fullfile(sumdir,compose("%s_%s_xscatcorr_%s.png", statname1, statname2, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_xscatcorr_%s.pdf", statname1, statname2, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_xscatcorr_%s.fig", statname1, statname2, savestr)))
end

function h = xscatter_plot(StatsTab_sum, statname1, statname2, masks, labels, titstr, savestr)
% Scatter 2 variables in a table. 
% with special treatments for theta, phi and others. 
global sumdir
h = figure;clf;hold on;set(h,'pos',[1106         327         560         526]);
% statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
scatter(StatsTab_sum.(statname1)(masks{i}),StatsTab_sum.(statname2)(masks{i}))
end
if any(strcmp(statname1,["phi","theta"])), xlim([-pi/2,pi/2]); end
if any(strcmp(statname2,["phi","theta"])), ylim([-pi/2,pi/2]); 
    pbaspect([1 1 1]);daspect([1 1 1]);end % set aspect ratio
% FStat = anova_cells(statscol);
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
title_str = compose("Scatter of %s - %s for\n %s channels %s",...
    statname1,statname2,join(labels),titstr);
if any(strcmp(statname1,["phi","theta"])),
for i = 1:numel(masks)
title_str = title_str +compose("\n%s",labels(i));%+
[~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}));
title_str = title_str + compose(" %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname1,CI(1),CI(2),TST.tstat,P);
[~,P,CI,TST] = ttest(StatsTab_sum.(statname2)(masks{i}));
title_str = title_str + compose(" %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname2,CI(1),CI(2),TST.tstat,P);
end
end
if strcmp(statname1,"kappa") && strcmp(statname2,"beta") % if kappa-beta then draw y=x/2 line
    YLIM = ylim(); XLIM=xlim();
    LB = max(2*YLIM(1),XLIM(1));
    UB = min(2*YLIM(2),XLIM(2));
    plot([LB,UB],[LB,UB]/2,'k-.')
    for i = 1:numel(masks)
    title_str = title_str +compose("\n%s",labels(i));
%     [~,P,CI,TST] = ttest(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i}) - 2);
    [~,P,CI,TST] = ttest(log(StatsTab_sum.(statname1)(masks{i}) ./ StatsTab_sum.(statname2)(masks{i})) - log(2));
    title_str = title_str + compose("%s / %s: [%.2f,%.2f] (t=%.2f,%.1e)",statname1,statname2,2*exp(CI(1)),2*exp(CI(2)),TST.tstat,P);
    end
end
Ns = cellfun(@sum,masks); % Num of points within each musk
xlabel(statname1);ylabel(statname2); %
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % Num of points marked on label
saveas(h,fullfile(sumdir,compose("%s_%s_xscat_%s.png", statname1, statname2, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_xscat_%s.pdf", statname1, statname2, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_xscat_%s.fig", statname1, statname2, savestr)))
end

function h = stripe_paired_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, 
    if numel(masks)==2, Tpairs = {[1,2]}; 
    elseif numel(masks)==3, Tpairs = {[1,2],[2,3],[1,3]}; end
end
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
xjitter = 0.1*randn(numel(statscol{1}),1);
Means = []; Stds = [];
for i = 1:numel(masks)
scatter(i+xjitter, statscol{i})
Means = [Means, mean(statscol{i})];
Stds = [Stds, std(statscol{i})];
end
statsarr = cell2mat(statscol)';
plot([1:numel(masks)]'+xjitter',statsarr,'color',[0,0,0,0.1])
xticks(1:numel(masks));xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@numel,statscol);
ylabel(statname)
title(title_str)
legend(compose("%s(%d) %.3f(%.3f)",labels',Ns',Means',Stds')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_strippaircmp.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_strippaircmp.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_strippaircmp.fig", statname, savestr)))
end

function h = stripe_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
mean_i = mean(statscol{i});
sem_i = sem(statscol{i});
N_i = numel(statscol{i});
legstr = compose("%s %.3f(%.3f) (%d)",labels(i),mean_i,sem_i,N_i);
xjitter = 0.15*randn(numel(statscol{i}),1);
scatter(i+xjitter, statscol{i},'DisplayName',legstr);
% medval = median(StatsTab_sum.(statname)(masks{i}));
% vline(medval,"-.",labels(i)+" "+num2str(medval))
end
xticks(1:numel(masks));xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(statname)
title(title_str)
legend('Location','best') % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = stripe_minor_plot(StatsTab_sum, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs)
% Stripe plot comparing a scaler data w.r.t. 2 variables, signified by masks, labels  and  minormsks, minorlabels
%  masks, labels are plotted in different columns; minormsks, minorlabels are plotted in same columne by different color.
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = []; Stds = []; Means = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    xjitter = 0.15*randn(sum(M),1);
    scatter(i+xjitter, StatsTab_sum.(statname)(M),'Marker',marker(j))
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    Stds = [Stds, std(StatsTab_sum.(statname)(M))];
    Means = [Means, mean(StatsTab_sum.(statname)(M))];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname)
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol); % F stats is across the major variables/ msks
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str)
% Ns = cellfun(@sum,masks);
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend(compose("%s(%d) %.2f (%.2f)",labelcol',Ns',Means',Stds'),'Location','best') % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end

function h = hist_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, normstr, Tpairs)
if nargin<8, Tpairs = {}; end
global sumdir
h = figure;clf;hold on 
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
histogram(StatsTab_sum.(statname)(masks{i}),20,'norm',normstr)
medval = median(StatsTab_sum.(statname)(masks{i}));
vline(medval,"-.",labels(i)+" "+num2str(medval))
end
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(normstr)
xlabel(statname) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.png", statname, savestr, normstr)))
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.pdf", statname, savestr, normstr)))
savefig(h,fullfile(sumdir,compose("%s_%s_cmp_%s.fig", statname, savestr, normstr)))
end

function h = hist_cmp_plot(StatsTab_sum, statname, masks, labels, titstr, namestr, normstr)
global sumdir
h = figure;clf;hold on 
for i = 1:numel(masks)
histogram(StatsTab_sum.D1E_Dpr(masks{i}),20,'norm',normstr)
end
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks);
ylabel(normstr)
xlabel(statname)
title(compose("Comparison of %s for\n %s channels in %s Experiments\n"+...
        "t=%.2f (p=%.1e (%d))",statname,join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
legend(compose("%s(%d)",labels',Ns'))%["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.png", statname, namestr, normstr)))
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.pdf", statname, namestr, normstr)))
savefig(h,fullfile(sumdir,compose("%s_%s_cmp_%s.fig", statname, namestr, normstr)))
end
