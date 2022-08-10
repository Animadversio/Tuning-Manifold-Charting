%% Manifold Paper Supplementary Figure. Single unit, Multi unit correlation. Compared across areas
% Separate script for calculating SU-MU tuning map similarity in more in
% depth fashion.
sumdir = "O:\CortiDistCorr\summary";
SUMUdir = "O:\Manif_SUHash\summary";
%% Load the tuning map similarity matrix 
Animal = "Both";Set_Path;
mat_dir = "O:\Mat_Statistics";
if strcmp(Animal,"Both") % load stats
A = load(fullfile(mat_dir, "Alfa"+'_CortiDisCorr.mat'), 'CortiDisCorr');
B = load(fullfile(mat_dir, "Beto"+'_CortiDisCorr.mat'), 'CortiDisCorr');
CortiDisCorr = [A.CortiDisCorr, B.CortiDisCorr];
else
load(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr') 
end

%% Collect pairs of SU MU correlation into 
areacorrvec_col = {}; % nExp by 5 cell array, containing vectorized corr coef for all pairs
SUMUcorrvec_col = {}; % nExp by 5 cell array, containing vectorized corr coef for SU-MU
for Expi=1:numel(CortiDisCorr) % 2 monkeys are included
% Fmsk = struct2table(CortiDisCorr(1).FStats).F_P<1E-2;
unit_num_arr = CortiDisCorr(Expi).units.unit_num_arr;
spikeID = CortiDisCorr(Expi).units.spikeID;
Fmsk = struct2table(CortiDisCorr(Expi).FStats).F_P<1E-2;
valmsk = unit_num_arr>0;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
ITmsk =  spikeID < 33;
corrmat = CortiDisCorr(Expi).avgsph_corrmat;
% corrmat = CortiDisCorr(1).res_corrmat;
msk_col = {V1msk, V4msk, ITmsk, true};
label_col = ["V1","V4","IT","all"];
for mi = 1:numel(msk_col)
areamsk = msk_col{mi}; arealab = label_col{mi};
pair_ids = get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&areamsk);
SUMUpair_linidx = sub2ind( size(corrmat), pair_ids(:,1), pair_ids(:,2));
SUMUcorrvec = corrmat(SUMUpair_linidx);
in_area_pair_corrvec = get_submat_value(corrmat, areamsk & valmsk & Fmsk, SUMUpair_linidx);
ttest2corr_print(SUMUcorrvec, in_area_pair_corrvec, compose("%s SU-MU",label_col(mi)), compose("%s val pairs",label_col(mi)));
SUMUcorrvec_col{Expi,mi} = SUMUcorrvec;
areacorrvec_col{Expi,mi} = in_area_pair_corrvec;
end

pair_ids = get_SUMUpair_ids(unit_num_arr',spikeID, Fmsk);
SUMUpair_linidx = sub2ind( size(corrmat), pair_ids(:,1), pair_ids(:,2));
SUMUcorrvec = corrmat(SUMUpair_linidx);
all_area_pair_corrvec = cat(1, get_submat_value(corrmat, V1msk & valmsk & Fmsk, SUMUpair_linidx)...
                             , get_submat_value(corrmat, V4msk & valmsk & Fmsk, SUMUpair_linidx)...
                             , get_submat_value(corrmat, ITmsk & valmsk & Fmsk, SUMUpair_linidx));
ttest2corr_print(SUMUcorrvec, all_area_pair_corrvec, "All SU-MU", "All same area pairs");
SUMUcorrvec_col{Expi,5} = SUMUcorrvec;
areacorrvec_col{Expi,5} = all_area_pair_corrvec;
end

%%
diary(fullfile(SUMUdir,"stat_summary.txt"))
ttest2corr_print(cat(1,SUMUcorrvec_col{:,1}), cat(1,areacorrvec_col{:,1}), "All Exp SU-MU (in V1)", "All Exp same area pairs (in V1)");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,2}), cat(1,areacorrvec_col{:,2}), "All Exp SU-MU (in V4)", "All Exp same area pairs (in V4)");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,3}), cat(1,areacorrvec_col{:,3}), "All Exp SU-MU (in IT)", "All Exp same area pairs (in IT)");
% ttest2corr_print(cat(1,SUMUcorrvec_col{:,4}), cat(1,areacorrvec_col{:,4}), "All Exp SU-MU", "All Exp all pairs");
ttest2corr_print(cat(1,SUMUcorrvec_col{:,5}), cat(1,areacorrvec_col{:,5}), "All Exp SU-MU", "All Exp same area pairs");
diary off
label_col = ["V1","V4","IT","all","all_same_array"];
save(fullfile(sumdir,Animal+"SUMU_MapCorr_cmp.mat"), "SUMUcorrvec_col", "areacorrvec_col", "label_col");
save(fullfile(matdir,Animal+"SUMU_MapCorr_cmp.mat"), "SUMUcorrvec_col", "areacorrvec_col", "label_col");

%% Visualize and Summary
h=stripe_simple_plot({cat(1,SUMUcorrvec_col{:,5}), cat(1,areacorrvec_col{:,5})},...
        "corr",["SUMU","All other"], sumdir, "SUMU_bsl_cmp", "markeredgealpha", 0.2);
%cat(1,SUMUcorrvec_col{:,1}), cat(1,areacorrvec_col{:,1})

%%
h=figure;set(h,'pos',[1000  412   484   555]);
cc_col = {cat(1,SUMUcorrvec_col{1:46,5}),...
          cat(1,areacorrvec_col{1:46,5}),...
          cat(1,SUMUcorrvec_col{47:end,5}),...
          cat(1,areacorrvec_col{47:end,5})};
label_arr = ["SU-MU A", "All A", "SU-MU B", "All B"];
violinplot_cell(cc_col, label_arr,'showData',true)
ylabel("Correlation")
saveallform(sumdir,"SUMU_other_cmp_violin_scatter",h)

%%
h=figure;set(h,'pos',[1000  412   484   555]);
violinplot_cell(cc_col, label_arr,'showData',false)
ylabel("Correlation")
saveallform(sumdir,"SUMU_other_cmp_violin_pure",h)
%%
cc_col = {cat(1,SUMUcorrvec_col{1:46,1}),...
          cat(1,areacorrvec_col{1:46,1}),...
          cat(1,SUMUcorrvec_col{47:end,1}),...
          cat(1,areacorrvec_col{47:end,1}),...
          cat(1,SUMUcorrvec_col{1:46,2}),...
          cat(1,areacorrvec_col{1:46,2}),...
          cat(1,SUMUcorrvec_col{47:end,2}),...
          cat(1,areacorrvec_col{47:end,2}),...
          cat(1,SUMUcorrvec_col{1:46,3}),...
          cat(1,areacorrvec_col{1:46,3}),...
          cat(1,SUMUcorrvec_col{47:end,3}),...
          cat(1,areacorrvec_col{47:end,3})};
label_arr = ["SU-MU V1 A", "All V1 A", "SU-MU V1 B", "All V1 B",...
             "SU-MU V4 A", "All V4 A", "SU-MU V4 B", "All V4 B",...
             "SU-MU IT A", "All IT A", "SU-MU IT B", "All IT B"];
h=figure;set(h,'pos',[1000  412   684   555]);
violinplot_cell(cc_col, label_arr,'showData',true,'GroupOrder',cellstr(label_arr))
ylabel("Correlation")
xtickangle(30)
saveallform(sumdir,"SUMU_other_area_anim_sep_cmp_violin_scatter",h)
clf(h)
violinplot_cell(cc_col, label_arr,'showData',false,'GroupOrder',cellstr(label_arr))
ylabel("Correlation")
xtickangle(30)
saveallform(sumdir,"SUMU_other_area_anim_sep_cmp_violin_pure",h)
%%
cc_col = {cat(1,SUMUcorrvec_col{:,1}),...
          cat(1,areacorrvec_col{:,1}),...
          cat(1,SUMUcorrvec_col{:,2}),...
          cat(1,areacorrvec_col{:,2}),...
          cat(1,SUMUcorrvec_col{:,3}),...
          cat(1,areacorrvec_col{:,3})};
label_arr = ["SU-MU V1", "All V1",...
             "SU-MU V4", "All V4",...
             "SU-MU IT", "All IT"];
h=figure;set(h,'pos',[1000         370         515         597]);
violinplot_cell(cc_col, label_arr,'showData',true,'GroupOrder',cellstr(label_arr))
ylabel("Correlation")
saveallform(sumdir,"SUMU_other_area_sep_cmp_violin_scatter",h)
clf(h)
violinplot_cell(cc_col, label_arr,'showData',false,'GroupOrder',cellstr(label_arr))
ylabel("Correlation")
saveallform(sumdir,"SUMU_other_area_sep_cmp_violin_pure",h)


%%
Expi = 11;
unit_num_arr = CortiDisCorr(Expi).units.unit_num_arr;
spikeID = CortiDisCorr(Expi).units.spikeID;
Fmsk = struct2table(CortiDisCorr(Expi).FStats).F_P<1E-2;
valmsk = unit_num_arr>0;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
ITmsk =  spikeID < 33;
get_SUMUpair_ids(unit_num_arr', spikeID, Fmsk&V4msk)
%% Tune map pair comparison 
Animal="Alfa"; Expi = 15;
unit_num_arr = CortiDisCorr(Expi).units.unit_num_arr;
spikeID = CortiDisCorr(Expi).units.spikeID;
Fmsk = struct2table(CortiDisCorr(Expi).FStats).F_P<1E-2;
valmsk = unit_num_arr>0;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
ITmsk =  spikeID < 33;
pairs = get_SUMUpair_ids(unit_num_arr', spikeID, valmsk&Fmsk);%&ITmsk
for rowi = 1:size(pairs,1)
Manif_Map_show_fun(MapVarStats, "Alfa", Expi, pairs(rowi,:))
pause
end
%%
saveallform(SUMUdir, compose("SUMUpair_cmp_%s_Exp%02d_%s",Animal,Expi,strrep(num2str(pairs(rowi,:)),"  ","_")))


function corrvec = get_submat_value(corrmat, msk, exclude_linidx)
linidxmat = reshape(1:numel(corrmat),size(corrmat)); % matrix formed by linear index
linidx_submat = linidxmat(msk,msk);
linidx_submat_tril = tril(linidx_submat,-1);
linidx_final = nonzeros(linidx_submat_tril);
if nargin > 2
    linidx_final = setdiff(linidx_final, exclude_linidx);
end
corrvec = corrmat(linidx_final); % use these linear indices to fetch values.
end

function pairs = get_SUMUpair_ids(unit_num_arr,chan_num_arr,msk)
U1idxs = strfind(unit_num_arr,[1,2]);
pairs = [U1idxs',U1idxs'+1];
chan_ids = chan_num_arr(pairs);
assert(all(chan_num_arr(U1idxs')==chan_num_arr(U1idxs'+1)),"the spike id doesn't match ")
if nargin >2
validunit = msk(pairs);
if size(pairs,1)==1, validunit=validunit'; end
valid_row = find(all(validunit,2));
pairs = pairs(valid_row,:);
assert(all(validunit(valid_row,:),'all'),"the spike id doesn't match ")
end
end

function h = violinplot_cell(value_col, label_arr, varargin)
Yvec = [];
labvec = [];
for i=1:numel(value_col)
    Yvec = cat(1,Yvec, reshape(value_col{i},[],1));
    labvec = cat(1, labvec, repmat(label_arr(i),numel(value_col{i}),1));
end
violinplot(Yvec,labvec,varargin{:})
end