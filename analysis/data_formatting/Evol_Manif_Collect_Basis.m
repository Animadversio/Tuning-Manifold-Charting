%% Collect latent codes
Set_Path;

for Animal=["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats')
%%
basisStats = repmat(struct(), 1, numel(Stats));
%% Load the source basis data for Manifold Experiments from python
% basis_path = fullfile(Stats(Expi).meta.stimuli,"PC_vector_data.npz");
for Expi = 1:numel(Stats)
basis_path = fullfile(EStats(Expi).meta.stimuli,"PC_imgs","PC_vector_data.npz");
f = py.numpy.load(basis_path);
PC_Vec = f.get('PC_vecs').double;
Rand_Vec2 = f.get('rand_vec2').double;
sphere_norm = f.get('sphere_norm').double;
f.close();
% Get the mat containing all the codes of the last generation. To see
% whether we should inverse PC1 
matfns = string(ls(fullfile(EStats(Expi).meta.stimuli,"*.mat")));
mean_codes = get_mean_codes(EStats(Expi).meta.stimuli, matfns);
proj_coord = mean_codes(end,:) * PC_Vec';
if proj_coord(1)>0
    PC1_sign = 1;
else
    fprintf("The evolution direction is inverse to the PC1 direction of PCA. Inverse PC1 as basis\n")
    PC1_sign = -1;% Note the final PC may need to reverse! not always the same dir!
end
basisStats(Expi).basis_23 = [PC1_sign,1,1]' .* PC_Vec(1:3,:);
basisStats(Expi).basis_4950 = [PC1_sign,1,1]' .* [PC_Vec(1,:);PC_Vec(49:50,:)];
basisStats(Expi).basis_RND = [PC1_sign,1,1]' .* [PC_Vec(1,:);Rand_Vec2];
basisStats(Expi).evoltraj = mean_codes;
basisStats(Expi).sphere_norm = sphere_norm;
end
%%
save(fullfile(matdir, Animal+"_Basis_stats.mat"), 'basisStats')
end

function mean_codes = get_mean_codes(stimpath, matfns_sorted)
% Load the mean codes from each blocks in a stimuli folder.
mean_codes = [];
for i = 1:numel(matfns_sorted)
    code_data = load(fullfile(stimpath,matfns_sorted(i)));
    mean_codes = [mean_codes; mean(code_data.codes,1)];
end
assert(size(mean_codes,1)==numel(matfns_sorted))
assert(size(mean_codes,2)==4096)
end