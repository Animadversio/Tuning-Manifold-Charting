%% 
% Experimental Code to generate samples in selected PC space from an experiment
%%
G = FC6Generator();
% Input the experimental backup folder containing the mat codes files.
imgN_per_arc = 11;
threadi = 1;
backup_dir = "N:\Stimuli\2021-ProjectPFC\2021-Evolutions\2021-06-30-Caos-01\2021-06-30-10-52-27";
newimg_dir = fullfile(backup_dir,"PC_imgs");
mkdir(newimg_dir)
fprintf("Save new images to folder %s\n", newimg_dir)
fprintf("Loading the codes from experiment folder %s\n", backup_dir)
[codes_all, img_ids, generations] = load_codes_all(backup_dir,threadi);
fprintf("Shape of code (%d, %d)\n", size(codes_all,1), size(codes_all,2))
%%%
final_gen_norms = vecnorm(codes_all(generations==max(generations), :),2,2);
final_gen_norm = mean(final_gen_norms);
fprintf("Average norm of the last generation samples %.2f\n", final_gen_norm)
sphere_norm = final_gen_norm;
fprintf("Set sphere norm to the last generations norm!\n")
%
fprintf("Computing PCs\n")
[PC_vectors,PC_Proj_codes,~] = pca(codes_all,'NumComponents',50);
if PC_Proj_codes(end, 1) < 0  % decide which is the positive direction for PC1
    inv_PC1 = true;
    PC1_sign = -1;
else
    inv_PC1 = false;
    PC1_sign = 1;
end

basis23 = PC_vectors(:, 1:3);
basis4950 = PC_vectors(:, [1,49:50]);
rand_vec2 = randn(2, 4096);
rand_vec2 = rand_vec2 - (rand_vec2 * PC_vectors) * PC_vectors';
rand_vec2 = rand_vec2 ./ vecnorm(rand_vec2,2,2);%np.sqrt((rand_vec2**2).sum(axis=1))[:, np.newaxis]
basis_rnd = [PC_vectors(:, 1), rand_vec2'];

basis_col = {basis23, basis4950, basis_rnd};
space_str = ["PC1, PC2, PC3", "PC1, PC49, PC50", "PC1, Random vector1, Random vector2"];
axes_str = ["PC2", "PC3";"PC49","PC50";"RND1","RND2"];
space_shortstr = ["PC23","PC4950","RND12"];
% %% Spherical interpolation
PC2_ang_step = 180 / (imgN_per_arc-1);
PC3_ang_step = 180 / (imgN_per_arc-1);
for spi = 1:numel(basis_col)
basis = basis_col{spi};
basis = [PC1_sign; 1; 1] .* basis';
fprintf("Generating images on %s sphere (radius = %.1f)\n", space_str(spi), sphere_norm)
img_list = {};
for j = -5:5
    for k = -5:5
        theta = PC2_ang_step * j / 180 * pi;
        phi = PC3_ang_step * k / 180 * pi;
        code_vec = [cos(theta)*cos(phi),...
                    sin(theta)*cos(phi),...
                    sin(phi)] * basis;
        code_vec = code_vec / norm(code_vec) * sphere_norm;
        img = G.visualize(code_vec);
        img_list{end+1} = img;
        imwrite(img, fullfile(newimg_dir, compose("norm_%d_%s_%d_%s_%d.jpg", ...
            sphere_norm, axes_str(spi,1), PC2_ang_step * j, axes_str(spi,2), PC3_ang_step* k)));
    end
end
mtg = imtile(img_list,'GridSize',[11,11], 'BorderSize',4,'ThumbnailSize',[256,256]);
imwrite(mtg,fullfile(backup_dir,compose("Norm%d_%s_Montage.png",sphere_norm,space_shortstr(spi))));
end
save(fullfile(newimg_dir, "PC_vector_data.mat"), ...
    "PC_vectors", "rand_vec2", "sphere_norm", "PC2_ang_step", "PC3_ang_step");

%%%

% % %% Spherical interpolation
% PC2_ang_step = 180 / 10
% PC3_ang_step = 180 / 10
% print("Generating images on PC1, PC49, PC50 sphere (rad = %d)" % sphere_norm)
% PC_nums = [0, 48, 49]  % can tune here to change the selected PC to generate images
% img_list = []
% for j in range(-5, 6):
%     for k in range(-5, 6):
%         theta = PC2_ang_step * j / 180 * np.pi
%         phi = PC3_ang_step * k / 180 * np.pi
%         code_vec = np.array([[PC1_sign* np.cos(theta) * np.cos(phi),
%                                         np.sin(theta) * np.cos(phi),
%                                         np.sin(phi)]]) @ PC_vectors[PC_nums, :]
%         code_vec = code_vec / np.sqrt((code_vec**2).sum()) * sphere_norm
%         % img = generator.visualize(code_vec)
%         % img_list.append(img.copy())
%         img = G.render(code_vec)[0]
%         img_list.append(img.copy())
%         plt.imsave(os.path.join(newimg_dir, "norm_%d_PC%d_%d_PC%d_%d.jpg" % (sphere_norm,
%                                                                             PC_nums[1] + 1, PC2_ang_step * j,
%                                                                             PC_nums[2] + 1, PC3_ang_step * k)), img)
% fig2 = utils.visualize_img_list(img_list)
% 
% % %% Spherical interpolation Random 
% PC2_ang_step = 180 / 10
% PC3_ang_step = 180 / 10
% % sphere_norm = 300
% print("Generating images on PC1, Random vector1, Random vector2 sphere (rad = %d)" % sphere_norm)
% % Random select and orthogonalize the vectors to form the sphere
% rand_vec2 = np.random.randn(2, 4096)
% rand_vec2 = rand_vec2 - (rand_vec2 @ PC_vectors.T) @ PC_vectors
% rand_vec2 = rand_vec2 / np.sqrt((rand_vec2**2).sum(axis=1))[:, np.newaxis]
% vectors = np.concatenate((PC_vectors[0:1, :], rand_vec2))
% img_list = []
% for j in range(-5, 6):
%     for k in range(-5, 6):
%         theta = PC2_ang_step * j / 180 * np.pi
%         phi = PC3_ang_step * k / 180 * np.pi
%         code_vec = np.array([[PC1_sign* np.cos(theta) * np.cos(phi),
%                                         np.sin(theta) * np.cos(phi),
%                                         np.sin(phi)]]) @ vectors
%         code_vec = code_vec / np.sqrt((code_vec**2).sum()) * sphere_norm
%         % img = generator.visualize(code_vec)
%         % img_list.append(img.copy())
%         img = G.render(code_vec)[0]
%         img_list.append(img.copy())
%         plt.imsave(os.path.join(newimg_dir, "norm_%d_RND1_%d_RND2_%d.jpg" % (sphere_norm, PC2_ang_step * j, PC3_ang_step * k)), img)
% fig3 = utils.visualize_img_list(img_list)
% 
% fig1.savefig(os.path.join(backup_dir,"Norm%d_PC12_Montage.png"%sphere_norm))
% fig2.savefig(os.path.join(backup_dir,"Norm%d_PC4950_Montage.png"%sphere_norm))
% fig3.savefig(os.path.join(backup_dir,"Norm%d_RND12_Montage.png"%sphere_norm))
% np.savez(os.path.join(newimg_dir, "PC_vector_data.npz"), PC_vecs=PC_vectors, rand_vec2=rand_vec2,
%          sphere_norm=sphere_norm, PC2_ang_step=PC2_ang_step, PC3_ang_step=PC3_ang_step)
function [codes_all, img_ids, code_geni] = load_codes_all(stim_path, threadi, loadblocks)
% Parameters: 
% stim_path: backuped stimuli path to the evolution experiment 
%            e.g. "N:\Stimuli\2019-12-Evolutions\2020-03-10-Alfa-01\2020-03-10-13-50-57";
% threadi: thread number to load; default as 1. 
%            Note threadi starts from 1, so 2nd thread is threadi=2
% loadblocks: default to be empty, means every block. 
%            "last" means the final block 
%            it can be a list of block id or a single number
%            For numbering convention, it usually starts from block001,
%            unless in some case there is a block000 then starts from 0,
%            offset=-1
%            There is test to make sure the file name you are loading is 
%            "block%03d_thread%03d_code.mat"%( block_k + offset, threadi - 1)
if nargin == 1, threadi = 1; loadblocks=[]; end
if nargin == 2, loadblocks=[]; end
data_fn  = ls(fullfile(stim_path, sprintf("*_thread%03d_code.mat", threadi - 1)));
data_fn = sort(string(data_fn)); 
if contains(data_fn{1}, "block000"), offset = -1; fprintf("Saved codes starts from 000. offset -1\n");
elseif contains(data_fn{1}, "block001"), offset = 0; else error; end
% note older version starts from block000
codes_all = [];
code_geni = [];
img_ids = {};
if isempty(loadblocks)
blocks = 1:length(data_fn);
elseif strcmp(loadblocks,"last")
blocks = length(data_fn);
else
blocks = loadblocks; % blocks is basically the list of blocks to load.
end
for block_k = blocks
D = load(fullfile(stim_path, data_fn{block_k}),"codes", "ids");
assert(contains(data_fn{block_k}, sprintf("block%03d_thread%03d_code.mat", block_k + offset, threadi - 1)),...
    sprintf("Missing code mat file %s", sprintf("block%03d_thread%03d_code.mat", block_k, threadi - 1)))
codes_all = [codes_all; D.codes];
code_geni = [code_geni, block_k * ones(1, size(D.codes, 1))];
img_ids = [img_ids, D.ids];
end
img_ids = string(img_ids);
end