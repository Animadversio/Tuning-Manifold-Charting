%% Find channel organization
% Encode the array information! the relative coordinates of each channel.

chan_dist = 400; % micron
% IT Channels. The array is a 36 pin triagular grid. 
% X axis goes posterior to anterior, Y axis goes dorsal to ventral
pin36_IT_XY = zeros(36, 2);
pin36_IT_L = zeros(36);
pin36_IT_XY([1:8, 18], 1) = 200 + [0:400:400*8]; 
pin36_IT_XY([1:8, 18], 2) = 0;
pin36_IT_XY(9:17, 1) = 0:400:400*8;
pin36_IT_XY(9:17, 2) = 400 * sin(pi/3);
pin36_IT_XY(20:28, 1) = 200 + [0:400:400*8];
pin36_IT_XY(20:28, 2) = 2 * 400 * sin(pi/3);
pin36_IT_XY([19, 29:36], 1) = 0:400:400*8;
pin36_IT_XY([19, 29:36], 2) = 3 * 400 * sin(pi/3);

IT_chan2pin = [20:28, 29:35, 2:8, 9:17]; % this numbering is not in order. 
IT_chan_XY = pin36_IT_XY(IT_chan2pin, :);

% V4 Channels 49:64. The V1 V4 array is a 18 pin triagular grid. 
% X axis goes posterior to anterior, Y axis goes dorsal to ventral
pin18_V4_XY = zeros(18, 2);
pin18_V4_L = zeros(18);
pin18_V4_XY(1:5, 1) = [0:400:400*4]; 
pin18_V4_XY(1:5, 2) = 0;
pin18_V4_XY(6:9, 1) = 200 + [0:400:400*3];
pin18_V4_XY(6:9, 2) = 400 * sin(pi/3);
pin18_V4_XY(10:14, 1) = [0:400:400*4];
pin18_V4_XY(10:14, 2) = 2 * 400 * sin(pi/3);
pin18_V4_XY(15:18, 1) = 200 + [0:400:400*3];
pin18_V4_XY(15:18, 2) = 3 * 400 * sin(pi/3);

V4_chan2pin = [2:17];
V4_chan_XY = pin18_V4_XY(V4_chan2pin, :);

% V1 Channels 33:48
% X axis goes posterior to anterior, Y axis goes dorsal to ventral
pin18_V1_XY = zeros(18, 2);
pin18_V1_L = zeros(18);
pin18_V1_XY(1:5, 1) = [0:400:400*4]; 
pin18_V1_XY(1:5, 2) = 0;
pin18_V1_XY(6:9, 1) = 200 + [0:400:400*3];
pin18_V1_XY(6:9, 2) = 400 * sin(pi/3);
pin18_V1_XY(10:14, 1) = [0:400:400*4];
pin18_V1_XY(10:14, 2) = 2 * 400 * sin(pi/3);
pin18_V1_XY(15:18, 1) = 200 + [0:400:400*3];
pin18_V1_XY(15:18, 2) = 3 * 400 * sin(pi/3);

V1_chan2pin = [2:17];
V1_chan_XY = pin18_V1_XY(V1_chan2pin, :);


% %%
% figure(1);clf
% scatter(IT_chan_XY(:,1),IT_chan_XY(:,2), 225, 1:32)
% set(gca,'Ydir','reverse')
% axis equal tight
% xlabel("posterior-anterior (mu)");ylabel("dorsal-ventral (mu)")
% colorbar
% %% Layout for IT channels 
% figure(2);clf
% t = tiledlayout(8,20);
% for i = 1:32
%     init_idx = calc_tile_init_idx(IT_chan_XY(i, 1), IT_chan_XY(i, 2), 8, 20);
%     ax = nexttile(init_idx, [2, 2]);
%     axis image
%     title(num2str(i))
% end
% %%
% figure(2);clf
% t = tiledlayout(8,20);
% for i = 1:32
%     init_idx = calc_tile_init_idx(IT_chan_XY(i, 1), IT_chan_XY(i, 2), 8, 20);
%     ax = nexttile(init_idx, [2, 2]);
%     axis image
%     title(num2str(i))
% end
% %%
% figure(3);clf
% t = tiledlayout(8,20);
% for i = 1:32
%     init_idx = calc_tile_init_idx(IT_chan_XY(i, 1), IT_chan_XY(i, 2), 8, 20);
%     ax = nexttile(init_idx, [2, 2]);
%     axis image
%     title(num2str(i))
%     pause
% end
% %% Compute the index in tiled array from X,Y coordinates
% % calc_tile_init_idx(IT_chan_XY(:, 1), IT_chan_XY(:, 2), 8, 20)
% function idx = calc_tile_init_idx(X, Y, rown, coln)
%     Yi = int32(Y / 400 / sin(pi/3));
%     Xi = int32(X / 200) ;
%     idx = coln * Yi * 2 + Xi + 1 ;
% end
% [centers,radii] = imfindcircles(rgb,[20 25],'ObjectPolarity','dark')
% imshow(rgb)
% h = viscircles(centers,radii);