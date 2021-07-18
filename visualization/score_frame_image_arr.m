function frame_img_list = score_frame_image_arr(img_list, score_mat, clim, cmap, LineWidth)
% Use the cmap and clim to map values in hl_mat to color, and form color frame
% for corresponding image in the img_list. 
% 
% Arguments 
% img_list is a image cell array, same shape as score_mat. 
%           Now support image file path (string array) as well, the images
%           will be read and padded. @Feb.4th
% score_mat is a score matrix, with nan is images with no observation.
% clim is the limit of value mapped to color Can be automatized
% cmap is a K-by-3 matrix coding the RGB values in 1:64. e.g. `jet`,
%           `parula` will return colormap in these formats s
% LineWidth control the width of padding
%
% Returns
%   frame_img_list: Cell array of framed images
if nargin <= 4 % default parameters
LineWidth = 10;
if nargin <= 3
cmap = parula;
if nargin == 2
clim = [min(score_mat,[],'all'), max(score_mat,[],'all')];
end;end;end
Cmin = clim(1); Cmax = clim(2);
assert(all(size(img_list) == size(score_mat)), ...
    "Score matrix and image cell array size doesn't matach")
% LineWidth = 50; % Key parameter controlling the width of margin, can be different for different
frame_img_list = cell(size(img_list)); % the list storing the padde image
for j = 1:size(img_list, 2)
    for i = 1:size(img_list, 1)
        if isnan(score_mat(i,j)) % check if there is any score
            frame_img_list{i,j} = [];
        else
            scale_val = (score_mat(i,j) - Cmin) / (Cmax - Cmin);
            c = interp1(cmap, scale_val * (size(cmap, 1) - 1) + 1); % Note, interpolation can be done from 1-64, not from 0
            if isstring(img_list{i,j}) || ischar(img_list{i,j})
                img = imread(img_list{i,j});
            else
                img = img_list{i,j};
            end
            pad_img = padarray(img, [2*LineWidth, 2*LineWidth], 0);
            tmp_img = insertShape(pad_img, ...
                'Rectangle', [LineWidth+1,LineWidth+1,...
                            size(img,2)+2*LineWidth,...
                            size(img,1)+2*LineWidth], ...
                'LineWidth', 2 * LineWidth, 'Color', 256 * c);
            frame_img_list{i,j} = tmp_img;
        end
    end
end
end