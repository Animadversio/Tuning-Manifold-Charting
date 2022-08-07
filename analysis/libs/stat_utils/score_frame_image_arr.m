function [frame_img_list, Clim] = score_frame_image_arr(img_list, score_mat, clim, cmap, LineWidth)
% Use the cmap and clim to map values in hl_mat to color, and form color frame
% for corresponding image in the img_list. 
% Signature
%   frame_img_list = score_frame_image_arr(img_list, score_mat, clim, cmap, LineWidth)
% 
% Example 
%   img_frame = score_frame_image_arr(code_mean_imgs(:,:,:,1:4:end),actmean(1:4:end));
%   mtg = imtile(img_frame, 'GridSize',[3,7],'Thumb',[296,296]);
%   montage(img_frame)
% 
% Arguments 
% img_list is a image cell array, same shape as score_mat. 
%           Now support image file path (string array) as well, the images
%           will be read and padded. @Feb.4th
%           Now support 4d image tensor assume to be [H,W,C,B], will be
%           converted to 1d image cell array in the parsing, and then
%           reshaped to match the score_mat. 
% score_mat is a score matrix, with nan is images with no observation.
% clim is the limit of value mapped to color Can be automatized
% cmap is a K-by-3 matrix coding the RGB values in 1:64. e.g. `jet`,
%           `parula` will return colormap in these formats s
% LineWidth control the width of padding
%
% Returns
%   frame_img_list: Cell array of framed images
% 
if nargin <= 4 % default parameters
LineWidth = 10;
if nargin <= 3
cmap = parula;
if nargin == 2
clim = [min(score_mat,[],'all'), max(score_mat,[],'all')];
end;end;end

Cmin = clim(1); Cmax = clim(2);
Clim = [Cmin, Cmax];
fprintf("Current color limit [%.2f,%.2f]\n",Cmin,Cmax)
if ndims(img_list)==4 && ~iscell(img_list)
    Nimgs = size(img_list,4);
    img_list = arrayfun(@(imgi)img_list(:,:,:,imgi),1:Nimgs,'uni',0);
    img_list = reshape(img_list, size(score_mat));
end
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
            scale_val = max(0,min(1,scale_val)); % added to clip the value
            c = interp1(cmap, scale_val * (size(cmap, 1) - 1) + 1); % Note, interpolation can be done from 1-64, not from 0
            if isstring(img_list{i,j}) || ischar(img_list{i,j})
                img = imread(img_list{i,j});
                color_scale = 255;
            else
                img = img_list{i,j};
                if max(img,[],'all') > 1.2 || isa(img,"uint8")
                color_scale = 255;
                else
                color_scale = 1;
                end
            end
            pad_img = padarray(img, [2*LineWidth, 2*LineWidth], 0);
            tmp_img = insertShape(pad_img, ...
                'Rectangle', [LineWidth+1,LineWidth+1,...
                            size(img,2)+2*LineWidth,...
                            size(img,1)+2*LineWidth], ...
                'LineWidth', 2 * LineWidth, 'Color', color_scale * c);
            frame_img_list{i,j} = tmp_img;
        end
    end
end
end