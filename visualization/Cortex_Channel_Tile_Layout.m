function [t,ax_arr] = Cortex_Channel_Tile_Layout(arr_name, figh)
% This is obsolete partial file. Use Cortex_Channel_Tile_Layout_All
% instead!
Extract_Channel_Organization;
ax_arr = {};
if contains(arr_name,"IT")
figure(figh);set(figh,'Visible','off')
t = tiledlayout(8,20,'TileSpacing','Compact');
for i = 1:32
    init_idx = calc_tile_init_idx(IT_chan_XY(i, 1), IT_chan_XY(i, 2), 8, 20);
    ax = nexttile(init_idx, [2, 2]);
    axis image
    %title(num2str(i))
    ax_arr{end+1} = ax;
end
end
end
function idx = calc_tile_init_idx(X, Y, rown, coln)
    Yi = int32(Y / 400 / sin(pi/3));
    Xi = int32(X / 200) ;
    idx = coln * Yi * 2 + Xi + 1 ;
end