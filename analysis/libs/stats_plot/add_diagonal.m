function add_diagonal(ax,varargin)
XLIM = xlim(gca);
YLIM = ylim(gca);
MIN = max(XLIM(1),YLIM(1));
MAX = min(XLIM(2),YLIM(2));hold on
plot([MIN,MAX],[MIN,MAX],varargin{:})%,'handlevisi',false
hold off
end