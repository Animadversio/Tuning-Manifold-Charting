function axs = AlignAxisCLimits(axs)
% Given a group of axis, make their Color limit the same as each other
    if iscell(axs), axs = [axs{:}]; end
    CLIM = caxis(axs(1));
    for i = 2:numel(axs)
        CLIM_tmp = caxis(axs(i));
        CLIM = [min(CLIM(1), CLIM_tmp(1)), max(CLIM(2), CLIM_tmp(2))];
    end
    for i = 1:numel(axs)
        caxis(axs(i), CLIM);
    end
end