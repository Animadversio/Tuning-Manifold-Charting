function axs = AlignAxisLimits(axs)
% Given a group of axis, make their YLim and XLim the same as each other
% Signature: axs = AlignAxisLimits(axs)
%     YLIM = [min([axs(1).YLim, axs(2).YLim]), max([axs(1).YLim, axs(2).YLim])];
%     XLIM = [min([axs(1).XLim, axs(2).XLim]), max([axs(1).XLim, axs(2).XLim])];
    if iscell(axs), axs = [axs{:}]; end
    XLIM = axs(1).XLim; YLIM = axs(1).YLim;
    for i = 2:numel(axs)
        XLIM = [min(XLIM(1), axs(i).XLim(1)), max(XLIM(2), axs(i).XLim(2))];
        YLIM = [min(YLIM(1), axs(i).YLim(1)), max(YLIM(2), axs(i).YLim(2))];
    end
    for i = 1:numel(axs)
        axs(i).XLim = XLIM; axs(i).YLim = YLIM;
    end
end