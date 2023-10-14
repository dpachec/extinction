
function [] = plot_TG_map(m1, m2, h, t, tRes, f2sav, clim)

f2t = strsplit(f2sav,'_');

figure(); tcl = tiledlayout(1,3);
nexttile
%imagesc(m1); hold on; axis square;colorbar
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; colorbar
if tRes == 1
    plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
else
    plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
end
set(gca, 'clim', clim)

nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square;hold on; colorbar 
if tRes == 1
    plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
else
    plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
end
set(gca, 'clim', clim)

nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
if tRes == 1
    plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
else
    plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
end
set(gca, 'clim', [-3 3])


axesHandles = findall(0, 'type', 'axes');
%set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 150], 'ylim', [1 150]); 
if tRes == 1
    set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', []); 
elseif tRes == 10
    if strcmp(f2t{3}, 'C')
        set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 17], 'ylim', [1 17]); 
    elseif strcmp(f2t{3}, 'V')
        set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 23], 'ylim', [1 23]); 
    end
end

%colorbar
title (tcl, f2sav, 'Interpreter', 'none')

