
function [] = plot_TG_map(m1, m2, h, t, f2sav, clim)

f2t = strsplit(f2sav,'_');
tRes  = strsplit(f2t {7}, '-'); slidTW  = double(string(tRes {2}));
winSize = double(string(tRes{1}));


figure(); tcl = tiledlayout(1,3);
nexttile
%imagesc(m1); hold on; axis square;colorbar
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; colorbar
if winSize == 50
    if slidTW  == 1
        plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
    else
        plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
elseif winSize == 20
    if slidTW  == 1
        plot(get(gca,'xlim'), [10 10],'k', 'linewidth', 1); plot([10 10], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
end
set(gca, 'clim', clim)

nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square;hold on; colorbar 
if winSize == 50
    if slidTW  == 1
        plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
    else
        plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
elseif winSize == 20
    if slidTW  == 1
        plot(get(gca,'xlim'), [10 10],'k', 'linewidth', 1); plot([10 10], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
end
set(gca, 'clim', clim)

nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
if winSize == 50
    if slidTW  == 1
        plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
    else
        plot(get(gca,'xlim'), [3 3],'k', 'linewidth', 1); plot([3 3], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
elseif winSize == 20
    if slidTW  == 1
        plot(get(gca,'xlim'), [10 10],'k', 'linewidth', 1); plot([10 10], get(gca,'ylim'),'k', 'linewidth', 1); 
    end
end
set(gca, 'clim', [-4 4])


axesHandles = findall(0, 'type', 'axes');

if winSize == 50
    if slidTW  == 1
        if strcmp(f2t{3}, 'C')
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 125], 'ylim', [1 125]); 
            %set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 200], 'ylim', [1 200]); 
        elseif strcmp(f2t{3}, 'V')
            %set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 125], 'ylim', [1 125]); 
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 225], 'ylim', [1 225]); 
        end
    elseif slidTW  == 10
        if strcmp(f2t{3}, 'C')
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 17], 'ylim', [1 17]); 
        elseif strcmp(f2t{3}, 'V')
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 23], 'ylim', [1 23]); 
        end
    end
elseif winSize == 20
    if slidTW  == 1
        if strcmp(f2t{3}, 'C')
            %set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 125], 'ylim', [1 125]); 
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 185], 'ylim', [1 185]); 
        elseif strcmp(f2t{3}, 'V')
            %set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 125], 'ylim', [1 125]); 
            set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 200], 'ylim', [1 200]); 
        end
    end
end

 

%colorbar
title (tcl, f2sav, 'Interpreter', 'none')

