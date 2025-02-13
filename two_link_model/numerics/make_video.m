function make_video(positions, phis, ts, params, skip, repeats, fname, plotPositions)
    figh = clf;
    hold on
    if nargin < 5
        skip = 1;
    end
    if nargin < 6
        repeats = 1;
    end
    if nargin < 7
        fname = "ani";
    end
    if nargin < 8
        plotPositions = false;
    end
    v = VideoWriter(fname,'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 80;
    open(v);

    xlims = [min(positions(:,1)), max(positions(:,1))] + [-2,2]*4e-6 + params.r * [-1,1];
    ylims = [min(positions(:,2)), max(positions(:,2))] + [-2,2]*4e-6 + params.r * [-1,1];
    box on
    xlabel('$x$ (m)','Interpreter','latex')
    ylabel('$y$ (m)','Interpreter','latex')
    
    for i = 1 : skip : length(ts) * repeats
        try 
            delete(h)
        end
        ind = mod(i,length(ts))+1;
        h = draw_robot(positions(ind,:), phis(ind), ts(ind), params);
        if plotPositions
            h = [h; plot(positions(1:ind,1),positions(1:ind,2),'LineWidth',2,'Color','black')];
            h = [h; scatter(positions(1,1),positions(1,2),80,'Filled','MarkerFaceColor','black')];
        end
        xlim(xlims)
        ylim(ylims)
        drawnow
        f = getframe(figh);
        writeVideo(v,f)
    end
    close(v)
end