function my_mediation_path_diagram(stats)
    
    wh = which('intersectLinePolygon');
    if isempty(wh)
        disp('Warning: To create a mediation path diagram, you must have the external');
        disp('geom2d toolbox on your path.  I can''t find it, so the path diagram will be skipped.');
        return
    end
    
    setup_figure();

    %%MD adaptations
    %what type of plot
    if length(stats.names)>5 % tested more than one mediator.
        disp('NOTE: plotting multi-level version (Mdavidson tampered with...)');
        mlevel=1;
    else
        mlevel=0;
    end
    %% plot markers

    Xname = strrep(stats.inputOptions.vnames{1}, '_', '');
    Mname = strrep(stats.inputOptions.vnames{3}, '_', '');
    Yname = strrep(stats.inputOptions.vnames{2}, '_', '');
    
    % X
    Xloc = [0 0];
    Xpts = plot_region_circ(Xname, Xloc,stats.plotCols(1,:));

    % M
    Mloc = [1.5 1];
    Mpts = plot_region_circ(Mname, Mloc,stats.plotCols(3,:));

    % Y
    Yloc = [3 0];
    Ypts = plot_region_circ(Yname, Yloc, stats.plotCols(2,:));

    %multilevel options
    if mlevel
        M2name =  strrep(stats.inputOptions.vnames{4}, '_', '');
        %M2
        if stats.Cov==1
            M2loc = [1.5 2];
        else
            M2loc = [1.5 -1];
        end
            M2pts = plot_region_circ(M2name, M2loc,stats.plotCols(4,:));
        
    end
    %% path arrows

    %% X -> M
    stats_idx = 1;
    [astart aend] = compute_line_points(Xloc, Xpts, Mloc, Mpts);
    if stats.p(1, stats_idx) < .05
        plot_arrow(astart, aend);
    end
    loc_text = Xloc + (Mloc - Xloc) .* .5 + [.05 .1];
    plot_stats_string(loc_text, stats, stats_idx, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline');
%add path ID
    text(loc_text(1)-.1, loc_text(2)+.1, 'a', 'FontSize', 16,'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    
    %%
    
    % plot a blue line if we have level 2 moderators as well
    if size(stats.p, 1) > 1  
        if stats.p(2, stats_idx) < .05
            hh = arrow(astart+.05, aend+.05, 'LineWidth', 2, 'Color', [.3 .5 1]);
        end
    end
    
    %% M -> Y
    stats_idx = 2;
    [bstart bend] = compute_line_points(Mloc, Mpts, Yloc, Ypts);
    if stats.p(1, stats_idx) < .05
        plot_arrow(bstart, bend);
    end
    loc_text = Mloc + (Yloc - Mloc) .* .5 + [-.05 .1];
    plot_stats_string(loc_text, stats, stats_idx, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');

    %add path ID
    text(loc_text(1)+.1, loc_text(2)+.1, 'b', 'FontSize', 16,'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    
    
    % plot a blue line if we have level 2 moderators as well
    if size(stats.p, 1) > 1  
        if stats.p(2, stats_idx) < .05
            hh = arrow(bstart+.05, bend+.05, 'LineWidth', 2, 'Color', [.3 .5 1]);
        end
    end

    %% X -> Y direct
    stats_idx = 3;
    [cstart cend] = compute_line_points(Xloc, Xpts, Yloc, Ypts);
    if stats.p(1, stats_idx) < .05
        plot_arrow(cstart, cend);
    else
        line([cstart(1) cend(1)], [cstart(2) cend(2)], 'LineWidth', .5, 'Color', [.7 .7 .7])
    end
    loc_text = Xloc + (Yloc - Xloc) .* .5 + [0 +.15];
    plot_stats_string(loc_text, stats, stats_idx, 'VerticalAlignment', 'top');
%
    %add path ID
    text(loc_text(1)-.05, loc_text(2)+.05, stats.names{stats_idx}, 'FontSize', 16,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    %%
    % plot a blue line if we have level 2 moderators as well
    if size(stats.p, 1) > 1  
        if stats.p(2, stats_idx) < .05
            hh = arrow(cstart+.05, cend+.05, 'LineWidth', 2, 'Color', [.3 .5 1]);
        end
    end

    %% AB
    stats_idx = 5;
    if stats.p(1, stats_idx) < .05
        plot_arrow(aend, bstart);
    end
    loc_text = aend + (bstart - aend) .* .5 + [0 .05];
    plot_stats_string(loc_text, stats, stats_idx, 'VerticalAlignment', 'bottom');
    
    %add path ID
    text(loc_text(1)-.05, loc_text(2)-.2, stats.names{stats_idx}, 'FontSize', 16,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    %%
    % plot a blue line if we have level 2 moderators as well
    if size(stats.p, 1) > 1  
        if stats.p(2, stats_idx) < .05
            hh = arrow(aend-.05, bstart-.05, 'LineWidth', 2, 'Color', [.3 .5 1]);
        end
    end
    
    %%
    %add multilevel results.
    if mlevel
        %% path a2
         %% X -> M2
    stats_idx = 6;
    [astart aend] = compute_line_points(Xloc, Xpts, M2loc, M2pts);
    if stats.p(1, stats_idx) < .05
        plot_arrow(astart, aend);
    end
    if stats.Cov~=1
    loc_text = Xloc + (M2loc - Xloc) .* .5 + [.1 -.15];
    else % move higher
        loc_text = Xloc + (M2loc - Xloc) .* .5 + [.1 .55];
        a2_text= loc_text; % use to mirror next (b2).
    end
%     plot_stats_string(loc_text, stats, stats_idx, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline');
%%
% text(loc_text(1)-.1, loc_text(2)-.15, 'a2', 'FontSize', 16,'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');    
    %% path b2
         %% M2 -> Y
     %% M -> Y
    stats_idx = 7;
    [bstart bend] = compute_line_points(M2loc, M2pts, Yloc, Ypts);
    if stats.p(1, stats_idx) < .05
%         plot_arrow(bstart, bend);
    end
    if stats.Cov~=1
        loc_text = M2loc + (Yloc - M2loc) .* .5 + [-.05 -.15];
    else %mirror the position used previously:
        loc_text = [M2loc(1)+[M2loc(1)-a2_text(1)],  M2loc(2)-[M2loc(2)-a2_text(2)]];
    end
    %%
    plot_stats_string(loc_text, stats, stats_idx, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');
    
    text(loc_text(1), loc_text(2)-.15, 'b2', 'FontSize', 16,'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');    
    %% path ab2
     stats_idx = 8;
    if stats.p(1, stats_idx) < .05
        plot_arrow(aend, bstart);
    end
    loc_text = aend + (bstart - aend) .* .5 + [0 -.15];
%     plot_stats_string(loc_text, stats, stats_idx, 'VerticalAlignment', 'bottom');
    
    %path def
    text(loc_text(1)-.1, loc_text(2)+.2, 'ab2', 'FontSize', 16,'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');    
    end
    cleanup_axis();
end



function setup_figure()
%     fh = create_figure('Path_Diagram');
        figure();
    
% %     fh = findobj('Tag', 'Path_Diagram');
% %     if isempty(fh)
% %         figure('Color', 'w', 'Tag', 'Path_Diagram');
% %         set(gca, 'DefaulttextFontSize', 18, 'DefaulttextHorizontalAlignment', 'center');
% %     else
% %         set(0, 'CurrentFigure', fh);
% %         cla();
% %     end

set(gca, 'DefaulttextFontSize', 18, 'DefaulttextHorizontalAlignment', 'center');
    hold on;
    axis off;
    axis equal;
end

function cleanup_axis()
    axis tight
    limits = axis;
    limits = limits * 1.05;
    axis(limits);
end

function pts = plot_region_circ(name, loc,col)
    %pts = circleToPolygon([loc(1) loc(2) .5], 200);
    pts = circleAsPolygon([loc(1) loc(2) .5], 200);
    fillPolygon(pts, [.8 .8 .8])
    h = drawPolygon(pts, 'k');
    set(h, 'LineWidth', 2, 'Color', [col]);
    text(loc(1), loc(2), name);
end

function [st en] = compute_line_points(loc1, pts1, loc2, pts2)
    l = createLine(loc1, loc2);
    pi1 = intersectLinePolygon(l, pts1);
    pi2 = intersectLinePolygon(l, pts2);
    dist = distancePoints(pi1, pi2);
    [i j] = find(dist == min(dist(:)));
    st = pi1(i,:);
    en = pi2(j,:);
end

function [st en] = plot_arrow(varargin)
    if(length(varargin) == 4)
        [st en] = compute_line_points(varargin{:});
    elseif(length(varargin) == 2)
        [st en] = deal(varargin{:});
    end
    arrow(st, en, 'LineWidth', 2);
end

function plot_stats_string(loc_text, stats, stats_idx, varargin)
    if stats.p(1, stats_idx) < .001, sigstr = '***';
    elseif stats.p(1, stats_idx) < .01, sigstr = '**';
    elseif stats.p(1, stats_idx) < .05, sigstr = '*';
    else sigstr = ''; end
    if stats.mean(stats_idx) < 0.01
        str = sprintf('%5.4f (%5.4f) %s', stats.mean(stats_idx), stats.ste(stats_idx), sigstr);
    else
        str = sprintf('%3.2f (%3.2f) %s', stats.mean(stats_idx), stats.ste(stats_idx), sigstr);
    end
    text(loc_text(1), loc_text(2), str, 'FontSize', 12, varargin{:});
end