function out_h = box_and_scatter(data, params)
%% plots box, whiskers, and scatter, 
% data = points x groups
% params = cols etc

out_h=[];
% test params. 
useCols = params.cols;
    if size(useCols,1) ~= size(data,2);
        disp(['Warning: incorrect colour information'])
        useCols = repmat(params.cols(1,:), size(data,2),1);
    end
    
%
if ~isfield(params, 'scatSize');
    params.scatSize= 100;
end


% sep col for scatter?
if ~isfield(params, 'scatCol');
    params.scatCol= useCols;    
end
scatCols = params.scatCol;
% if un matched sizes, repeat the scatCols.
if size(scatCols,1) ~= size(useCols,1);
    scatCols = repmat(scatCols,size(useCols,1),1);
end


b= boxplot(data, 'color', 'k', 'whisker', 10);
out_h.box= b;

% for each col, adjust appearance:
for ib = 1:size(data,2)
    h=findobj(gca, 'tag','Box');
    % for some reason, the handle populates backwards (from last plotted),
    %so flip the handle, to keep cols consistent.
    h=flipud(h);
    h(ib).LineWidth = 2;
    h(ib).Color = useCols(ib,:);
    h(ib).MarkerFaceColor = useCols(ib,:);
    
    % fill the boxes:
    
    ph=patch(h(ib).XData, h(ib).YData, useCols(ib,:));
    ph.FaceAlpha=.2;
    
    h=findobj(gca, 'tag','Median');
    h(ib).LineWidth = 2;
    h=findobj(gca, 'tag','Upper Whisker');
    h(ib).LineWidth = 2;
    h=findobj(gca, 'tag','Lower Whisker');
    h(ib).LineWidth = 2;
    h=findobj(gca, 'tag','Lower Adjacent Value');
    h(ib).LineWidth = 2;
    h=findobj(gca, 'tag','Upper Adjacent Value');
    h(ib).LineWidth = 2;
    
    out_h.ibox(ib).boxnwhisk = h;

end
 %% add individual data points using scatter
if params.plotScatter        
    hold on;
    
    for id=1:size(data,2)
        useD = data(:,id);
        jitter = zscore(rand(length(useD), 1))/40;
        sc=scatter(id+jitter, useD);
        sc.LineWidth = 2;
        sc.SizeData = params.scatSize;
        sc.MarkerFaceColor = scatCols(id,:);
        sc.MarkerEdgeColor = scatCols(id,:);
        sc.MarkerFaceAlpha = .4;
        sc.MarkerEdgeAlpha=.4;
        %
    end
end
      
end

