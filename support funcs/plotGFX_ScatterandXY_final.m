% function plotGFX_ScatterandXY_final(homedir,datadir)
cd(homedir)
dbstop if error
%plots all in one figure, outputs in figure directory.
%we are plotting present vs absent as our main distinction.
%%
cmapB= cbrewer('qual', 'Paired', 12);
%plot parameters (H,M,FA,CR):

col1 = cmapB(8,:);
col2 = cmapB(10,:);

markertype ='o';
coltype = {col1,col1,col2,col2};
linewtype=2;%
%
errorsonly=0; % if wanting to plot data for misses and FA only.
plotXYdistributions=0; % =1 if want subplots to captre distribution of each response class.

plotXYcorrelation =1;
plotRhosummary = 1;

   fontsize=25;
for iExp=1%:2
    if iExp==1
        Axis_X = -100:100; % confidence measure.
        Axis_Y = Axis_X; % Attention measure.
        xlabis='Conf';
        xlabeluse = 'Confidence';
        xlimslabels = {'Sure absent', 'Sure present'};
        
    else
        Axis_X = 0:200; % confidence measure.
        %             Axis_X = 0:200; % PAS measure.
        Axis_Y= -100:100;% as above.
        xlabis='PAS';
        xlabeluse = 'Visibility';
        xlimslabels = {'None', 'All'};
    end
    cd(datadir)
    
    
    loadvarname= [ xlabis '_Attn_Allpp.mat'];   
    load(loadvarname);
    
    nppants = size(allSDTindex,1);
   
 
    %%%%%%%%
    %% plot for all participants, overlaying scatter.
    
    figure(iExp); clf
    set(gcf, 'units', 'normalized', 'position', [ 0 0 1 1])
    sc_lg=[];
    sc_counts=zeros(1,4);
    
    % use these for plot legends on X and Y axes:
    [lgAx_Y,lgAx_X]=deal([]);
    %
    if errorsonly==0
        iSDTplots=1:4;
    else        
        iSDTplots = 2:3; % errors only
    end
    %%
   
    for iSDT= iSDTplots
        
        for ippant=1:nppants
            
            SDTindex = squeeze(allSDTindex(ippant,:,:));
            XY_RESP_data= squeeze(alltoScatter(ippant,:,:));
            %% scatterplot first.            
            plotresps = find(logical(SDTindex(iSDT,:)));
            hold on
%             subplot(2,3,2:3)
            hold on
            sc=scatter(XY_RESP_data(plotresps ,1),XY_RESP_data(plotresps ,2), '.');
            sc.Marker = markertype;
            sc.MarkerEdgeColor= coltype{iSDT};
            sc.MarkerFaceColor= coltype{iSDT};
            sc.MarkerFaceAlpha= .5;
            sc_lg(iSDT)=sc;
            sc_counts(iSDT)= sc_counts(iSDT)+length(plotresps);

            
            if iSDT==4 && ippant == size(alltoScatter,1) % fill legend and axis features at end.
           
             
                %%

                %trick to increase size of legend markers:
                %create copy of objs.
            
                lg=legend([sc_lg(2),sc_lg(3)],{['target present'], ['target absent']'});                
                set(lg, 'location', 'SouthEast', 'fontsize', fontsize)

                ylim([-101 101])                
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                
                axis square
                
                gt= get(gca);
                %Y
                Yticklabeln = cell(1, length(gt.YTick));
                Yticklabeln{1} ='Less focused';
                Yticklabeln{end} ='More focused';               
                % X
                  %%
%                 set(gca, 'YtickLabel', Yticklabeln, 'fontsize', fontsize);
%                 xlabel(xlabeluse, 'fontweight', 'bold','fontsize', fontsize*1.5)                
%                 ylabel('Attention', 'fontweight', 'bold',  'rotation', 0,'fontsize', fontsize*1.5) 
              
                
                %% matlab sometimes squeezes the last x tick off, so reset:
                gt= get(gca);
                               
                Xticklabeln = cell(1, length(gt.XTick));
                Xticklabeln{1} =xlimslabels{1};
                Xticklabeln{end} =xlimslabels{2};    % - 1to account for axis square forcing symmetry below.
%               
%                 set(gca, 'YtickLabel', Yticklabeln, 'XtickLabel', Xticklabeln);
                set(gca, 'YtickLabel', [], 'XtickLabel', []);
            end
            
          
        end
%         
   
    end
    %% plot the resulting click maps.
    
    set(gcf, 'color', 'w')
    %%
    cd(homedir)
    cd ../
    %%
%     cd(['Documents/Figures/Exp' num2str(iExp+1) ' BEH resp per participant'])
    printnametmp = ['Click map All Participants_' xlabis '_final'];
    if errorsonly==1
    printname= [printnametmp 'errorsonly'];
    else
        printname=printnametmp;
    end    
    print('-dpng',printname);
    
    %%
        
    if plotXYdistributions==1
        
          %% now show density across axes in new plots.
        %define new coordinate positions.                
        pos = get(gca,'position');      
        figure(iExp+2); clf;
        set(gcf, 'units', 'normalized', 'position', [ 0 0 1 1])
        %transpose to new figure
        % % we want to plot all present vs absent, on X and Y axes.
        %and show the variance across ppants.
        
        % might need to bin the result:
        binsare = ceil(linspace(1,size(allSDT_on_Xdimension,3),10));
         
            [binsforPlot,errsforPlot]= deal(nan(4,10));            
            count=1;
        
            for XYax=[1,2]
            
            switch XYax
                case 1                    
                    plotme=allSDT_on_Ydimension;
                    plotx='Attention';
                    axis_is=Axis_Y;
                    
                case 2
                    plotme=allSDT_on_Xdimension;
                    plotx= xlabis;
                    axis_is=Axis_X;
            end
           
              for iTargClass = 1:2 % present vs absent.                
                %first collapse within ppant types.
                switch iTargClass
                    case 1
                Pdata = squeeze(nansum(plotme(:,[1,2], :),2)); % present
                colis= col1;
                    case 2
                Pdata = squeeze(nansum(plotme(:,[3,4], :),2)); % absent
                colis= col2;
                end
                
                figure(iExp+2)
                subplot(2,2,XYax); hold on;
%%%%% plot shaded error bar (mean across pps)            
%             % mean to plot
            plotM = squeeze(nanmean(Pdata,1));            
            % errorbar to plot
            wSEM = CousineauSEM(Pdata);
            
            %now plot
            sh=shadedErrorBar(axis_is, plotM, wSEM, [], 1);
            
            %adjust appearance
            sh.mainLine.Color= colis;
            sh.patch.FaceColor= colis;
            sh.mainLine.LineWidth= linewtype(1);

            ylabel('Count')
            xlabel(plotx);
            set(gca,'fontsize', fontsize)
              
            axis tight
            
            %% now prepare to plot binned version:
            binsare = ceil(linspace(1,size(allSDT_on_Xdimension,3),10));
            plotBins = zeros(size(Pdata,1),length(binsare));
            for ibin = 1:(length(binsare)-1)
%                 plotBins(:,ibin) = squeeze(nansum(Pdata(:, [binsare(ibin):binsare(ibin+1)-1]),2));            
                plotBins(:,ibin) = squeeze(nanmean(Pdata(:, [binsare(ibin):binsare(ibin+1)-1]),2));            
            end            
            %add final column
%             plotBins(:,length(binsare)) = squeeze(nansum(Pdata(:, [binsare(length(binsare)-1):binsare(length(binsare))]),2));            
            plotBins(:,length(binsare)) = squeeze(nanmean(Pdata(:, [binsare(length(binsare)-1):binsare(length(binsare))]),2));            
            
%%
            mB = squeeze(nanmean(plotBins,1));            
            stEB = CousineauSEM(plotBins);
            %% store to prep plot:
            binsforPlot(count,:) = mB;
            errsforPlot(count,:) = stEB;
            count= count+1;
%              subplot(2,2,XYax+2);
%             bh= bar(mB);hold on
%             bh.FaceColor = colis;
%             bh.FaceAlpha= .1;
%             eh=errorbar(1:10, mB(1:10), stEB(1:10));
%             eh.LineStyle = 'none';            
              end
              % now plot bars (stacked)
              %yaxis
                subplot(2,2,3);
                tmp = [binsforPlot(1,:);binsforPlot(2,:)];
                tmpE = [errsforPlot(1,:);errsforPlot(2,:)];
            bh= bar(tmp');hold on
            bh(1).FaceColor = col1;
            bh(2).FaceColor = col2;
            errorbar_groupedfit(tmp', tmpE');

 ylabel('Count')
            xlabel('Attention bins');
            set(gca,'fontsize', fontsize)
              
            axis tight
% 
%%
                subplot(2,2,4); % xaxis
                tmp = [binsforPlot(3,:);binsforPlot(4,:)];
                tmpE = [errsforPlot(3,:);errsforPlot(4,:)];
            bh= bar(tmp');hold on
            bh(1).FaceColor = col1;
            bh(2).FaceColor = col2;
            errorbar_groupedfit(tmp', tmpE');
             ylabel('Count')
            xlabel([ xlabeluse ' bins']);
            set(gca,'fontsize', fontsize)
              
            axis tight
            %%     
          legend('Present', 'absent')
        end
        
        
    
    end % plot XY distrib
    
    
    %% plot heat map, folding on Y axis.
    if plotXYcorrelation ==1 % (heatmap)
      %% now show density across axes in new plots.
        %define new coordinate positions.                
        pos = get(gca,'position');      
        figure(iExp+4); clf;
        set(gcf, 'units', 'normalized', 'position', [ 0 0 1 1])
        %transpose to new figure
        % % we want to plot all present vs absent, on X and Y axes.
        %and show the variance across ppants.
        
        
            %%
            clf
            
            cmapC= cbrewer('seq', 'YlOrRd', 15);
            
            interpval = 10;
            % all targets for now:
            allX = squeeze(alltoScatter(:,:,1));
            allY = squeeze(alltoScatter(:,:,2));
            
            if iExp==1 % take abso confidence.
                allX = abs(allX);
                xbins=0:interpval:100;
                xlimslabels = {'Unsure' , 'Sure'};
                xlabeluse = '| Confidence |';
            else
                xbins=0:interpval:200;
            end
                
            ybins = Axis_Y(1):interpval:Axis_Y(end);
            %
            % map vals to inds
            % over-range vals are truncated instead of extrapolated
            xi = interp1(xbins, 1:numel(xbins), allX, 'nearest', numel(xbins));
            yi = interp1(ybins, 1:numel(ybins), allY, 'nearest', numel(ybins));
            
            % make 2-d hist (heatmap) using accumarray
            H = accumarray([yi(:), xi(:)], 1, [numel(ybins), numel(xbins)]);
            
            %
            % generate plot
            
            ha = axes;
          
            axis([xbins(1), xbins(end), ybins(1), ybins(end)]);
            grid(ha, 'on');
            % plot heat map
            hi = imagesc(xbins, ybins, H);
            set(ha, 'YDir', 'normal'); % correct y-axis direction after imagesc
            
            hold on;
            % plot circles
            [Xbins, Ybins] = meshgrid(xbins, ybins);
            ind = H(:) > 0;
            sc=scatter(Xbins(ind), Ybins(ind), H(ind).^1.5, 'k');
            sc.LineWidth = 2;
            colormap(cmapC);
            c=colorbar;
            ylabel(c, 'Click count')
            ylim([-101 101])
            plot(xlim,[0 0], 'k-')
            plot([0 0], ylim, 'k-')
            
            gt= get(gca);
            %Y
            Yticklabeln = cell(1, length(gt.YTick));
            Yticklabeln{1} ='Less focused';
            Yticklabeln{end} ='More focused';
            % X
            Xticklabeln = cell(1, length(gt.XTick));
            Xticklabeln{1} =xlimslabels{1};
            Xticklabeln{end} =xlimslabels{2};
            %                 %%
%             set(gca, 'YtickLabel', [],'fontsize', fontsize, 'XtickLabel', Xticklabeln);
            
%             xlabel(xlabeluse, 'fontweight', 'bold','fontsize', fontsize*1.5)
%             ylabel('Attention', 'fontweight', 'bold',  'rotation', 0,'fontsize', fontsize*1.5)
            axis square
            set(gcf, 'color', 'w');
            %% matlab sometimes squeezes the last x tick off, so reset:
            gt= get(gca);
            
            Xticklabeln = cell(1, length(gt.XTick));
            Xticklabeln{1} =xlimslabels{1};
            Xticklabeln{end} =xlimslabels{2};    % - 1to account for axis square forcing symmetry below.
            %
%             set(gca, 'YtickLabel', [], 'XtickLabel', Xticklabeln);
            set(gca, 'YtickLabel', [], 'XtickLabel', [], 'fontsize', fontsize*1.5);
            
            
            
            
           %%
           
           if plotRhosummary 
           plot_RHOsummary;
           end
    end % plot correlation
end % iExp
    
% end %function
%%


