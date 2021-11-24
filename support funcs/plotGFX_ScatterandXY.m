function plotGFX_ScatterandXY(homedir,datadir)
cd(homedir)

%plots all in one figure, outputs in figure directory.
%%

errorsonly=0; % if wanting to plot data for misses and FA only.
plotXYdistributions=1; % =1 if want subplots to captre distribution of each response class.
plotRTsin3D =0; % plot surf map of RTs also
for iExp=1:2
    if iExp==1
        Axis_X = -100:100; % confidence measure.
        Axis_Y = Axis_X; % Attention measure.
        xlabis='Conf';
        xlabeluse = 'Confidence';
        
        
    else
        Axis_X = 0:200; % confidence measure.
        %             Axis_X = 0:200; % PAS measure.
        Axis_Y= -100:100;% as above.
        xlabis='PAS';
        xlabeluse = 'Visibility';
    end
    cd(datadir)
    
    
    loadvarname= [ xlabis '_Attn_Allpp.mat'];   
    load(loadvarname);
    
    %plot parameters:
    markertype ={'o','x','x','o'};
    coltype = {'b','b','r','r'};
    linewtype=[2,1,1,2]; % to differentiate subplots.
    fontsize=15;
    %%%%%%%%
    %% plot for all participants, overlaying scatter.
    clf
    set(gcf, 'units', 'normalized', 'position', [ 0 .7 .4 .6])
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
    
    for iSDT= iSDTplots
        
        for ippant=1:size(alltoScatter,1)
            
            SDTindex = squeeze(allSDTindex(ippant,:,:));
            XY_RESP_data= squeeze(alltoScatter(ippant,:,:));
            %% scatterplot first.
            
            plotresps = find(logical(SDTindex(iSDT,:)));
            hold on
            subplot(2,3,2:3)
            hold on
            sc=scatter(XY_RESP_data(plotresps ,1),XY_RESP_data(plotresps ,2), '.');
            sc.Marker = markertype{iSDT};
            sc.MarkerEdgeColor= coltype{iSDT};
            sc.MarkerFaceColor= coltype{iSDT};
            sc_lg(iSDT)=sc;
            sc_counts(iSDT)= sc_counts(iSDT)+length(plotresps);
%             title('All participants')
            
            
            if iSDT==4 % fill legend and axis features at end.
                
                alpha(.15);
                xlabel(xlabeluse, 'fontweight', 'bold')
%                 lg=legend([sc_lg(1),sc_lg(2),sc_lg(3),sc_lg(4)],...
%                     {['Hits (' num2str(sc_counts(1)) ')'], ['Misses (' num2str(sc_counts(2)) ')']', ...
%                     ['FA (' num2str(sc_counts(3)) ')'], ['CR (' num2str(sc_counts(4)) ')']}, 'Autoupdate', 'off');
                
%legend without the counts:
                lg=legend([sc_lg(1),sc_lg(2),sc_lg(3),sc_lg(4)],...
                    {['Hits'], ['Misses']', ...
                    ['FA'], ['CR']}, 'Autoupdate', 'off');
                
                set(lg, 'location', 'SouthEast')
                yt=ylabel('Attention','fontweight', 'bold', 'rotation', 0);
                ylim([-101 101])
                %                 xlim([min(Axis_X) max(Axis_X)]);
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                set(gca,'fontsize', fontsize)
            end
            %  debugging (plotting errors only)
            if errorsonly==1 && iSDT==3
                xlabel(xlabeluse, 'fontweight', 'bold')
                  lg=legend([sc_lg(2),sc_lg(3)],...
                    {['Misses (' num2str(sc_counts(2)) ')']', ...
                    ['FA (' num2str(sc_counts(3)) ')']}, 'Autoupdate', 'off');
                
                set(lg, 'location', 'SouthEast')
                yt=ylabel('Attention','fontweight', 'bold', 'rotation', 0);
                ylim([-101 101])
%                 xlim([min(Axis_X) max(Axis_X)]);
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                set(gca,'fontsize', fontsize)
            end
            xlim([Axis_X(1) Axis_X(end)])
            %%
            if iExp==1 % label x axis limits
  
                set(gca,'Ytick',[-100 100],'Xtick', [Axis_X(1) Axis_X(end)], 'xticklabel', {'Sure absent', 'Sure present'})
            else
                set(gca,'Ytick',[-100 100],'Xtick', [Axis_X(1) Axis_X(end)], 'xticklabel', {'None', 'All'})                
            end
            set(gca, 'Yticklabel', {'Less Focused', 'More Focused'});
        end
        
      
        %
        %%
        
    if plotXYdistributions==1
          %% now show density across axes.
        %define new coordinate positions.
        pos = get(gca,'position');
        q = [pos(1) pos(2) 0 0];
        s = [pos(3) pos(4) pos(3) pos(4)];
        
        for XYax=1:2
            switch XYax
                case 1
                    
                    plotme=squeeze(allSDT_on_Xdimension(:,iSDT,:));
                    plotx= xlabis;
                    axis_is=Axis_X;
                    
                    %X axis subplot:
                    %define new subplot position.
                    
                    if mod(iSDT,2) ==1 % plot the responded as 'present' types together
                        
                        % [TL xdev, TL ydev, %width, %height]
                        newpos=[0 -.55 1 .45].*s  +q;
                        
                    else
                        %plot target absent types together
                        %increase starting dist.
                        newpos=[0 -1.15 1 .45].*s  +q;
                        
                    end
                    
                    %if errors, just throw on some subplot:
                    if errorsonly==1; newpos=[0 -.55 1 .45].*s  +q; end

                    
                case 2
                    plotme=squeeze(allSDT_on_Ydimension(:,iSDT,:));
                    plotx='Attention';
                    axis_is=Axis_Y;
                    
                    %Y axis subplot:
                    if  mod(iSDT,2) ==1 % plot the responded as 'present' types together
                        
                        % [TL xdev, TL ydev, %width, %height]
                        newpos=[-.3 0 .25 1].*s  +q;
                        
                    else
                        %plot target absent types together
                        %increase starting dist.
                        newpos=[-.65 0 .25 1].*s  +q;
                        
                    end
                    %if errors only ,same plot:
                    if errorsonly==1; newpos=[-.3 0 .25 1].*s  +q; end
                    
                    
            end
            %%
            %draw subplot
            h(XYax) = subplot('Position',newpos);
            hold on;
            
%%%%% plot nansum:
plotM = squeeze(nansum(plotme,1));
sh.mainLine= plot(axis_is, plotM, 'color', coltype{iSDT}, 'linew',linewtype(iSDT));


%%%%% plot shaded error bar (mean across pps)            
%             % mean to plot
%             plotM = squeeze(nansum(plotme,1));
            
            % errorbar to plot
%             wSEM = CousineauSEM(plotme);
            
            %now plot
%             sh=shadedErrorBar(axis_is, plotM, wSEM, [], 1);
            
            %adjust appearance
%             sh.mainLine.Color= coltype{iSDT};
%             sh.patch.FaceColor= coltype{iSDT};
%             sh.mainLine.LineWidth= linewtype(iSDT);
            
            
            %flip if necessary.
            if XYax==2
                view([90 -90]);
                %store properties of mainLines for legend.
                
%                 lgAx_Y=[lgAx_Y, sh.mainLine];
            else
                lgAx_X=[lgAx_X, sh.mainLine];    
            end
            ylabel('Count')
            xlabel(plotx);
            set(gca,'fontsize', fontsize)
            
            
              ylabel('Count')
            xlabel(plotx);
            set(gca,'fontsize', fontsize)
              
            axis tight
            
        end
        %% add legend to bottom axes
        
        if iSDT==4 && errorsonly==0
           legend(h(1), [lgAx_X(2) lgAx_X(4)], {dimsare{2}, dimsare{4}})
        elseif iSDT==3 && errorsonly==0
            
            legend(h(1), [lgAx_X(1) lgAx_X(3)], {dimsare{1}, dimsare{3}})
        end
        %%
    else % not XY distributions , instead simple scatter
        
        
    end
    %%
%     legend off
    set(gcf, 'color', 'w')
    %%
    cd(homedir)
    cd ../
    %%
    cd(['Documents/Figures/Exp' num2str(iExp+1) ' BEH resp per participant'])
    printnametmp = ['Click map All Participants_' xlabis];
    if errorsonly==1
    printname= [printnametmp 'errorsonly'];
    else
        printname=printnametmp;
    end    
    print('-dpng',printname);
    
    %% %finally, also draw complementary scatter plot showing RTs:
    if plotRTsin3D 
    clf;
    
    %define matrix for RTs based on X and Y vals.
    %combine into one huge vector:
    nppants= size(alltoScatter,1);
    ntrials = size(alltoScatter,2);
    
    allXY=reshape(alltoScatter, [nppants*ntrials, 2]);
    allRT=reshape(allRTtoScatter, [nppants*ntrials,1]);
%% unless we want to look just at errors:
if errorsonly==1
    allXY=[]; allRT=[];
    for iSDT= 2:3
        for ippant=1:size(alltoScatter,1)
            SDTindex = squeeze(allSDTindex(ippant,:,:));
            XY_RESP_data= squeeze(alltoScatter(ippant,:,:));
            %% scatterplot first.
            plotresps = find(logical(SDTindex(iSDT,:)));
            %concatenate per ppant (only subset of responses).
            allXY=[allXY; XY_RESP_data(plotresps,:)];
            allRT = [allRT; allRTtoScatter(ippant,plotresps)'];
        end
    end
end
%%


                dotsize=60;
                c=colormap('jet');
                c=flipud(c);
                
                %just display RTs under 10s:
                allRT_sh= allRT(allRT<10000);
                allXY_sh= allXY(1:length(allRT_sh),:);
                
                %resample original length to new length.               
                origl= length(c); newl = length(allRT_sh);
                c_r=resample(c,newl,origl);
                %sort for colour scaling based on zvalue (RT)
                [nrt,indrt]=sort(allRT_sh);
                
                scatter3(allXY_sh(indrt,1), allXY_sh(indrt,2), nrt ,dotsize, c_r,'filled')
%% set plot info
                cb=colorbar;
                ylabel(cb, 'RT speed');
                crange = caxis;
                set(cb, 'ytick', [crange(1), crange(2)], 'yticklabel', {'slowest', 'fastest'},...
                    'fontsize', 15)
                set(gca, 'fontsize', 15);
                xlabel(xlabis); ylabel('Attention'); zlabel('RT(ms)')

            %%
            view([0 90]) ; % reorient to match view:
            %save adjustable file:
            savenametmp=['Click map All Participants_' xlabis '_wRT'];
            if errorsonly==1
                savename=[savenametmp '_errorsonly'];
                title('errors RT')
            else
                savename=savenametmp;
                title('clicks with RT')
            end
            savefig(savename)
    end
end % iExp
    
end %function
%%


