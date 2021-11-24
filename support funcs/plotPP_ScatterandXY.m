function plotPP_ScatterandXY(homedir,datadir)
cd(homedir)    

% plots all in one figure, outputs in figure directory.
%%
errorsonly=0; % 1 to plot only misses and FA

for iExp=1:2
    
        if iExp==1
            Axis_X = -100:100; % confidence measure.
            Axis_Y = Axis_X; % Attention measure.
            xlabis='Conf';
            loadnppants=[1:12];
        else
             Axis_X = 0:200; % confidence measure.
%             Axis_X = 0:200; % PAS measure.
            Axis_Y= -100:100;% as above.
            xlabis='PAS';
            loadnppants=[2:10];
            
        end
        cd(datadir)    
    %how many data types/ppants this folder.
    nppants = length(loadnppants);        
    
%     dimsare = ['nppants', 'X_Attn', 'SDTm(H,M,FA,CR)', 'axis'];
    for ippant=loadnppants
        cd(datadir)
          % open p_table
          loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
          load(loadvarname);
             
        %plot parameters:
        markertype ={'o','x','x','o'};
        coltype = {'b','b','r','r'};
        linewtype=[2,1,1,2];
        fontsize=15;
        %%%%%%%%
        %plot for this participant
         figure(ippant); clf; 
 %set half screen figure size, and white colour:
 set(gcf, 'Units', 'normalized','Position', [0 0 .56 1], 'color', 'w');  
        sc_lg=[];
        sc_counts=zeros(1,4);
        %for each metric, plot scatter, and SDT vectors:
       if errorsonly==0           
        iSDTplots =1:4;
       else
           iSDTplots=2:3;
       end
           
        for iSDT= iSDTplots
            
            %% scatterplot first.
            
            hold on
            subplot(2,3,2:3)
            hold on
            plotresps=find(logical(SDTindex(iSDT,:)));
            sc=scatter(XY_RESP_data(plotresps,1),XY_RESP_data(plotresps,2), '*');          
            sc.Marker = markertype{iSDT};
            sc.MarkerEdgeColor= coltype{iSDT};
            sc.MarkerFaceColor= coltype{iSDT};
            sc_lg(iSDT)=sc;
            sc_counts(iSDT)=  length(plotresps); 
            title(['Participant ' num2str(ippant)])
            set(gca,'Ytick',[],'Xtick', [])
            
            if iSDT==4 % pop at end.
                
                alpha(.15);
                xlabel(xlabis)
                 lg=legend([sc_lg(1),sc_lg(2),sc_lg(3),sc_lg(4)],...
                    {['Hits (' num2str(sc_counts(1)) ')'], ['Misses (' num2str(sc_counts(2)) ')']', ...
                    ['FA (' num2str(sc_counts(3)) ')'], ['CR (' num2str(sc_counts(4)) ')']}, 'Autoupdate', 'off');
                
                set(lg, 'location', 'SouthEast')
                ylabel('Attention')
                ylim([-101 101])
%                 xlim([min(Axis_X) max(Axis_X)]);
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                set(gca,'fontsize', fontsize)
            end
%  debugging
            if iSDT==3 && errorsonly==1
                xlabel(xlabis)
                   lg=legend([sc_lg(2),sc_lg(3)],...
                    {['Misses (' num2str(sc_counts(2)) ')']', ...
                    ['FA (' num2str(sc_counts(3)) ')']}, 'Autoupdate', 'off');
                
                set(lg, 'location', 'SouthEast')
                ylabel('Attention')
                ylim([-101 101])
%                 xlim([min(Axis_X) max(Axis_X)]);
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                set(gca,'fontsize', fontsize)
            end
            xlim([Axis_X(1) Axis_X(end)])
        %% now show density across axes.
        %define new coordinate positions.
        pos = get(gca,'position');
        q = [pos(1) pos(2) 0 0];
        s = [pos(3) pos(4) pos(3) pos(4)];
       
%        
        %%
        
        
        for XYax=1:2
            switch XYax
                case 1
                    plotme=(SDTmetrics_on_Xaxis(iSDT,:));                    
                    plotx= xlabis;
                    axis_is=Axis_X;
                    
                    %X axis subplot:
                    %define new subplot position.
                  
                    if mod(iSDT,2) ==1 % plot the reponded 'present' types together
            
                        % [TL xdev, TL ydev, %width, %height]          
                        newpos=[0 -.55 1 .45].*s  +q;    
                        
                    else %plot the responded as 'absent' types together.
                        
             %increase starting dist.
                        newpos=[0 -1.15 1 .45].*s  +q; 
                        
                    end
                
                case 2
                    plotme=(SDTmetrics_on_Yaxis(iSDT,:))';
                    plotx='Attention';
                    axis_is=Axis_Y;
                    
                    %Y axis subplot:                                       
                     if mod(iSDT,2) ==1 % plot the responded as 'present' types together
            
                        % [TL xdev, TL ydev, %width, %height]          
                         newpos=[-.3 0 .25 1].*s  +q;  
                        
                    else
                        %plot target absent types together
             %increase starting dist.
                        newpos=[-.65 0 .25 1].*s  +q; 
                        
                    end
                    
            end
            %%       
            %draw subplot
            h(XYax) = subplot('Position',newpos); 
            hold on;
            %draw data            
             plot(axis_is, smooth(plotme),'color', coltype{iSDT}, 'linew', linewtype(iSDT))
            
             %flip if necessary.
             if XYax==2
                view([90 -90]);               
             end
            ylabel('Count')
            xlabel(plotx);
            set(gca,'fontsize', fontsize)
            
         if iSDT==3 && errorsonly==0
%           legend(h(1), dimsare{1}, dimsare{2})        
          legend(dimsare{1}, dimsare{3}, 'Location', 'NorthEast')        
        elseif iSDT==4
            legend(dimsare{2}, dimsare{4}, 'Location', 'NorthEast')        
        end
         axis tight
        end
        %% add legend to bottom axes
        
            %%
        end % after each
        %% add total correlation.
         subplot(2,3,2:3)
         if iExp==1; xvaltext=-80; else, xvaltext=20; end
         
        text(xvaltext, -110, ['rho=' sprintf('%.2f',rho(1,2)) ',p=' sprintf('%.3f',corr_pval(1,2)) ], 'fontsize', fontsize);
        %%
        set(gcf, 'color', 'w')
        %%
        cd(homedir); cd ../
        cd(['Documents/Figures/Exp' num2str(iExp+1) ' BEH resp per participant'])
        if errorsonly==0
        print('-dpng',['Click map Participant ' num2str(ippant) '_' xlabis '_smoothed']);
        else
            print('-dpng',['Click map Participant ' num2str(ippant) '_' xlabis '_smoothed_errorsonly']);
        end

        end
         
    end