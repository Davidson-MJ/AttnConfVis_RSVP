% function plotPP_ScatterandXY(homedir,datadir)
cd(homedir)

% plots all in one figure, outputs in figure directory.

%% set colours and marker types (match main fig.).
cmapB= cbrewer('qual', 'Paired', 12);

%plot parameters (H,M,FA,CR):

col1 = cmapB(8,:);
col2 = cmapB(10,:);

markertype ='o';
coltype = {col1,col1,col2,col2};
linewtype=2;
fontsize= 10;
for iExp=2%:2
    
    % Plot figure with subplots per ppant.
    figure(1); clf
    set(gcf, 'Units', 'normalized','Position', [0 0 1 1], 'color', 'w');
    
    if iExp==1
        Axis_X = -100:100; % confidence measure.
        Axis_Y = Axis_X; % Attention measure.
        xlabis='Conf';
        xlabeluse = 'Confidence';
        xlimslabels = {'Sure absent', 'Sure present'};
        loadnppants=[1:12];
        
        %where to plot the axes (which ppants).        
        xAX_on=[9:12];
        yAX_on=[1,5,9];
        
        
            %create tight subplots to deal with space between axes:
            ha = tight_subplot(3,4,[ .03, .03], [.1 .01], [.1 .1]);
    else
        Axis_X = 0:200; % confidence measure.
        %             Axis_X = 0:200; % PAS measure.
        Axis_Y= -100:100;% as above.
        xlabis='PAS';
        xlabeluse= 'Visibility';
        xlimslabels = {'None', 'All'};
        loadnppants=[2:10];
        
         %create tight subplots to deal with space between axes:
            ha = tight_subplot(3,3,[ .03, .03], [.1 .01], [.1 .1]);
         %where to plot the axes (which ppants).        
        xAX_on=[7:9];
        yAX_on=[1,4,7];
    end
    cd(datadir)
    %how many data types/ppants this folder.
    nppants = length(loadnppants);
    
    %
    
    
    for ippant=1:length(loadnppants)
        cd(datadir)
        % open p_table
        loadvarname= [ xlabis '_Attn_participant' num2str(loadnppants(ippant)) '.mat'];
        load(loadvarname);
        
       
        %%%%%%%%
        
        sc_lg=[];
        sc_counts=zeros(1,4);
        %for each metric, plot scatter, and SDT vectors:
        axes(ha(ippant));
        
        for iSDT= 1:4% plot each response type (H,M,FA,CR)
            
            %% scatterplot first.
            
            %create subplot
            
            hold on
            % create tight subplots to deal with the space between axes.
            
            
            plotresps=find(logical(SDTindex(iSDT,:)));
            sc=scatter(XY_RESP_data(plotresps,1),XY_RESP_data(plotresps,2), '*');
            
            alpha(.15);
            
            sc.Marker = markertype;
            sc.MarkerEdgeColor= coltype{iSDT};
            sc.MarkerFaceColor= coltype{iSDT};
            sc_lg(iSDT)=sc;
            sc_counts(iSDT)=  length(plotresps);
            
            %labelling:
            tt=    title(['Participant ' num2str(ippant)], 'VerticalAlignment', 'bottom');
            
               
                if iExp==1
                set(gca, 'XTick', [-100:25:100]) % add increments to space labels.
                else
                    set(gca, 'XTick', [0:25:200]) % add increments to space labels.
                end
                set(gca, 'YTick', [-100:25:100]) % add increments to space labels.
                
                gt= get(gca);
                %Y
                Yticklabeln = cell(1, length(gt.YTick));
                Yticklabeln{1} ='Less focused';               
                Yticklabeln{end} ='More focused';               
                % X
                Xticklabeln = cell(1, length(gt.XTick));
                Xticklabeln{2} =xlimslabels{1};
                Xticklabeln{end-1} =xlimslabels{2};     
                
                % declutter the axes across PP
                if ismember(ippant, xAX_on);         
                   xlabel(xlabeluse, 'fontweight', 'bold', 'VerticalAlignment', 'middle','fontsize', fontsize*1.5)
                set(gca, 'XtickLabel', Xticklabeln, 'XTickLabelRotation', 0);
                h=gca; h.XAxis.TickLength = [0 0];
                else
                    set(gca, 'XtickLabel', [], 'Xtick',[]);
                end
                
                if ismember(ippant, yAX_on);
                  ylabel('Attention', 'fontweight', 'bold',  'rotation', 0, 'HorizontalAlignment', 'left','fontsize', fontsize*1.5);
                  h=gca; h.YAxis.TickLength = [0 0];
                   
                set(gca, 'YtickLabel', Yticklabeln);
                else 
                    set(gca, 'YtickLabel', [],'Ytick', []);
                end
                box on
                axis square
                xlim([Axis_X(1) Axis_X(end)])
                 
                plot(xlim,[0 0], 'k-')
                plot([0 0], ylim, 'k-')
                %%
                if iSDT==4 && ippant ==4 % % legend once.
%                     lg=legend([sc_lg(2),sc_lg(3)],{['Target'], ['No Target']'});
                    %
%                     set(lg, 'location', 'SouthEastOutside')
                    
                end
 
          
             
          
        end % after each
       
       
    end % ppant
     %%
        cd(homedir); cd ../
        cd(['Documents/Figures/Exp' num2str(iExp+1) ' BEH resp per participant'])
        
        print('-dpng',['Click map All Participants Supp Fig_' xlabis ]);
        
        
end % exp