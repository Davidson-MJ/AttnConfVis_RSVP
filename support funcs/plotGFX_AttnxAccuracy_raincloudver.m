function plotGFX_AttnxAccuracy_raincloudver(homedir,datadir)
cd(homedir)

%plots all in one figure, outputs in figure directory.
cd(homedir )
cd(datadir)
xlabsare={'Conf', 'PAS'};
xprints={'confidence', 'visibility'};
AUCtype={'ROC', 'Precision-Recall'};
plotROCorRecall=1;

plotnback=0; % change to 1 to plot equivalent data using attention value from n-trials back.
%%
figure(1);    clf;
set(gcf, 'Units', 'normalized','Position', [.5 0 .5 1], 'color', 'w');

counter=1;
for iExp=1:2  % for both experiments
    xlabis = xlabsare{iExp};
    
    cd(datadir)
    
    loadvarname= [ xlabis '_Attn_Allpp.mat'];
    load(loadvarname, 'all_AUC*', 'all_Accuracy*', 'nshift');
    
    
    %%
    
    for ianalysistype=1:2
        
        switch ianalysistype
            case 1  % use Accuracy measures
                useMedianData=all_Accuracybymedian;
                useQdata= all_Accuracybyquantile;
                %                     ylabis=  {['Accuracy'];['(TP+TN)/(P+N)']};
                titletouse = 'Accuracy';
            case 2 % use AUC measures
                useMedianData=all_AUCbymedian;
                useQdata= all_AUCbyquantile;
                %                     ylabis=  {['AUC, ' AUCtype{plotROCorRecall}]};
                titletouse = 'AUROC';
                
        end
        for icat=2%1:2 % median split or quantile
            switch icat
                case 1
                    plotd= useMedianData;
                    
                    
                case 2
                    
                    plotd= useQdata;
                    
            end
            
            subplot(2,2,counter)
            %         title(titletouse); hold on
            
            %% raincloud plots:
            cmap = cbrewer('qual', 'Paired',12);
            colormap(cmap)
            %%
            [~,plotd_adj] =CousineauSEM(plotd);
            dataX=[];
            dataX{1} = plotd_adj(:,3); % note the inverted order for plotting
            dataX{2} = plotd_adj(:,2);
            dataX{3} = plotd_adj(:,1);
            
            hold on
            %%
%             clf
            h=rm_raincloud(dataX', cmap(8+ianalysistype,:), 1, 'ks');
            
            set(gcf, 'color', 'w');
            %% Label axes.
            set(gca, 'yticklabel', {'Lowest', 'Medium', 'Highest'}); shg
            ylabel('Attention')
            xlabel(titletouse)
            axis tight
            set(gca, 'fontsize', 25)
            %%
            counter=counter+1;
            
            if ianalysistype==1 % set xlim based on accuracy level.
                axis tight;
                xlimsare = get(gca, 'xlim');
            end
        end
        
        %set xaxis to match first plot
        xlim([xlimsare(1) .95]);
    end
    cd(homedir); cd ../
    
    cd(['Documents/Figures/Exp' num2str(iExp+1) ' Attention and Accuracy per participant'])
    
    printname=['Accuracy by attention, all participants, ' xlabis ' exp, raincloud ver'];
    
    
    
    
    
    %         print('-dpng',printname );
    
end