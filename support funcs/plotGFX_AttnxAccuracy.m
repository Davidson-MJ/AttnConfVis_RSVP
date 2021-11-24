function plotGFX_AttnxAccuracy(homedir,datadir)
cd(homedir)

%plots all in one figure, outputs in figure directory.
cd(homedir )
cd(datadir)
xlabsare={'Conf', 'PAS'};
xprints={'Confidence', 'visibility'};
AUCtype={'ROC', 'Precision-Recall'};
plotROCorRecall=1;

cmap = cbrewer('qual', 'Paired',12);
colormap(cmap)
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
    
    %      set(gcf, 'Units', 'normalized','Position', [.2 .4 .8 .4], 'color', 'w');
    
    %%
    
    for ianalysistype=1:2
        switch ianalysistype
            case 1  % use Accuracy measures
                useMedianData=all_Accuracybymedian;
                useQdata= all_Accuracybymultiforpred;
                ylabis=  {['Accuracy']}%;['(TP+TN)/(P+N)']};
                titletouse = 'Accuracy by attention';
                col='b';
            case 2 % use AUC measures
                useMedianData=all_AUCbymedian;
                useQdata= all_AUCbymultiforpred;
                %                     ylabis=  {['AUC, ' AUCtype{plotROCorRecall}]};
                ylabis=  {['AUROC']};
                titletouse = 'AUROC by attention';
                col='m';
                
        end
        %match colour to raincloud plots.
        
%         col = cmap(8+ianalysistype,:);
        
        % 
        
        for icat=2% 1:2 % median split or quantile
            switch icat
                case 1
                    plotd= useMedianData;
                    
                    catsare ={['Attention < median'],  [ 'Attention > median']};
                case 2
                    
                    plotd= useQdata;
                    
                    catsare ={['Lowest'],  ['...'], ['...'],['...'],['...'], ['Highest']};
                    %             case 3
                    
            end
            
            subplot(2,2,counter)
            %         title(titletouse); hold on
            % mean and stE for each category split
            mD = nanmean(plotd,1);
            stE=CousineauSEM(plotd);
            bH= bar(mD); hold on;
            %         if ianalysistype==1
            bH.FaceColor = col; %[.6 .5 .3];
            %         end
            
            errorbar(1:length(mD), mD,stE, 'Color', 'k', 'linewidth', 2,'LineStyle', 'none')
            ylim([.65 .95]);
            hold on;
            plot(xlim, [.5 .5], ['k-'])
            set(gca, 'XTickLabel', catsare, 'xtick', 1:length(mD), 'XTickLabelRotation', -15, 'Fontsize', 15);
            
            ylabel(ylabis);
            set(gcf, 'color', 'w');
            xlabel('Attention')
            
            
            %% test linear or quad to fit to data.
            
            %compare linear and quadratic GoF statistics:
            
            [F1,G1] = fit([1:size(useQdata,2)]', mD', 'poly1');
            [F2,G2] = fit([1:size(useQdata,2)]', mD', 'poly2');
            [F3,G3] = fit([1:size(useQdata,2)]', mD', 'poly3');           
            [F4,G4] = fit([1:size(useQdata,2)]', mD', 'exp1');
            [F5,G5] = fit([1:size(useQdata,2)]', mD', 'exp2');
            [F6,G6] = fit([1:size(useQdata,2)]', mD', 'power1');
            [F7,G7] = fit([1:size(useQdata,2)]', mD', 'power2');
%             
            %Goodness of fits:
            goodnessfitsare = [G1.adjrsquare, G2.adjrsquare,...
                G3.adjrsquare, G4.adjrsquare, ...
                G5.adjrsquare, G6.adjrsquare, G7.adjrsquare];
                
            
            %Output of Fits
            Fitsare  = {F1, F2, F3, F4, F5, F6, F7};
                
            
            %whichever was a better fit, now select that to plot on top of figure:
            [~,polyu] = max(goodnessfitsare);
            
            lgprint = ['best fit, adj. R^2 =' num2str(round(goodnessfitsare(polyu),2))];
            
            %fit and evaluate best polynomial:
            Fplot = Fitsare{polyu};
            hold on;
            pl=  plot(Fplot, [1:size(plotd,2)], mD);            
            pl(2).LineWidth = 3;            
             
            legend(pl(2), lgprint);
            
            %% replace axis labels:            
            ylabel(ylabis);
            xlabel('Attention')
        end
        counter=counter+1;
    end
    
    
    
    
    
    
end
cd(homedir); cd ../

cd(['Documents/Figures/Exp' num2str(iExp+1) ' Attention and Accuracy per participant'])

printname=['Fits by attention, both exps.'];


%         print('-dpng',printname );
