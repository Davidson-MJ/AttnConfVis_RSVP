function plotGFX_AttnxAccuracy_useboxplot(homedir,datadir)
cd(homedir)

%plots all in one figure, outputs in figure directory.
cd(homedir )
cd(datadir)
xlabsare={'Conf', 'PAS'};
titlelabsare={'Confidence', 'Visibility'};
xprints={'Confidence', 'visibility'};
AUCtype={'ROC', 'Precision-Recall'};
plotROCorRecall=1;

cmap = cbrewer('qual', 'Paired',12);
colormap(cmap)
plotnback=0; % change to 1 to plot equivalent data using attention value from n-trials back.
%%


figsize = [ 0 0 1 .4];

%have already done the stats:


figure(1);    clf;
set(gcf, 'Units', 'normalized','Position',figsize, 'color', 'w');
counter=1;

%grab colours for consistency between plots.
RespColours = brewCOLOURS;
%

%
%analysis types, differ in what to plot:
%     1) Accuracy
%     2) AUC
%     3) HRate
%     4) FArate
%     5) Dprime
%     6) criterion

ylabsare = {'Accuracy', 'AUROC2', 'HRate', 'FArate', 'd''', 'Crit.'};

for ianalysistype=1:2
    
    for iExp=1:2  % for both experiments
        
        xlabis = xlabsare{iExp};
        
        cd(datadir)
        
        loadvarname= [ xlabis '_Attn_Allpp.mat'];
        load(loadvarname, 'all_AUC*', 'all_Accuracy*', '*bymulti', 'nshift');
        switch ianalysistype
            case 1  % use Accuracy measures
                useQdata= all_Accuracybymultiforpred;
%                 SCcol='k';
                SCcol=squeeze(RespColours(iExp,1,:))';
            case 2 % use AUC measures
                useQdata= all_AUCbymultiforpred;
                SCcol=squeeze(RespColours(iExp,1,:))';
            case 3 % HR
                useQdata= all_HRbymulti;
                SCcol=squeeze(RespColours(iExp,1,:))';
            case 4 % FA
                useQdata= all_FArbymulti;
                
                
            case 5 % Dpr
                useQdata= all_Dprimebymulti;
                
            case 6 % Crit
                useQdata= all_Critbymulti;
                
        end
        %%
      
        ylabis = ylabsare{ianalysistype};
        titletouse = [ylabis ' by attn'];
        
        
        subplot(1,4,counter)
        %         title(titletouse); hold on
        % mean and stE for each category split
        
        remm = isinf(useQdata);
        useQdata(remm) = nan;
        
        plotd= useQdata;
        %%
        params.cols = flipud(squeeze(RespColours(3,:,:)));
%         params.cols = [0 0 1];
        params.plotScatter=1;
        params.scatSize = 40;
        params.scatCol =squeeze(SCcol);
        %         clf
        h= box_and_scatter(plotd, params);
        
        %%
        
        fontsize= 10;
        ylim([.45 1.05]);
        hold on;
        plot(xlim, [.5 .5], ['k--'])
        
        set(gca, 'XTickLabel', [], 'xtick', 1:size(plotd,2),...
            'Fontsize', fontsize*1.5);
        
        hold on;
        %add sig line and asterisk? (if sig, confirm in JASP).
        
        
        %         plot([1, 5], [4.01,4.01], ['k-'], 'linew', 1.1)
        %         tc=text(3, 4.01, '*','fontsize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
        plot([1, 5], [1.01,1.01], ['k-'], 'linew', 1.1)
        tc=text(3, 1.01, '***','fontsize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        ylabel(ylabis);
        set(gcf, 'color', 'w');
        %             xlabel('Attention')
        
        
        %% test linear or quad to fit to data.
        %             mD= squeeze(nanmean(plotd,1));
        %             plotGoF(mD);
        %
        
        
        %% replace axis labels:
        ylabel(ylabis);
        title(['Experiment ' num2str(iExp)], 'fontweight', 'normal')
        
        counter=counter+1;
        
        c=colorbar('Ticks', [0,1], 'Ticklabels', {'Low', 'High'}, 'location', 'southoutside');
        ylabel(c, 'Attention rating', 'fontsize', fontsize*1.5)
        colormap(flipud(squeeze(RespColours(3,:,:))));
        
        set(gca, 'fontsize', fontsize*1.5)
        %     axis tight
    end
    %% export for stats?
    %     output as matrix:
    %     remove infinite values
    remm = isinf(useQdata);
    useQdata(remm) = nan;
    %     TableforExport= splitvars(table(useQdata));
    TableforExport= array2table(useQdata);
    
    for icol = 1:size(TableforExport,2)
        varN = ['AttQ' num2str(icol) ];
        TableforExport.Properties.VariableNames(icol)={varN};
    end
    cd(homedir)
    cd('JASP stats')
    writetable(TableforExport, ['GFX ' xlabis '-Exp,' ylabis 'by_Attention.csv'])
    %
    %
    %
    
    
    
    
    
    
    
end
cd(homedir); cd ../

cd(['Documents/Figures/Exp' num2str(iExp+1) ' Attention and Accuracy per participant'])

% printname=['Fits by attention, both exps.'];


%         print('-dpng',printname );
