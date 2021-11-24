function     plotPP_Attncalibration(homedir,datadir,plotROCorRecall)
cd(homedir)

xlabsare = {'Conf', 'PAS'};

% plotROCorRecall=; % 1 to use ROC characteristic, or 2 for Precision-Recall space for calculating AUC
                    % note that P-R space might be best for uneven
                    % quantiles, check.

fontsize=15;

for iExp=1:2
    xlabis = xlabsare{iExp};
	cd(datadir)    
    %how many data types/ppants this folder.
%     nppants = length(dir([pwd filesep xlabis '_Attn_p*']));        
    
    if iExp==1
        loadnppants=[1:12];
    else
        loadnppants=[2:10];
    end
    
%     dimsare = ['nppants', 'X_Attn', 'SDTm(H,M,FA,CR)', 'axis'];
    for ippant=loadnppants
        cd(datadir)
        % open p_table
        loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
        
        %load vars of interest.        
        load(loadvarname, 'RespxAttn*', 'Qp', 'attnMedian');
        
%% and plot as we go. 
% close all
 figure(ippant); clf; 
 %set half screen figure size, and white colour:
 set(gcf, 'Units', 'normalized','Position', [.5 0 .5 1], 'color', 'w');
 %%
 linest='-';
 
 for iRespCat=1:2 % median split and quantile split
        
     lgprint={}; % reset legend characters.
     counter=1;  % for legend entries.
     switch iRespCat
         case 1
              catsare = {'Attn < Median', 'Attn > Median'};
              
              if plotROCorRecall==1
              plotRespData = RespxAttn_mediansplit_AUC;
              plotRespDataxy = RespxAttn_mediansplit_AUCxy;
              
              %for type ROC data:
              ylabpl=['True positive rate'];
              xlabpl=['False positive rate'];
              
              lgloc = 'SouthEast';
              else
                  plotRespData = RespxAttn_mediansplit_AUC_PRver;
              plotRespDataxy = RespxAttn_mediansplit_AUCxy_PRver;
              
              
              %for type Precision recall data:
              ylabpl=['Precision'];
              xlabpl=['Recall'];
              lgloc = 'SouthWest';
              end
         case 2
               catsare = {'Attn Q1', 'Attn Q2', 'Attn Q3', 'Attn Q4'};
               
                   if plotROCorRecall==1
               plotRespData = RespxAttn_Qtlsplit_AUC;   
               plotRespDataxy = RespxAttn_Qtlsplit_AUCxy;
                   else
                       plotRespData = RespxAttn_Qtlsplit_AUC_PRver;   
               plotRespDataxy = RespxAttn_Qtlsplit_AUCxy_PRver;
                   end
               
     end
     pH=[];
    for icat = 1:length(plotRespData) %for each iteration this category.
     
        %data to plot:
        x= plotRespDataxy(icat).xdata;
        y= plotRespDataxy(icat).ydata;
        AUC=plotRespData(icat);
        
        %plot
        subplot(2,2,iRespCat)
        hold on
       pH(icat)= plot(x, y, 'linew', 2, 'linestyle', [linest]); hold on;
        lgprint{counter} = [catsare{icat} ', AUC = ' num2str(round(AUC,2))];
        counter=counter+1;
    end
    

           %legend on subplot properties
    plot(nan(1,1), 'color', 'w')
    if iRespCat==1
        % add the median value to our legend
        lgprint{counter} = ['(median Attn=' num2str(attnMedian) ')'];
    else
        % add the quartile splits
        lgprint{counter} = ['(quartiles Attn=0 ' num2str(round(Qp)) ')'];
    end              
    %legend on subplot properties
    lg=legend(lgprint);
       xlabel(xlabpl); ylabel(ylabpl);
       set(gca, 'fontsize', fontsize)
       %note that recall ~= FPR:
        % prec = tp/(tp+fp);
        % FP =  fp/(tn+fp);
    
        
    
    %for type ROC data:
    xlabel(xlabpl); ylabel(ylabpl);
    title(['Attention and ' xlabis ' calibration']);
    
    set(lg, 'Location', lgloc)
    set(gcf,'color', 'w');
    set(gca,'fontsize', fontsize)
        
    
    %%
    %also Bar plots for quick Reference.
    subplot(2,2,iRespCat+2);
   
        bH=bar(plotRespData);
        bH.FaceColor='flat'; % need to define so that we use CData property.
        
        %recolour individual bars:
        for ibar=1:length(plotRespData)
            colis=get(pH(ibar), 'color');
            bH.CData(ibar,:) =colis;
        end
    %%


    %what were the previous object colours?
    

    ylabel('AUC'); ylim([0 1]);
    title('AUC comparison');
    set(gca, 'XTickLabel', catsare, 'XTickLabelRotation', -15, 'fontsize', 15)
 
 end
 
 % now print to figdir:
 cd(homedir); cd ../
        cd(['Documents/Figures/Exp' num2str(iExp+1) ' Attention Calibration per participant'])
        
        printname=['Attention Calibration Participant ' num2str(ippant) '_' xlabis ];
        if plotROCorRecall==2
            printname=[printname ' PrecRecall ver'];
        end
        
        print('-dpng',[printname]);
 
    end
end