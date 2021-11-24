function     plotPP_Attn_x_accuracy(homedir,datadir)

   cd(homedir)

xlabsare = {'Conf', 'PAS'};
xprint={'Confidence', 'visibility'};
% plotROCorRecall=; % 1 to use ROC characteristic, or 2 for Precision-Recall space for calculating AUC
                    % note that P-R space might be best for uneven
                    % quantiles, check.

fontsize=15;

for iExp=1:2
    xlabis = xlabsare{iExp};
	cd(datadir)    
    %how many data types/ppants this folder.
    nppants = length(dir([pwd filesep xlabis '_Attn_p*']));        
    
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
        load(loadvarname, 'RespxAttn*', 'Y_Qp', 'Y_median');
        
%% and plot as we go. 
% close all
 figure(ippant); clf; 
 %set half screen figure size, and white colour:
 set(gcf, 'Units', 'normalized','Position', [.5 .5 .5 .5], 'color', 'w');
 %%
 linest='-';
 
 for iRespCat=1:2 % median split and quantile split
        
     lgprint={}; % reset legend characters.
     counter=1;  % for legend entries.
     switch iRespCat
         case 1
              catsare = {['Attention < Median'],  [ 'Attention > Median']};
              plotRespData = RespxAttn_mediansplit_Accuracy;                                      
         case 2
             
               catsare = {[  'Attention Q1'], [ 'Attention Q2'], [ 'Attention Q3'], [ 'Attention Q4']};                   
               plotRespData = RespxAttn_Qtlsplit_Accuracy;                  
     end
     pH=[];
   
%plotting
        % Bar plots for quick Reference.
        subplot(1,2,iRespCat);
                  
    %legend on subplot properties
%     plot(nan(1,1), 'color', 'w'); hold on;
    
    title({['Discrimination performance, by Attention,'];[ xlabsare{iExp} ' exp, p' num2str(ippant)]});
    hold on
        bH=bar(plotRespData);
        bH.FaceColor = [.6 .5 .3];
                    
    ylabel({['Accuracy'];['(TP+TN)/(P+N)']}); ylim([0 1]);    
    set(gca, 'XTickLabel', catsare, 'xtick' , 1:length(plotRespData), 'XTickLabelRotation', -15, 'fontsize', 15)
    
    if iRespCat==1
        % add the median value to our legend
       legend([  'Attention median=' num2str(Y_median) ]);
    else
        % add the quartile splits
    legend(['Attention quartiles=' num2str(round(Y_Qp))]);
    end
        
 end
 
 % now print to figdir:
 cd(homedir); cd ../
        cd(['Documents/Figures/Exp' num2str(iExp+1) ' Attention and Accuracy per participant'])
        
        printname=['Accuracy by attention, Participant ' num2str(ippant) ' ' xlabis ' exp'];
       
        
        print('-dpng',printname);
 
    end
end
end