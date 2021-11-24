% calcAlpha_vs_Behaviour
%script called from Data_explore_B_EEG.

% Having already calculated alpha per trial, and split 
% various metrics by alpha quintiles 
% (previous job: calcAlpha_vs_Behaviour_table). 
% Here we will plot.

% which features to plot?

% % % now split by paper figures.
dbstop if error
job1.plotIndividualAlpha_byAttention=0;                 % plots attention ratings by alpha bins.
job1.plotIndividualAlpha_byXrating=1;                   % plot alpha bins vs xaxis ratings (visibility or confidence).
useAUC=0;                                               % will convert the above to AUC.

%note that there are different options for detection output, spec below:
job1.plotIndividualAlpha_byDetectionAccuracy =0;       % plots accuracy, HR, FAR, c, d', all together.     
combineExp1n2 =0; % to increase stat power.


%what to do with the output?
addGoF = 1;                                              % for plotting final fit:
normOUTPUT = 0;                                          % normalize per ppant, to the max  across alpha bins


% which plot level to show?
job1.plotPPantFX = 0;                                    % plot ppant level.
job1.plotGFX  =1;                                        % plot GFX
% job1.exporttoJasp =0;                                    % export table to JASP for further analysis.




statsOUTPUT=[]; % this collects as we go,
% and outputs to command window upon completion
%note that detection can be subdivided into objectively present or absent

%labels for plots:
TargClasses ={'Target-present', 'Target-absent', 'All trials'};

% labels for table columns:
ttypeCOLS = {'trgPres', 'trgAbs', 'trgAll'};

%% for combining across experiments:
% GFX_alphavsBEH = [GFX_alphavsBEH_exp1; GFX_alphavsBEH_exp2];
%%
% setpaths;
mydirs;
%
%% load channe data and specify channels to plot.
elocs=getelocs(2); 



nquintiles=5;
clf
% DETECTclass=4%1:3 % used to plot acc metrics.
    
for plotTargs= 1:3;%[1:3]; % 1:3 % present, absent all
    
    % 
    plotcombinedEXPs = []; % we'll concatenate and plot afterwards.

    for iExp=2%:2
    %% load GFX, will either plot ppant (each row), or groupdata.
    
    xlabis = xlabsare{iExp};
    
    cd(eegdir)
    cd( [xlabis ' alpha data']);
      
    load(['GFX_' xlabis '_ratingsandalpha.mat'], ...
        'GFX_behdata_byalphaTable',...
        'GFX_alphaTOPO',...
        'GFX_OccSpec', 'f_OccSpec'); % these two for plotting GFX of topo/spec
    %%
    nppants = size(GFX_behdata_byalphaTable,1);
    ppantFITS=nan(nppants ,2);
    GFX_alphavsBEH = nan(nppants , 2,nquintiles) ; % splitting into quintiles.

    dataTable = GFX_behdata_byalphaTable;
    
    plotOutput= nan(size(GFX_behdata_byalphaTable,1), 3, 5); % [ppants, targs, quintiles].
  
        %% which data from table?
        if job1.plotIndividualAlpha_byAttention;
            %extract data from table.            
            extractDV = 'Attention';
            yDV='z(Attention ratings)';
            compwas = 'Attention';
            usecol='r';
            
        elseif job1.plotIndividualAlpha_byXrating
            extractDV = xratings{iExp};
            compwas = extractDV;
            yDV=['z(' compwas ')'];
            usecol='b';
            
        
        elseif job1.plotIndividualAlpha_byDetectionAccuracy 
            
            usecol='k';
            switch DETECTclass
                case 1
                    extractDV='accuracy';
                    yDV = 'Accuracy';                   
                case 2
                    extractDV='Hitrate';
                    yDV = 'Hit rate';
                case 3
                    extractDV='FARate';
                    yDV = 'False alarm rate';
                case 4
                    extractDV='crit';
                    yDV='crit';
                case 5
                    extractDV='dprime';
                    yDV='d-prime';
                
            end
           
            
        end 
        
        %%
                    
            for iTarg = plotTargs
                ttype = ttypeCOLS{iTarg}; % string target type.
                
                %% note that ttype is only for subjective ratings.
                if ~job1.plotIndividualAlpha_byDetectionAccuracy                     
                colname = [extractDV '_' ttype];
                else
                    colname = [extractDV '_'];
                end
                
                if useAUC
                    colname= [colname '_AUC'];
                    
                    usecol='m';
                end
                
                %for each quantile
                for iq=1:5;
                    
                    plotOutput(:,iTarg,iq)= dataTable.([colname 'q' num2str(iq)]);
                end % end
            end % for each targ type.
            
    
            
            %% Which level to plot (ppant or group?).
            
    if job1.plotPPantFX % plot ppant level.
        plotPPant_alphaxBeh;        
    end
    
    %% after all ppants, plot GFX:
    if job1.plotGFX==1 && combineExp1n2~=1
        plotGFX_alphavsBeh;
    end
    
    if combineExp1n2
        % combine both into a matrix.
        plotcombinedEXPs = cat(1, plotcombinedEXPs , plotOutput);
    end
%     
% if job1.exporttoJasp==1
%     %% export for JASP stats?
%     % output as matrix:
%     pM= plotData; % the last plotted bardata.
%     TableforExport= splitvars(table(pM));
%     %name columns for stats:
%     Aquint = 1:size(pM,2);
%     
%     for icol = 1:size(TableforExport,2)
%         varN = ['AlphaQ' num2str(icol) ];
%         
%         TableforExport.Properties.VariableNames(icol)={varN};
%         
%     end
%     
%     %
%     cd('/Users/mdavidson/Desktop')
%     cd('JASP stats')
%     writetable(TableforExport,['GFX, alpha x ' ptitle ' ' xlabis '-Exp, ' TargClass '.csv'])
%     
%     %%
%     %also rotate for regr ession analysis
%     
%     conds = repmat(1:size(pM,2), size(pM,1), 1);
%     subjs = repmat([1:size(pM,1)]', 1, size(pM,2));
%     
%     %convert to columns
%     condsr = conds(:);
%     pMr= pM(:);
%     subjsr = subjs(:);
%     
%      TableforExport= splitvars(table([condsr,subjsr,pMr]));
%     varnames = {'alphabin', 'subjs', 'data'};
%     for icol = 1:size(TableforExport,2)        
%         TableforExport.Properties.VariableNames(icol)={varnames{icol}};        
%     end
%     cd('/Users/mdavidson/Desktop')
%     cd('JASP stats')
%     writetable(TableforExport,['GFX, alpha x ' ptitle ' ' xlabis '-Exp, ' TargClass '_rotated.csv'])
% end
% 
% % to concat across exps.
% if iExp==1
%     GFX_alphavsBEH_exp1 = GFX_alphavsBEH;
% else
%     GFX_alphavsBEH_exp2 = GFX_alphavsBEH;
% end


end % datatype


%% plot combined if requested.
if combineExp1n2
    plotOutput= plotcombinedEXPs;
    plotGFX_alphavsBeh
end


end
%% 
% print('-dpng', ['GFX, alpha x ' yDV ' ' xlabis '-Exp']);
print('-dpng', ['GFX, alpha x Obj BOTH-Exp']);
% saveas(gcf, ['GFX, alpha1 x ' yDV ' ' xlabis '-Exp'], 'png');

%% now display results in CW
% %
% clc;
% for it=3%[1:3]
%    %gather values for text output:
%     LRstat = statsOUTPUT.TARGtype(it).compBase_v_Linear.LRStat(2);
%     pval = statsOUTPUT.TARGtype(it).compBase_v_Linear.pValue(2);
%     
%     LRstat2 = statsOUTPUT.TARGtype(it).compBase_v_Quad.LRStat(2);
%     pval2 = statsOUTPUT.TARGtype(it).compBase_v_Quad.pValue(2);
%     
%     LRstatq= statsOUTPUT.TARGtype(it).compLinear_v_Quad.LRStat(2);  
%     pvalq= statsOUTPUT.TARGtype(it).compLinear_v_Quad.pValue(2);
% 
% disp([TargClasses{it} ]);
% disp(['>>>>>> base vs linear LRtest = ' num2str(LRstat) ',p=' num2str(pval)]); 
% disp(['>>>>>> base vs quad LRtest = ' num2str(LRstat2) ',p=' num2str(pval2)]); 
% disp(['>>>>>> linear vs quadratic LRtest = ' num2str(LRstatq) ',p=' num2str(pvalq)]); 
% end
% %%
%also evalue difference in regression slopes.
%%
% %targ present:
% try
% B1 = statsOUTPUT.TARGtype(1).MDLs.linearMDL.fixedEffects;
% SE1= statsOUTPUT.TARGtype(1).MDLs.linearMDL.Coefficients.SE;
% %targ absent:
% B2 = statsOUTPUT.TARGtype(2).MDLs.linearMDL.fixedEffects;
% SE2= statsOUTPUT.TARGtype(2).MDLs.linearMDL.Coefficients.SE;
% 
% %test for the equivalence of regression coefficients, as per:
% %Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. (1998)
% %. Using the correct statistical test for the equality of regression coefficients. 
% %Criminology, 36(4), 859-866.
% b1=B1(2);
% b2=B2(2);
% se1 = SE1(2);
% se2 = SE2(2);
% % Z= (b1-b2) / sqrt(se1^2 + se2^2)
% Z= (b1-b2) / sqrt(se1^2 + se2^2);
% pZ= normcdf(Z);
% disp(['equiv or regression coefs (present- absent): ' num2str(Z) ', p=' num2str(pZ)]);
% catch
% end