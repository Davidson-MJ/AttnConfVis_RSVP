% Behavioural Data  Analysis
% Script to extract the data per participant (mouse clicks). 
% Then many separate behavioural analysis and plotting jobs, listed below.


%set dirs
mydirs;
%% %% clear to begin 
clearvars -except homedir datadir ;

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % Analysis scripts: % % % % % %% % % % % % % % 
%% Parameter analysis
job.dispMeanContrast =0;
job.dispTargPosResults = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % CLICK analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first extract and save participant data for all clicks in XY for scatter
% plots. Also save the distribution across X and Y axes for comparison.

job.extractMouseClickdata                            =1;% performed per participant
job.combineacrossparticipants                        =1;% collates across (no averaging).
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Calibration x attention analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does attention impact Confidence calibration (AUC /subjective)?
job.calcGrand_Accuracy_and_Calib                   =0;

job.calcPpant_AttnCalib_accuracy                   =0;% Calculates various metrics by attention quintile.
job.calcPpant_AttnCalib_Mcog                       =0;% calculates AUC, by attention split (median or tercile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINING ACROSS Participants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
job.combine_Attncalib                                =0; % combines all at group level.
job.combine_Xcalib                                   =0;% combines all at group (after x-axis splits).




% % % % % % % % % % % % Plotting: % % % % % % % % % % % % % % % % % %

% Plots scatter data per SDT metric, along with collapsed counts of
% responses along each axis.
%%%%%%%%%%%%%%%%%%%%%%%
% clicks
%%%%%%%%%%%%%%%%%%%%%%%
job.plot_Clicksperppant                             =1;% output in figdir.
job.plot_Clicksacross                               =0; %Oldver also plots RTs across.
job.plot_Clicksacross_final                         =1; %new ver, condensed and pretty.
job.plot_Clicks_perppant_comb                       =1; % supplementary figure.

% job.plot_ClicksCorr_final                           =1; % across participants, plot the correlation between measures.
%%%%%%%%%%%%%%%%%%%%%%%
% y axis
%%%%%%%%%%%%%%%%%%%%%%%
job.plot_Attn_x_Accuracyperppant                    =0;
job.plot_AttnCalibperppant                          =0;

plotROCorRecall                                     =1; % [1,2] = ROC,Precision-Recall
%%%%%%%%%%%%%%%%%%%%%%%
% across participants:
%%%%%%%%%%%%%%%%%%%%%%%
%Y-axis split
job.plot_Attn_Accuracyacross                        =0; % This plots bar graphs, and fits (polyfit).
job.plot_Attn_Accuracyacross_wboxplot               =1; % Accuracy as well as other metrics (dprime ,crit, FAR, etc).  
job.plot_Attn_Accuracyacross_wraincloud             =0; % as above, but uses raincloud (horizontal) on tercile split.
job.plot_AttnCalibacross                            =0; % will also plot combined version (comparing ROC to P-Recall).

%% %%%%%%% %%%%%%% %%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%                        Analysis 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%% parameter analysis

if job.dispMeanContrast ==1
    dispContrastresults(homedir);
end
   
if job.dispTargPosResults==1
    dispTargPosResults_onAccuracy(homedir);
    dispTargPosResults_onRatings(homedir);
end



%% perform data extraction per participant
if job.extractMouseClickdata==1    
    % also saves per participant.   
    extractMouseClicksperppantandSave(homedir, datadir);       
end

%% collate click data across participants
if job.combineacrossparticipants==1
    combineMouseClicksacrossppantandSave(homedir,datadir);
end

%% Accuracy and calibration, without segmentation
if job.calcGrand_Accuracy_and_Calib==1
    % calcs per ppant, also concatenates and plots GFX >>>>>>>>>>>>>>>>>>>>>>>(MS plot)
    % using 'plotBOXES_AccandROC.m' 
    calc_ppantAcc_ROC(homedir,datadir);
end

%% confidence calibration (segmenting by attention (y-axis bins)

if job.calcPpant_AttnCalib_accuracy==1
   calc_ppantAttnCalib_objectiveAccuracy(homedir, datadir) 
end

if job.calcPpant_AttnCalib_Mcog==1
    calc_ppantAttnCalib_ROC(homedir, datadir)    
end


 %% COMBINING across particiapnts
if job.combine_Attncalib==1
    
    combineAttnCalibacrossppantandSave(homedir,datadir);
end

%% %%%%%%% %%%%%%% %%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%                        Plotting 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% 

%% % % % click data:

if job.plot_Clicksperppant==1    
    plotPP_ScatterandXY(homedir,datadir)
end

if job.plot_Clicksacross==1
    plotGFX_ScatterandXY(homedir,datadir);
end
if job.plot_Clicksacross_final==1
    % >>>>>>>>>>>>>>>>>>>>>>>(MS plot)
    plotGFX_ScatterandXY_final(homedir,datadir);
end

if plot_Clicks_perppant_comb
    % >>>>>>>>>>>>>>>>>>>>>>>(MS plot)
    plotPP_ScatterandXY_suppfig(homedir, datadir);
end
%% % % % Y-axis SPLITS
%%%%%%%%%%%%%%%%%%%%%%%
% Y-axis split ppant level
%%%%%%%%%%%%%%%%%%%%%%%
if job.plot_AttnCalibperppant==1
    plotPP_Attncalibration(homedir,datadir, plotROCorRecall);   
end
if job.plot_Attn_x_Accuracyperppant    ==1
    plotPP_Attn_x_accuracy(homedir,datadir)    
end
%%%%%%%%%%%%%%%%%%%%%%%
%  Y-AXIS group data:
%%%%%%%%%%%%%%%%%%%%%%%

if job.plot_Attn_Accuracyacross==1
   plotGFX_AttnxAccuracy(homedir,datadir); 
end
% can also plot using boxplots (prettier)
if job.plot_Attn_Accuracyacross_wboxplot==1
    % >>>>>>>>>>>>>>>>>>>>>>>(MS plot)
plotGFX_AttnxAccuracy_useboxplot(homedir, datadir);
end
%can also plot horizontal version, using rainclouds:
if job.plot_Attn_Accuracyacross_wraincloud==1
plotGFX_AttnxAccuracy_raincloudver(homedir, datadir);
end

