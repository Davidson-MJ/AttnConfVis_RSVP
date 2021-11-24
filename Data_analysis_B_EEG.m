% Data Exploration - EEG data
% Script to import, preprocess and analysis EEG data accompanying
% behavioural reports.
% for Behavioural analysis, see Behavioural_analysis_A_Beh.m




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % Analysis scripts: % % % % % %% % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % EEG preprocessing analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first extract and save participant data in .mat format, for analysis.

%> done for all.
job.extractEEGdata                            =0; % performed per participant, then saves in new fol.

%>done for (done for both exps.)
job.batch_preprocessEEGdata                   =0; % preprocess per participant (func for details) includes trial rejection

% alternate script, minimal preprocessing just ICA for eye (NY comments re: Eeglab).
job.EpochandICA_eeglab                        =0;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%analysis:
% %split by Hit, Miss, FA, CR, align to target ONSET, within the RSVP stream.
%all trials
job.define_targLockedEEG                      =0; % saves all trials, target locked Epochs.       
% by H,M,FA,CR
job.define_trialbasedonBEH                    =0; % separates EEG into trial conditions, for comparison (about target). (HIT, MISS, CR, FA)
%spit the above by quartile / tercile
job.define_trialbasedonAXISquantiles          =0; % median /quantile split of each category above, e.g. based on attention ratings.


%Calculate prestim per trial
% alpha power analyses. 
%calc alpha per trial, define alpha quintiles:
job.prestim_Alpha_power_pertrialFFT           = 0; % Calculate Alpha power per trial in the pre RSVP window (from complex values). 
%alpha phase analyses (complex value per trial)
job.Alphaphase_pertrial                       = 0; % Calculats Alpha phase per trial, at RSVP onset, and target 'onset'

% BEHAVIOUR by prestimulus alpha and phase:
%>>>>> PLOT! PLOT! these jobs also PLOT
job.calcBehaviour_byAlphaTerciles                = 0; % prepare participant tables of ALpha vs Behaviour for plotting.

% Inverted U / pre-stim alpha analysis:
%       Focusing on start of RSVP stream:
job.calc_ERPfeats_perQuintile                    =0; % for each quintile, we want measures of p1 and CPP.
job.calc_ERPfeats_perQuintileRestricted          =0; %, as above but sub-selecting by trial type (e.g. HITS only, for CPP).



%       Focusing on Target locked ERP 
job.calcBEHquintiles                             = 0; % saves with targetlocked data.

job.calc_TarglockedERPfeats_perQuintileRestricted = 0; % calculate TL erp, per Hit, Miss.



%%  NEW: Mediation analysus

job.performMediation                            =0;


%%%%%%%%%%%%%%%%%%%%%
%% Plots: the results of prestimulus alpha on various features:
% (MS plot)
job.plotinvertedUsummary_alpha                        =0; % stationary feats (peaks) etc (results of calc_TarglockedERPfeats_perQuintile)

%also works for the restricted case (e.g. Hits only).

job.plotinvertedUsummary_BEH                          =0; % ERP by Attn/Conf

%plot example singletrials 
job.plotExampletrialsbyPRESTIMULUSfeatures = 0;         % for presentations e.g. phase, power.

job.plotMSfigure_grandERP=0;


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% plotting:
% within these jobs,  whether to plot whole epoch or target locked
% ERPs is a variable set within each function.

job.plot_ppantbasedonBEH                       =0; % plots H,M,FA,COR, ppant level.
job.plot_ppantbasedonquantiles                 =0; % plots H,M,FA,COR, ppant level.

job.GFX_basedonQuantilesplits                  =1; % plots based on splits along each axis, for categories described above.

%% %%% %% set your home directory first!

% this determines which user is running the script:
% Query type depends on Mac vs Windows log.
mydirs;
%% %% clear to begin
clearvars -except job homedir datadir ;

%% %%%%%%% %%%%%%% %%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%                        Analysis
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%% perform data extraction per participant

Mycfg=[];
Mycfg.dsampleto            = 250;          % Hz, downsample the data to reduce comp and storage burden.
Mycfg.Epochlength_stim     = [.5, 3];   % epoch around both 'get ready' at trial start
Mycfg.Epochlength_resp     = [1.5, 1.5]; % epoch around response.



if job.extractEEGdata==1
    %import, epoch, demean, detrend, and downsample before saving per participant.

    % the old script, using field trip.
        extractEEGdataperppantandSave(homedir, datadir, Mycfg);    


end

if job.batch_preprocessEEGdata==1
    % artefact reject, reref, baseline subtract
    prepEEGdataperppantandSave(homedir, datadir, Mycfg);
end



if job.EpochandICA_eeglab==1 % new script, minimal preprocess.
    %import, epoch, ICA. Sanity check for NY comments re:eye mvmnts.    
    prep_run_ICA_eeglab;%(homedir, datadir);
    
end
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%% end preprocessing. Begin response based analysis:


% we will look at targetlocked ERPs, (hit-miss), based on target position in RSVP
% stream. Note that this uses a smaller epoch window.
if job.define_targLockedEEG==1       
    sortEEG_targetlocked(homedir,datadir, Mycfg); %
     % loads pEEG_stim_detr_dsamp_ICACOMPSr, 
     % saves TargLocked_EEGd
end

if job.define_trialbasedonBEH ==1
    
    %%%%%  Sorts target presentation by
    % HIT MISS, FA, CR. saves beh ratings per epoch (e.g. attention and
    % visibility/conf).
    sortEEG_percondition(homedir,datadir, Mycfg);
    % loads TargLocked_EEGd
    % saves HIT_EEGd, GFX_HIT_EEGd, GFX_MISS_EEGd,  etc
end



if job.define_trialbasedonAXISquantiles==1%     % as above, but sorting by value given by yaxis (attn ratings).
    sortEEG_percondition_byATTN(homedir,datadir, Mycfg);
    % loads HIT_EEGd
    % saves  GFX_HIT_EEGd_Attn; GFX_MISS_EEGd_XAXIS
 end


%% ALPHA power and phase


if job.prestim_Alpha_power_pertrialFFT              == 1 %Calculate Alpha power per trial in the pre RSVP window, using FFT
    calcAlphaP_pertrialFFT(homedir);
    % loads pEEG_stim_detr_dsamp
    % saves AlphaMean_FFT
end

if job.Alphaphase_pertrial              == 1 %Calculate Alpha phase per trial in the pre RSVP window.
    calcAlphaPHASE_pertrial(homedir);
    % loads pEEG_stim_detr_dsamp
    % saves complexEEG
end

%% %% P1 ERPs:
if job.calc_ERPfeats_perQuintile==1
    calcERPfeats_perQuintile(homedir) % uses predetermine alpha quintiles.    
    %saves f1_byalpha', 'Twof1_byalpha', 'P1_byalpha','P2_byalpha'
    
end
if job.calc_ERPfeats_perQuintileRestricted ==1 % as above, but only on subsets of trials (e.g. Hits only).
    calcERPfeats_perQuintile_restricted(homedir);
    %saves f1_bySDTandAlpha ...
end


if job.calcBEHquintiles==1          % used to align EEG with quintiles.
    calcBEH_quintiles(homedir);
end

if job.calc_TarglockedERPfeats_perQuintileRestricted==1
    calcTarglockedERPfeats_perQuintile_Restricted(homedir) % uses predetermine alpha quintiles.
    % loads TargLockedEEGData from [xlabis ' targlocked data lowpfilt'])
    % saves 'TL_bySDTandAlpha_erp' etc (used in plots for MS below).
end


%%

if job.performMediation                            ==1
    
    calc_EEG_BEH_mediation(homedir);
    %loads trial_table, 'AlphaMean_FFT, TargLocked_EEg
end





%% %%%%%%% %%%%%%% %%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
                        % PLOTTING
%%%%%%%%% %%%%%%% %%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
% BEHAVIOUR:
%% MAIN plot jobs (MS):


if job.plotinvertedUsummary_alpha   ==1    
   plotAlpha_vs_ERPfeats; % single exp at a time (CPP by alpha) (MS)
   % ERPs from : ['GFX_' xlabis '_ratingsandalpha.mat'], 
   % 'GFX_alphavsTL_erp';
   
   %alternatively (lite)
   plotAlpha_mainresults(homedir);  % MS plots of alpha x P1, and CPP (MS)
end

if job.plotinvertedUsummary_BEH   ==1    
   plotBEH_vs_ERPfeats;   % single Exp at a time (CPP by BEH quints). (MS)
   
   % pretty version for manusript:
   plotCPP_mainresults;
end


%%Plot Grand Average long ERP, for Manuscript:
if job.plotMSfigure_grandERP
GFX_wholetrialERPfigure(homedir,figuredir);
end


if job.calcBehaviour_byAlphaTerciles==1 
    % for each participants Alpha, split into terciles, and then collect
    % behavioural metrics (Attention, Conf, Vis, MCOG).
    
    calcAlpha_vs_Behaviour_table; % create across ppant table for easy plotting.
    
    prepareplotsAlpha_vs_Behaviour; % plots within this script (MS plot)
    
    
    
end


%%% other useful plots:

% ERPs by H, M, FA, CR:
if job.plot_ppantbasedonBEH    ==1    
    PFX_basedonCAT(homedir);
end
if job.plot_ppantbasedonquantiles==1    
    PFX_basedonquantiles(homedir);
end


if job.GFX_basedonBEH ==1
    % VAN (hit - miss).
    GFX_basedonCAT(homedir); % deprecated (high pfilt).    
end
