% function     calc_EEG_BEH_mediation(homedir);


%called from Data_explore_B_EEG.
% this script will prepare and perform a mediation analysis,
%investigating
%- prestimulus alpha (power)
%- target locked ERP (amplitude)
%- Subjective criteria (attention-confidence-visibility.


job1.wrangle_Across_participants =1;
%step 1. Data wrangle. For the analysis, we need 3 vectors (time-series).
% can be matrices.


job1.crunchMediation = 1;





%% set up params
%for a single subject, these can be single trial
%alpha,
% TL amp
% ratings.
ParietoOccChan = 29:32; %
CPchans= [15:17, 20:22];
elocs=getelocs(2);

%set up directories
% connect to behavioural data:
cd(homedir)
cd('Exp 2 and 3 mat files')
behdatadir = pwd;

cd(homedir);
cd('Exp 2 and 3 processed EEG');
eegdir = pwd;
%set up plot output.

xratings = {'Confidence' , 'Visibility'};
xlabsare = {'Conf', 'PAS'};

%deonte mediation trial subset.
trialtypes = {'target present', 'Hits only', 'all trials'};

usetrials=2;


if job1.wrangle_Across_participants ==1
    %%
    [X, Y, M, M2] = deal([]); % parameters for our mediation model.
    [Xmat, Ymat, Mmat, M2mat] = deal(zeros(936,12)); % parameters for our mediation model.
    %which trial types to retain? (see above).
    
    
    for idatatype=1:2
        xlabis = xlabsare{idatatype};
        
        % concat across ppants.
        
        %% for each participant,
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        
        for ippant =1:length(nfiles)
            %%
            cd(loaddir);
            %plot attention and state ratings. The process is:
            load(nfiles(ippant).name, 'trial_table', 'AlphaMean_FFT');%
            %also gather TL eeg.
            cd ../
            cd([ xlabis ' targlocked data lowpfilt']);
            load(nfiles(ippant).name, 'TargLocked_EEGd', 'time_new');
            
            %% collect data for this participant
            
            %first, which subset of trials shall we use?
            switch usetrials
                case 1 %use all target present trials.
                    allH = find(trial_table.TargPresent==1); % Hits
                case 2 %use hits only.
                    allH = find(trial_table.Outcome==1); % Hits
                case 3 %use all trials in mediation analysis.
                    allH=1:size(trial_table,1); % all trials
            end
            
            
            %Predictor (Alpha).
            X_tmp = zscore(squeeze(nanmean(AlphaMean_FFT(29:32,:),1))); % prestimulus alpha power.
            
            %Outcome
            %ERP for now. target locked EEG amp, sigle trials, averaged 250-550ms after onset.
            
            tmp = squeeze(nanmean(TargLocked_EEGd(CPchans, : ,:),1));
            timeav = dsearchn(time_new', [.25 .55]');
            %
            Y_tmp= nanmean(tmp(timeav(1):timeav(2),:),1);
            
            % Mediators or alternatively, use other behavioural resp as mediator.
            M_tmp = zscore(trial_table.Resp); % x axis
            M2_tmp= zscore(trial_table.AttResp);  % y axis (attention)
            
            
            
            
            
            %append to group data, use cell array, given possibility of
            %different trial counts.
            
            X{ippant} = X_tmp(allH)';
            Y{ippant} = Y_tmp(allH)';
            M{ippant} = M_tmp(allH);
            M2{ippant}= M2_tmp(allH);
            
            if usetrials== 3 % we have equal trial nums. so can also store
                %as matrix (used in debugging).
                Xmat(:,ippant)=X_tmp(allH);
                Ymat(:,ippant)=Y_tmp(allH);
                Mmat(:,ippant)=M_tmp(allH);
                M2mat(:,ippant)=M2_tmp(allH);
            end
            
        end % ppants.
        
        cd(eegdir);
        try cd([xlabis ' mediation results'])
        catch
            mkdir([xlabis ' mediation results'])
            cd([xlabis ' mediation results'])
        end
        %% just use cells for each type, easier.
        X_alpha = X;
        Y_TLERP = Y;
        M1_resp = M;
        M2_attresp = M2;
        %%
        usedtrials = trialtypes{usetrials};
        save(['GFX_' xlabis '_for_mediation'], 'X_alpha','Y_TLERP', 'M1_resp', 'M2_attresp', 'usedtrials');
    end % datatype
    
end % job1

if job1.crunchMediation
    %% perform mediation analysis.
    %% since we are limited to linear relationships, consider different examples
    %X = predictor, Y= outcome, M = mediator(s).
    %1. X = Alpha, Y=ERP, M1, M2 = resps.
    
    %% set mediation type
    %1)  multiple mediator example.
    
    %      this will use X= alpha, Y= ERP, mediators are subjective ratings.
    
    %2) %simple mediator model, subj ratings as mediator on each other
    
    % this will use X= ERP, Y=XRating, M1 = Attention
    
    
    %3) %simple mediator model, complement to the above.
    
    % this will use X= ERP, Y=Attention, M1 = Xrating
    
    %%
    
    
    
    %%
    close all;
    %     mmodel = [x,y, m, m1];
    
    
    
    
    for idatatype=2
        
        
        cd(eegdir)
        xlabis = xlabsare{idatatype};
        cd(eegdir)
        cd([xlabis ' mediation results'])
        load(['GFX_' xlabis '_for_mediation'])
        
        
        for MediationType = 1
            switch MediationType
              
               
                case 1 % simple model:
                    %%
                    namesAre = {'CPP', xratings{idatatype}, 'Attention'};
                    
                    
                    % perform mediation, save figure output.
                    %                 [paths,stats] = mediation(Y_TLERP',M1_resp', M2_attresp', 'boot', 'bootsamples', 1000, 'names', namesAre);
                    [paths,stats] = mediation(Y_TLERP',M1_resp', M2_attresp','covs', X_alpha,'verbose', 'boot', 'bootsamples', 1000, 'names', namesAre);
                    %add colours for plots
                    stats.plotCols = [.7 .7 .7; 0, 0, .9; .9,0,0];
                    my_mediation_path_diagram(stats);
                    % add titles to figures
                    %                mediation_path_diagram(stats, 'plots', 'verbose')
               
                case 2 % flip of above
                    %%
                    namesAre = {'CPP', 'Attention', xratings{idatatype}};
                    
                    % perform mediation, save figure output.
                    [paths,stats] = mediation(Y_TLERP',M2_attresp', M1_resp','verbose', 'boot',  'names', namesAre);
                    %add colours for plots
                    stats.plotCols = [.7 .7 .7; .9,0,0;0, 0, .9;0, 0, .9];
                    stats.Cov= 1; % covariate plotting option, shifts the path diagram.
                    my_mediation_path_diagram(stats);
                    % add titles to figures
            end
            num=gcf;
            figure(gcf);
            set(gcf,  'Position',[num.Position+(num.Number*10)])
            sgtitle(['Trials: ' usedtrials])
            set(gcf, 'color', 'w')
            %% perform mediation, save figure output.
            
        end %datatype
    end
end % job mediation