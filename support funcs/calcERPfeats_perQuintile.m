function  calcERPfeats_perQuintile(homedir) % uses predetermine alpha quintiles.
% need to adjust each epoch to use target locked responses.

% Focusing on the target-locked ERPs.
dbstop if error
%
% called from Data_explore_EEG
% MDavidson mjd070 dot gmail dot com
%%
%make sure ft functions are accesible.
% addpath(genpath('/Users/matthewdavidson/Documents/MATLAB/Matlab toolboxes/fieldtrip-20210920'));
%which experiment to analyze? Conf or PAS?

showPeakssaved=0; %this slows down the pipeline, but enables to debug whether correct peaks are being saved.
useChans = 29:32; % parieto-occipital subset (for alpha)

%%
%new addition.
MeasurePeakorAvg =1; % 0=raw peak amp, 1=measure avg
AtttopHalf=0; % median split, focus on only top half of attention ratings (sanity check)
%%
% for
for idatatype=2
    switch idatatype
        case 1
            xlabis = 'Conf';
        case 2
            xlabis = 'PAS';
    end
    
    
    dbstop if error
    
    %load channel data:
%     getelocs
    
    
    %set up directories
    % connect to behavioural data:
    cd(homedir)
    cd('Exp 2 and 3 mat files')
    behdatadir = pwd;
    
    cd(homedir);
    cd('Exp 2 and 3 processed EEG');
    eegdir = pwd;
    
    %% begin:
    %step through and perform preprocessing jobsL
    % 13:nfiles are the PAS subset
    
    
    %%
    
    
    %collect participant index of interest
    cd(eegdir)
    nfiles = dir([pwd filesep xlabis '*.mat']);
    
    %% load EEG
    %     for
    for ippant = 1:length(nfiles)
        
        
        disp([ '>>> loading participant ' num2str(ippant) ' exp ' num2str(idatatype)]);
        cd(eegdir)
        
        rej_trials =[];
        % load(nfiles(ippant).name, 'pEEG_stim_prepd', 'rej_trials');
        load(nfiles(ippant).name, 'pEEG_stim_detr_dsamp');
        
        
        %% save data per ppant, with ALPHA, to make it easy.
        
        cd(eegdir)
        cd([ xlabis ' alpha data'])
        searchf= [xlabis '_Attention_participant%d'];
        str=nfiles(ippant).name;
        ppantnum = sscanf(str, searchf);
        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
        %
        load(allf(1).name,  'AlphaMean_FFT', 'trial_table');
        %
        
        %% first correct for pre Stream baseline.
        
        % note that stream onset is at the 0.7s mark, AFTER
        onset = 0.7; % ms
        
        cfg=[];
        cfg.demean        = 'yes';
        cfg.baselinewindow = [onset-.1 onset]; % in seconds
        cfg.hpfilter = 'yes';
        cfg.hpfreq= 1 ; % .2 Hz
        cfg.hpfiltwintype = 'hann';
        cfg.lpfilter = 'yes';
        cfg.lpfreq= 25 ; % Hz
        cfg.lpfiltwintype = 'hann';
        
        % we can only look at top data if need be.
        if AtttopHalf==1
            attall=median(trial_table.AttResp);
            keeptrials = find(trial_table.AttResp > attall);
            cfg.trials = keeptrials;
        else
            keeptrials =1:length(pEEG_stim_detr_dsamp.trial);
            cfg.trials = 'all';
        end
        
        br_data = ft_preprocessing(cfg, pEEG_stim_detr_dsamp);
        nchans = 32;
        
        %% sanity check, plot the mean ERP, per ppant, POchans
        ParietoOccChan = 24:32; % as per MacDonald et al.
        POchans = ft_channelselection(ParietoOccChan, br_data);
        %
        cfg=[];
        cfg.channel = POchans;
        avERP = ft_timelockanalysis(cfg, br_data);
        % needs avERP as input >
        timeax=avERP.time;
%         showP1N1_SSVEPintrialERP;
        
       
        
        
        %take quintiles of alpha values,
        AlphaMean = squeeze(nanmean(AlphaMean_FFT(useChans,keeptrials),1));
        zAlpha = zscore(AlphaMean);
        zAttention = zscore([trial_table.AttResp]);
        zXrating= zscore([trial_table.Resp]);
        
        % also split by attention, or xratings.
        for iResp = 1:3
            disp([ '>>>  participant ' num2str(ippant) ' split ' num2str(iResp)]);
            
             % continue with analysis.
        % will perform the same, after first averaging within alpha quintiles.
        
        [P1_quintiles]=  deal(nan(nchans, 5));
        
            
            switch iResp
                case 1 % use alpha
                    [RespCats, lbins]= splitTrialsintoBins(zAlpha, 4); % n-1 bins.
                    
                    br_datatmp = br_data; % use all eeg trials (no restriction as below)
                case 2                                         
                    [RespCats, lbins] = splitTrialsintoBins(zAttention, 4); % n-1 bins.
                    
                    br_datatmp = br_data; % use all eeg trials (no restriction as below)
                case 3 % for conf/ vis ratings, restrict to Hits only.
                     restrange= trial_table.Outcome ==1;
                     %now take quintiles of this range:                     
                     tmp_Ad = [trial_table.Resp(restrange)];
                     
                     %also restrict EEG to this range
                     cfg=[];
                     cfg.trials= find(restrange);
                     cfg.keeptrials=1;
                     br_datatmp  = ft_timelockanalysis(cfg, br_data);
                     
                     
                     %>Continue by taking zscore of remaining relant trials, for
                     %alpha amplitudes.
                     tmp_Ad_z = zscore(tmp_Ad);
                     
                     [RespCats, lbins] = splitTrialsintoBins(tmp_Ad_z, 4); % n-1 bins.
            end
            
            
            %% beware, some ppants responses make it impossible to split into quintiles.
            if min(lbins)>0 % at least one trial per quintile.               
                for iquin = 1:5
                    trialsvector = RespCats{iquin};
                    
                    %% take average ERP over these channels.
                    
                    
                    %Keep trials separate.
                    cfg=[];
                    cfg.trials= trialsvector;
                    cfg.keeptrials=1;
                    allD = ft_timelockanalysis(cfg, br_datatmp);
                    %%
                    
                    
                    % now take P1,N1 per channel (POchans).
                    for ichan = 1:32 % don't waste time with frontal chans for vis response.
                        
                        % dive into a similar script as above...
                        
                        
                        tmpERP = squeeze(nanmean(allD.trial(:,ichan,:),1));
                        % find zero crossings:
                        
                        %appends to P1_byalpha within script.
                        measureBasetoPeak;
                        
                        if MeasurePeakorAvg==1
                            
                            %appends to P1_byalpha within script.
                            % also compute mean over 50ms window
                            if P1loct ~=301;
                                % take average over 20ms?
                                p1= mean(tmpERP(P1loct-10: P1loct+10)); % was 5
                                
                                P1_quintiles(ichan,iquin) = p1;
                            else
                                P1_quintiles(ichan,iquin) = nan;
                            end
                        end
                        
                        %%
                        
                        
                        %                  %show the stored peak(s).
                        if showPeakssaved
                            %%
                            P1t=P1_quintiles(ichan,iquin);
                            P2t=P1_quintiles(ichan,iquin);
                            
                            figure(10); clf
                            plot(timeax, tmpERP); hold on;
                            %addsearch window
                            %     Xs =[P1win(1), P1win(1), P1win(2), P1win(2)];
                            %     Ys = [-4 4 4 -4];
                            %     ph=patch(Xs, Ys, 'b');
                            %     Xs =[P1win(1), P1win(1), P1win(2), P1win(2)];
                            %
                            %     ph=patch(Xs, Ys, 'b');
                            %     ph.FaceAlpha = .1;
                            plot([onset onset], [ylim ], 'k:')
                            if ~MeasureProminence
                                % add p1 peak stored.
                                plot([timeax(P1loct) timeax(P1loct)], [0 tmpERP(P1loct)], 'b-', 'linew', 3)
                                % add p2 peak stored.
                                plot([timeax(P2loct) timeax(P2loct)], [0 tmpERP(P2loct)], 'b-','linew', 3)
                            else % show prominence.
                                % add p1 peak stored, and prominence.
                                
                                plot([timeax(P1loct) timeax(P1loct)], [tmpERP(P1loct)-P1t, tmpERP(P1loct)], 'b-', 'linew', 3)
                                % add p2 peak stored.
                                
                                plot([timeax(P2loct) timeax(P2loct)], [tmpERP(P2loct)-P2t, tmpERP(P2loct)], 'b-', 'linew', 3)
                            end
                            
                            title({['ippant ' (num2str(ippant)) ', peaks found (blue)'];...
                                ['chan ' num2str(ichan) ', quin ' num2str(iquin)]})
                            xlim([.5 1.5])
                            pause(.5);
                            
                        end
                        
                        
                        
                        
                    end % ichans
                    
                end % quin
            else
                
                disp(['>>>>>>> WARNING, EXP' num2str(idatatype) ' CATEGORY : ' num2str(iResp) 'pp ' num2str(ippant) ' skipped '])
                disp(['>>>>>>> insufficient trial range for ' allf(1).name ])
                disp(['>>>>>>> >>>>>>> >>>>>>> '])
            end % if length(quintiles) >5
            
        switch iResp
            case 1
                P1_byalpha=   P1_quintiles;
            case 2
                P1_byAttention=   P1_quintiles;
            case 3
                P1_byXrating=   P1_quintiles;
        end
                
        end % iResp

        
        save(allf(1).name, 'P1_byalpha', 'P1_byAttention', 'P1_byXrating', '-append');
        
    end%ppant
end % datatype
% fun