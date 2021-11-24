function sortEEG_targetlocked(homedir,~, ~)
% function to sort preprocessed data from fieldtrip, into categories for
% statistical comparison.


% Focusing on the target-locked ERPs.

% For ease of manipulation, the data is converted from ft structure to
% matlab n-d matrices, such that:

% [nchans, nsamps, ntrials]=  size(data);

% nchans = number of channels,
% nsamps = number of samples in epoch
% ntrials = number of trials (total).

% MDavidson mjd070 dot gmail dot com
%%

%which experiment to analyze? Conf or PAS?
for idatatype=1:2
    switch idatatype
        case 1
            xlabis = 'Conf';
        case 2
            xlabis = 'PAS';
    end
    lowpassfilter = 1; % change to 1 for low pass filter the data (below 10 Hz).
    
    dbstop if error
    
    %load channel data:
    getelocs
    
    
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
    for ippant = 1:length(nfiles)
        
        cd(eegdir)
        rej_trials = [];
%         load(nfiles(ippant).name, 'pEEG_stim_prepdICA')
        load(nfiles(ippant).name, 'pEEG_stim_detr_dsamp_ICACOMPSr');
        
      pEEG_touse= pEEG_stim_detr_dsamp_ICACOMPSr ; 
      
      
      nchans = length(pEEG_touse.label);
        nsamps = length(pEEG_touse.time{1});
        ntrials = length(pEEG_touse.trial);
        
        time_secs = pEEG_touse.time{1};
        %note that time_secs is actually zero'd at the offset of 'get ready' (so
        %300ms intostart of trial(!)
        
        %may want to lowp filter, but beware time-distortion:
        if lowpassfilter==1
            cfg=[];
            cfg.lpfilter= 'yes';
            cfg.lpfreq = 8 ; %Hz
            cfg.lpfilttype = 'firws';
            cfg.hpfilter= 'no';
%             cfg.hpfreq = .1 ; %Hz
%             cfg.hpfilttype = 'firws';
            
            dnow = ft_preprocessing(cfg, pEEG_touse);
            
            datatmp= cat(2, dnow.trial{:});
        else
            %no lowpfiltering.
%             cfg.hpfilter= 'yes';
%             cfg.hpfreq = .1 ; %Hz
%             cfg.hpfilttype = 'firws';
%             
%             dnow = ft_preprocessing(cfg, pEEG_touse);
%             
            datatmp= cat(2, pEEG_touse.trial{:});
            
        end
        
        %unpack: Channels x samples x trials
        EEGdata= reshape(datatmp, [nchans, nsamps,ntrials]);
        
        
        
        %% load beh
        %careful here - need to extract participant number (1 is missing).
        %%
        searchf= [xlabis '_Attention_participant%d'];
        str=nfiles(ippant).name;
        
        ppantnum = sscanf(str, searchf);
        cd(behdatadir)        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
%       

disp(['CHECK: pairing ' str '(EEG), with ' allf(1).name ',(Beh)']);
        cd(behdatadir)
        
        load(allf(1).name, 'p_table', 'SDTindex'); % load behavioural data in table form.
        % adjust trial info:
        p2= p_table;
        p2(rej_trials,:)=[];
        
        %% using new trial info, collect indices for all HIT, Miss, FA, CR:
        small_window = [-.2 1.5]; %seconds to realign epoch, about stim presentation
        % VAN, (visual awareness negativity), Hits vs Misses.
        % collect all like EEG index by when the stimulus was actually
        % present
        
        %note that original EEG time axis is given by:
        % EEGtime_ax = dnow.time{1}; time_secs
        %%
            
            %for these trials, when was the stimulus presented?
            %note that the variable 'Targchip', denotes the frame.
            % SSVEP began 1s into trial. 10 Hz rate.
            % so index = 1 + nchip * .050 ; seconds
            
            %% Epoch data:
            
            %What was target location per trial? targchip column = 11
            TargChips = table2array(p2(:,11));
            %what was the vis and attn rating per trial
            XY_ratings =table2array(p2(:,18:19));
            
            TargLocked_EEGd=zeros(35,426,936);
            

            %epoch
            for itrial = 1:size(EEGdata,3)
                
                targchip = TargChips(itrial);
                
                %define the '0' of stimulus onset:
                % At the 0.7 second mark in current epoch, RSVP stream begins
                % at 10 Hz. Each targ 'chip' represents the position in the stream.
                
                onset_s = .7 + (targchip-1)*.1;
                
                onset_epoch_s= onset_s+ small_window;
                %define epoch boundaries:
                tmp_window = dsearchn(time_secs', onset_epoch_s');
                
                %epoch all relevant data:
                
                tmp_EEGdata = squeeze(EEGdata(:,tmp_window(1):tmp_window(2),itrial));
                
                %baseline subtract (find onset t):
                time_new = [small_window(1):1/pEEG_touse.fsample:small_window(2)];
                f0 = dsearchn(time_new', [-.1,0]'); % 100ms around tg onset.
                
                %preallocate
                tmp_EEGdata2=zeros(size(tmp_EEGdata));
                
                for ichan = 1:size(tmp_EEGdata,1)
                    tmpch = squeeze(tmp_EEGdata(ichan, :));
                    basesub= mean(tmpch(1,f0(1):f0(2)));
                    tmp_EEGdata2(ichan,:) = tmpch - repmat(basesub, [1, length(tmpch)]);
                end
                
                
                %store
                TargLocked_EEGd(:, :, itrial)=tmp_EEGdata2;
                
            end
            
            
        
        
        %% Sanity check
%                
%                     figure(2); clf 
%                 h1=squeeze(mean(TargLocked_EEGd(1:32,:,p2.Outcome==1),3));
%                    m1=squeeze(mean(TargLocked_EEGd(1:32,:,p2.Outcome==2),3));
%         %
%                    VAN = h1-m1;
%                     figure(2); clf
%                     plottopo(VAN, 'chanlocs', elocs32);
%                    
% %                 
%         
        %% save data per ppant.
        
        cd(eegdir)
        if lowpassfilter==1
            cd([ xlabis ' targlocked data lowpfilt'])
        else
            cd([ xlabis ' targlocked data'])
        end
        %%
        save(allf(1).name, 'TargLocked_EEGd', 'time_new', 'XY_ratings', 'time_secs');
        
    end
    
   
end
end
