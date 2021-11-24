function sortEEG_percondition(homedir,~, ~)
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
    cd([xlabis ' targlocked data lowpfilt'])
    nfiles = dir([pwd filesep xlabis '*.mat']);
    
    %% load EEG
    for ippant = 1:length(nfiles)
        
%         cd(eegdir)
%         rej_trials = [];

        load(nfiles(ippant).name, 'TargLocked_EEGd',  'time_new');
        
        pEEGtouse = TargLocked_EEGd;
        %WAS  pEEGtouse = pEEG_stim_prepd;%
        
        [nchans, ~, ~]=size(TargLocked_EEGd);
        
        
        %% load beh
        %careful here - need to extract participant number (1 is missing).
        %%
        searchf= [xlabis '_Attn_participant%d'];
        str=nfiles(ippant).name;
        
        ppantnum = sscanf(str, searchf);
        cd(behdatadir)        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
%       

disp(['CHECK: pairing ' str '(EEG), with ' allf(1).name ',(Beh)']);
        cd(behdatadir)
        
        load(allf(1).name, 'p_table'); % load behavioural data in table form.
        % adjust trial info:
        p2= p_table;
        
        %what was the vis and attn rating per trial
        XY_ratings =table2array(p2(:,18:19));
        
        
        
        %%
        for  ioutcome = 1:4 %hit, miss, CR, FA
            
            id = find(p2.Outcome ==ioutcome); %HITs, M, CR, FA in data
            
            %for these trials, when was the stimulus presented?
            %note that the variable 'Targchip', denotes the frame.
            % SSVEP began 1s into trial. 10 Hz rate.
            % so index = 1 + nchip * .050 ; seconds
            
            
            
             
                
                VAN_EEGd= TargLocked_EEGd(:,:,id);
                VAN_BEHd=XY_ratings(id,:)   ;
                
            %save as appropr:
            switch ioutcome
                case 1 %HIT
                    HIT_BEHd= VAN_BEHd;
                    HIT_EEGd=VAN_EEGd;
                case 2 %MIss
                    MISS_BEHd= VAN_BEHd;
                    MISS_EEGd=VAN_EEGd;
                case 3 % CORrej
                    CR_BEHd= VAN_BEHd;
                    CR_EEGd=VAN_EEGd;
                case 4 % FA
                    FA_BEHd= VAN_BEHd;
                    FA_EEGd=VAN_EEGd;
            end
            
            
        end
        
        %% Sanity check
        
%                    h1=squeeze(mean(HIT_EEGd(1:32,:,:),3));
%                    m1=squeeze(mean(MISS_EEGd(1:32,:,:),3));
% %         %
%                    VAN = h1-m1;
%                     figure(2); clf
%                     plottopo(VAN, 'chanlocs', elocs32);
%                     figure(4); plot(time_new, VAN(32,:));
%         
        %% save data per ppant.
  cd(eegdir)
    cd([xlabis ' targlocked data lowpfilt'])
        save(nfiles(ippant).name, 'HIT_EEGd', 'HIT_BEHd', 'MISS_EEGd', 'MISS_BEHd', ...
            'FA_EEGd', 'FA_BEHd', 'CR_EEGd', 'CR_BEHd', 'time_new', '-append');
        
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Concatenate across participants
    %
    
    %% preallocate EEG
    if ippant == length(nfiles)
    nfiles = dir([pwd filesep xlabis '*.mat']);
    
    [GFX_HIT_EEGd,GFX_MISS_EEGd,GFX_CR_EEGd,GFX_FA_EEGd] = deal(zeros(length(nfiles),nchans, size(FA_EEGd,2)));
    
    
    
    for ippant = 1:length(nfiles)
        %change into correct directory:
        
        cd(eegdir)
        
        if lowpassfilter==1
            cd([ xlabis ' targlocked data lowpfilt'])
        else
            cd([ xlabis ' targlocked data'])
        end
        %%
        %load participant data:
        load(nfiles(ippant).name); % loads all EEG and also loads time_new (epochtime axis)
        
        for itype=1:4
            switch itype
                case 1
                    dIN = HIT_EEGd;
                case 2
                    dIN = MISS_EEGd;
                case 3
                    dIN = FA_EEGd;
                case 4
                    dIN = CR_EEGd;
            end
            
            %note that for this first test, we will simply average all
            mEEG= squeeze(mean(dIN,3));
            
            %store across participants.
            switch itype
                case 1
                    GFX_HIT_EEGd(ippant, :,:)=mEEG;
                case 2
                    GFX_MISS_EEGd(ippant, :,:)=mEEG;
                case 3
                    GFX_FA_EEGd(ippant, :,:)=mEEG;
                case 4
                    GFX_CR_EEGd(ippant, :,:)=mEEG;
            end
            
            
            
        end % itype
        
        
    end % ippants
    
    
    %save across all:
    save(['GFX_' num2str(xlabis) '_Attn_ERPs'], ...
        'GFX_HIT_EEGd', 'GFX_MISS_EEGd', 'GFX_FA_EEGd', 'GFX_CR_EEGd', 'time_new');
    end
end
end
