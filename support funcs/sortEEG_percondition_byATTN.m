function sortEEG_percondition_byATTN(homedir,~,~)
% function to sort preprocessed data from fieldtrip, into categories for
% statistical comparison.

% for ease of manipulation, the data is converted from ft structure to
% matlab n-d matrices, such that:

% [nchans, nsamps, ntrials]=  size(data);

% nchans = number of channels,
% nsamps = number of samples in epoch
% ntrials = number of trials (total).

% MDavidson mjd070 dot gmail dot com
%%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


uselowpfilt=1; % which epoch type to use - trying to replicate Jmac study,
for idatatype =1:2
    switch idatatype
        case 1
            xlabis = 'Conf';
        case 2
            xlabis = 'PAS';
    end
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
    %We can use the already sorted EEG data per participant:
    
    
    %collect participant index of interest
    cd(eegdir)
    nfiles = dir([pwd filesep xlabis '*.mat']);
    
    %% load EEG
    %load data per participant
    cd(eegdir)
    if uselowpfilt==1
    cd([ xlabis ' VAN data lowpfilt'])
    else
        cd([ xlabis ' VAN data'])
    end
    %%
    
    allppants = dir([pwd filesep '*_Attn_participant*']);
    
    for ippant = 1:length(allppants)
        cd(eegdir)
        if uselowpfilt==1
            cd([ xlabis ' VAN data lowpfilt'])
        else
            cd([ xlabis ' VAN data'])
        end
        load(allppants(ippant).name)
        
        nchans = size(MISS_EEGd,1);
        
        %now for each type of response, sort by attention (top/ bottom / quartiles)
        
        
        [outgoingEEG_Attnsplit, outgoingEEG_XAXISsplit ]= deal(zeros(4,4,nchans,size(MISS_EEGd,2))); % SDT x tercile, X chans x samps.
        
        [outgoingBEH_Attnsplit , outgoingBEH_XAXISsplit] = deal([]); % use structur since different trial counts.
        for id = 1:4
            
            switch id
                case 1%HITS, M, FA, CR
                    behData = HIT_BEHd;
                    eegd = HIT_EEGd;
                case 2
                    behData = MISS_BEHd;
                    eegd = MISS_EEGd;
                case 3
                    behData = FA_BEHd;
                    eegd = FA_EEGd;
                case 4
                    behData = CR_BEHd;
                    eegd = CR_EEGd;
            end
            
            %ok, sort trials based on quartile split along y-dimension (ATTENTION).
            %%
            for IAXISSPLIT = 1:2
                
                switch IAXISSPLIT
                    case 1      %Attnratings
                        useDATA= behData(:,2);
                    case 2 % XAXIS
                        useDATA= behData(:,1);
                        %take absolute of x-axis data, as the Misses and CR category
                        %will be negative, with larger numebrs indicating higher
                        %confidence.
                        %               useData = abs(useDATA);
                end
                
                    Qps = quantile(useDATA, [.25 .5 .75]);


                if length(Qps) ==3
                    A_split = [];
                    
                    %QP1
                    ind1 = find(useDATA <= Qps(1));
                    
                    ind2a = find(useDATA > Qps(1));
                    ind2b = find(useDATA <= Qps(2));
                    %QP2
                    ind2 = intersect(ind2a, ind2b);
                    
                    %QP3
                    ind3a = find(useDATA > Qps(2));
                    ind3b = find(useDATA <= Qps(3));
                    ind3 = intersect(ind3a, ind3b);
                    
                    %QP4
                    ind4 = find(useDATA> Qps(3));
                    
                    
                    %sanity check.
                    mBar = [mean(useDATA(ind1)), mean(useDATA(ind2)),mean(useDATA(ind3)),  mean(useDATA(ind4))];
                    %     bar(mBar); shg
                    
                    
                    
                    %% store index of relevant splits for saving eeg :
                    switch IAXISSPLIT
                        case 1
                            outgoingEEG_Attnsplit(id, 1, :, :) = squeeze(nanmean(eegd(:,:, ind1),3));
                            outgoingEEG_Attnsplit(id, 2, :, :) = squeeze(nanmean(eegd(:,:, ind2),3));
                            outgoingEEG_Attnsplit(id, 3, :, :) = squeeze(nanmean(eegd(:,:, ind3),3));
                            outgoingEEG_Attnsplit(id, 4, :, :) = squeeze(nanmean(eegd(:,:, ind4),3));
                            
                            outgoingBEH_Attnsplit(id, :)= mBar;
                        case 2
                            outgoingEEG_XAXISsplit(id, 1, :, :) = squeeze(nanmean(eegd(:,:, ind1),3));
                            outgoingEEG_XAXISsplit(id, 2, :, :) = squeeze(nanmean(eegd(:,:, ind2),3));
                            outgoingEEG_XAXISsplit(id, 3, :, :) = squeeze(nanmean(eegd(:,:, ind3),3));
                            outgoingEEG_XAXISsplit(id, 4, :, :) = squeeze(nanmean(eegd(:,:, ind4),3));
                            
                            outgoingBEH_XAXISsplit(id, :)= mBar;
                            
                    end
                    
                end
            end
        end
        
        %save at ppant level.
        cd(eegdir);
        if uselowpfilt==1            
                cd([ xlabis ' VAN data lowpfilt quartilesplit'])            
        else
                cd([ xlabis ' VAN data quartilesplit'])
            
        end
            
                save([allppants(ippant).name],'outgoingEEG_XAXISsplit', 'outgoingEEG_Attnsplit', ...
                    'outgoingBEH_Attnsplit', 'outgoingBEH_XAXISsplit', 'time_new');
            
        end
    
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Concatenate across participants
        
        %% preallocate EEG
        
        nfiles = dir([pwd filesep xlabis '*participant' '*.mat']);
        
        
        [GFX_HIT_EEGd_Attn,GFX_MISS_EEGd_Attn,GFX_CR_EEGd_Attn,GFX_FA_EEGd_Attn,...
            GFX_HIT_EEGd_XAXIS,GFX_MISS_EEGd_XAXIS,GFX_CR_EEGd_XAXIS,GFX_FA_EEGd_XAXIS] = deal(zeros(length(nfiles),length(Qps)+1, nchans,size(FA_EEGd,2)));
        
        for ippant = 1:length(nfiles)
            %change into correct directory:
            
            cd(eegdir)
            if uselowpfilt==1                
                cd([ xlabis ' VAN data lowpfilt quartilesplit'])                
            else                
                cd([ xlabis ' VAN data quartilesplit'])               
            end
            
            %load participant data:
            load(nfiles(ippant).name); % loads all EEG and also loads time_new (epochtime axis)
            
            for itype=1:4
                
                
                %store across participants.
                switch itype
                    case 1
                        GFX_HIT_EEGd_Attn(ippant, :,:,:)= squeeze(outgoingEEG_Attnsplit(1,:,:,:));
                        GFX_HIT_EEGd_XAXIS(ippant, :,:,:)= squeeze(outgoingEEG_XAXISsplit(1,:,:,:)); % first dimension is the SDt class (HIT etc.)
                    case 2
                        GFX_MISS_EEGd_Attn(ippant, :,:,:)=squeeze(outgoingEEG_Attnsplit(2,:,:,:));
                        GFX_MISS_EEGd_XAXIS(ippant, :,:,:)=squeeze(outgoingEEG_XAXISsplit(2,:,:,:));
                    case 3
                        GFX_FA_EEGd_Attn(ippant, :,:,:)=squeeze(outgoingEEG_Attnsplit(3,:,:,:));
                        GFX_FA_EEGd_XAXIS(ippant, :,:,:)=squeeze(outgoingEEG_XAXISsplit(3,:,:,:));
                    case 4
                        GFX_CR_EEGd_Attn(ippant, :,:,:)=squeeze(outgoingEEG_Attnsplit(4,:,:,:));
                        GFX_CR_EEGd_XAXIS(ippant, :,:,:)=squeeze(outgoingEEG_XAXISsplit(4,:,:,:));
                end
                
                
                
            end % itype
            
            
        end % ippants
        
        
        %save across all:
        save(['GFX_' num2str(xlabis) '_Quantilesplit_ERPs'], ...
            'GFX_HIT_EEGd_Attn', 'GFX_HIT_EEGd_XAXIS',...
            'GFX_MISS_EEGd_Attn','GFX_MISS_EEGd_XAXIS',...
            'GFX_FA_EEGd_Attn', 'GFX_FA_EEGd_XAXIS',...
            'GFX_CR_EEGd_Attn','GFX_CR_EEGd_XAXIS', 'time_new');
        
    end % data type
end % function end

