function extractEEGdataperppantandSave(homedir, datadir, Mycfg)
% Steps through the cnt files per participant, performs Epoching based on trigger codes
% - at the moment this script epochs based on the presentation of the
% target.

%
%
%
% EEG trigger codes:
% {'240', '243','244', '245', '246', '247', '248'} = trial start 'get
% ready', the last number indexes the type of trial.
% 240  = no target, 24* *=target location in sequence.

% '250' = responded absent 
% '251' = responded present
%

%%%%%% MD August 2019

        
dbstop if error
        
cd(homedir)
cd(datadir)
%%

for iExp=1:2
    
    cd(homedir); cd(datadir);
    if iExp==1        
        xlabis='Conf';   
    else        
        xlabis='PAS';        
    end
    
    cd(['Exp ' num2str(iExp+1) ' data EEG']);
    
    %how many data types/ppants this folder.
    allppants = (dir([pwd filesep 'pasrep' num2str(iExp+1) '*.cnt']));        
    nppants = length(allppants);
    
    for ippant=1:nppants

            cd(homedir); cd(datadir);
            
            cd(['Exp ' num2str(iExp+1) ' data EEG']);
    
        % Two types of epoching, start of trial based, and  response based
        % load cnt into EEGlab
        
        % careful, as matlab indexes start at 10 (then 1), so define name
        % per ppant explicitly. Also there is no "1" in exp3
        
        if iExp==1 % see above
            loadvarname=['pasrep' num2str(iExp+1) '_' num2str(ippant) '.1.cnt'];
        else
            loadvarname=['pasrep' num2str(iExp+1) '_' num2str(ippant+1) '.1.cnt'];
        end
        
        %EPOCH data based on trial start, or response:        
        for epochat = 1%:2
            
            switch epochat
                case 1          % trial start          
                    eventsare = [240, 243,244, 245, 246, 247, 248];
                    % note that 240 = target absent,
                    % 24* = target present at * (sorted later)
                    
                    EpochLength = Mycfg.Epochlength_stim;
                case 2
                    eventsare = [250 ,251];
                    %resp  absent, resp present.                    
                    EpochLength = Mycfg.Epochlength_resp;
                
                case 3 % epoch at 'GET READY', +- 
                    
                    
            end
                        
        %% configure dataset:
        cfg=[];       
        cfg.datafile = loadvarname;       
        %which events to extract?        
        cfg.trialdef.eventtype = 'trigger';
        cfg.trialdef.eventvalue = eventsare;        
        cfg.trialdef.prestim = EpochLength(1);
        cfg.trialdef.poststim = EpochLength(2);
        
        cfg = ft_definetrial(cfg);
        
        
 
        %note that cfg.trl now contains our timestamps for epoching (1st
        %and 2nd column), offset (3rd col), and original event value (4th
        %col).
        
        %% proceed with first round of preprocessing:
        
        %no LP filter so as not to temporally displace events. will demean
        %and detrend instead.
        cfg.lpfilter                    = 'no';

        cfg.demean = 'yes'; 
        cfg.detrend = 'no';
%       
        %run ! ~ 10 mins per subj.
        tic
        data_detr = ft_preprocessing(cfg);
        toc
        %%
        %downsample the data to speed up the next steps:
        cfg.resamplefs = 250;
        cfg.detrend = 'no';
        
        data_ds= ft_resampledata(cfg,data_detr);
        
        %save
        cd(homedir)
        cd('Exp 2 and 3 processed EEG');
        %save to which file?
        if iExp==1
        allf = dir([pwd filesep xlabis '*t' num2str(ippant) '.mat']);
        else
            allf = dir([pwd filesep xlabis '*t' num2str(ippant+1) '.mat']);
        end
            
        %append data.
        
        switch epochat
            case 1
                pEEG_stim_detr_dsamp = data_ds;
                save(allf(1).name,  'pEEG_stim_detr_dsamp',  'cfg', 'Mycfg');
            case 2
                
                pEEG_resp_detr_dsamp = data_ds;
                save(allf(1).name,  'pEEG_resp_detr_dsamp', 'cfg', '-append');
        end
        
        
       
        end % epoch at
    end % nppants
        
      
      
end %iEXP %nppants
    

end


%in eeglab

