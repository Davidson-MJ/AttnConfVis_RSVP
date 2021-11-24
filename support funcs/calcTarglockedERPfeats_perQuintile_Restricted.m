function  calcTarglockedERPfeats_perQuintile_Restricted(homedir) % uses predetermine alpha quintiles.
% need to adjust each epoch to use target locked responses.

% Focusing on the target-locked ERPs.
dbstop if error
%
% called from Data_explore_EEG
% MDavidson mjd070 dot gmail dot com
%%

%which experiment to analyze? Conf or PAS?
nquintiles =5;

ParietoOccChan = 29:32; % 

for idatatype=1:2
    switch idatatype
        case 1
            xlabis = 'Conf';
        case 2
            xlabis = 'PAS';
    end
    
    
    dbstop if error
    
    %load channel data:
    getelocs
    
    
    
    cd(homedir);
    cd('Exp 2 and 3 processed EEG');
    eegdir = pwd;
    
    %% begin:
    
    %collect participant index of interest
    cd(eegdir)
    cd([xlabis ' targlocked data lowpfilt'])
    %%
    nfiles = dir([pwd filesep xlabis '*.mat']);
    
    %% load EEG
    %     for
    for ippant = 1 :length(nfiles)
        
         cd(eegdir)
         cd([xlabis ' targlocked data lowpfilt'])
        
        load(nfiles(ippant).name, 'TargLocked_EEGd', 'time_new');
        
       [nchans, nsamps, ntrials] = size(TargLocked_EEGd);
        %% save data per ppant, with ALPHA, to make it easy.
        
        cd(eegdir)
        cd([ xlabis ' alpha data'])
         searchf= [xlabis '_Attn_participant%d'];
        str=nfiles(ippant).name;
        ppantnum = sscanf(str, searchf);
        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
        %
        load(allf(1).name,  'AlphaMean_FFT', 'trial_table');
        
        qtypes= {'alpha', 'xaxis', 'attention ratings'};
        for iqtype = [1,2,3]%2,3] % how should we define the quintiles? i.e. which type of data
            switch iqtype
                case 1
                    usedforquins= squeeze(nanmean(AlphaMean_FFT(ParietoOccChan,:),1));
                
                case 2
                     usedforquins= trial_table.Resp; %XrespQuintiles;
                     
                case 3
                    usedforquins = trial_table.AttResp; %AttrespQuintiles;
            end
            
        [TL_byquin_erp]=  deal(nan(nchans, nsamps, 4,5));
        
        %%
        
        %now for each type, we also want to subselect by OUTCOMES.
         for isub = 1:2 % H,M,CR,FA
            restrange= trial_table.Outcome ==isub;
            %now take quintiles of this range:
             tmp_Ad = usedforquins(restrange);
            %restrict EEG to same subset of trials
             tmp_EEG = squeeze(TargLocked_EEGd(:,:,restrange));
             
             
            %>Continue by taking zscore of remaining relevant trials, for
            %alpha amplitudes.
            tmp_Ad_z = zscore(tmp_Ad);
            %% now split alpha into quintiles.
            
            RespCats = splitTrialsintoBins(tmp_Ad_z,nquintiles-1);
            
            if length(unique(tmp_Ad_z))>=5 % atleast one trial in each quintile
                for iquin = 1:length(RespCats)
                    trialsvector = RespCats{iquin};
                    
                    %% take average ERP over these channels.
                    TLtmp = squeeze(nanmean(tmp_EEG(:,:,trialsvector),3));
                    
                    %whole target-locked erp
                    TL_byquin_erp(:, :, isub, iquin) = TLtmp;
                                                           
                    
                    
                end
            else
                
                disp(['>>>>>>> WARNING, EXP' num2str(idatatype) ' CATEGORY : ' qtypes{iqtype} 'pp ' num2str(isub) ' skipped '])
                disp(['>>>>>>> insufficient trial range for ' allf(1).name ])
                disp(['>>>>>>> >>>>>>> >>>>>>> '])
            end % if length(quintiles) >5
            
         end
        
        
        %now save back.
        switch iqtype
            case 1
                TL_bySDTandAlpha_erp = TL_byquin_erp;
            case 2
                
                TL_bySDTandXresp_erp = TL_byquin_erp;
            case 3
                
                TL_bySDTandAttresp_erp = TL_byquin_erp;
        end
        
        end
        %%
        time_TLerp = time_new;
        
        save(allf(1).name, 'time_TLerp',...
            'TL_bySDTandAlpha_erp', ...
            'TL_bySDTandAttresp_erp',...);
            'TL_bySDTandXresp_erp','-append'); % save each TL erp, by relevant quintile split
            
        disp(['>>>>>> Saved Target locked ERP, by alpha * restricted outcomes, ppant ' allf(1).name])
        
        
        
        
    end%ppant
end % datatype
 % function
    
 