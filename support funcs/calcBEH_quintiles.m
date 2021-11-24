function calcBEH_quintiles(homedir)
%

% MDavidson mjd070 dot gmail dot com
%%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nquintiles=5;
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
    
    %load data per participant
    cd(eegdir)
    cd([ xlabis ' targlocked data lowpfilt'])
    nfiles = dir([pwd filesep xlabis '*.mat']);
    %%
    
    allppants = dir([pwd filesep '*_Attn_participant*']);
    
    for ippant = 1:length(allppants)
        
        %% load behavioural data:
        cd(behdatadir)
        searchf= [xlabis '_Attn_participant%d'];
        str=nfiles(ippant).name;
        ppantnum = sscanf(str, searchf);
        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
        %
        load(allf(1).name,  'p_table');
        
        
        %%back to EEG
        cd(eegdir)
        cd([ xlabis ' targlocked data lowpfilt'])
        
        %now for each type of response, sort by attention and confidence.
        
        %         [outgoingEEG_Attnsplit, outgoingEEG_XAXISsplit ]= deal(zeros(4,4,nchans,size(MISS_EEGd,2))); % SDT x tercile, X chans x samps.
        %
        %         [outgoingBEH_Attnsplit , outgoingBEH_XAXISsplit] = deal([]); % use structur since different trial counts.
        
        
        
        %ok, sort trials based on quartile split along y-dimension (ATTENTION).
        %%
        for IAXISSPLIT = 1:2
            
            switch IAXISSPLIT
                case 1      %X
                    useDATA= p_table.Resp;
                    
                    if idatatype==2 % for visibility case
                        useDATA=useDATA(useDATA>0);
                    end
                case 2 %
                    useDATA= p_table.AttResp;
            end
            
            %take quintiles:
            %Continue by taking zscore of remaining relant trials, for
            %alpha amplitudes.
            tmp_Ad_z = zscore(useDATA);
            
            %% now split alpha into quintiles.
            
            RespCats = splitTrialsintoBins(tmp_Ad_z,nquintiles-1);
            
            
            if IAXISSPLIT==1
                XrespQuintiles= RespCats;
            else
                AttrespQuintiles=RespCats;
            end
            
            
        end
        
        disp(['saving pairing: ' nfiles(ippant).name ' with ' allf(1).name ]);
        
        save(nfiles(ippant).name, 'XrespQuintiles', 'AttrespQuintiles', 'p_table', '-append');
        
    end
    
    
end % data type
end % function end

