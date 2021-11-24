% First script to convert tsv files into a workable matlab format (table).
%
% This script imports tsv files from your desktop, and saves in new folder for matlab
% analyses

% Before beginning, set your home directory (homedir), as the folder
% location you have stored/downloaded the relevant tsv files at [see below]
 

%% %%% %% set your home directory first! 

%this determines which user is running the script:
% Query type depends on Mac vs Windows log.
mydirs;
%% %% clear to begin 
clearvars -except homedir ;
close all; clc;

%% where are we?
cd(homedir); 
cd('Exp 2 and 3 tsv files')
datadir=pwd;
%% collect experiment directories here
if ~exist('nexpFiles', 'var')
    nexpFiles = dir([pwd filesep 'Exp *']);
end
for ifol = 1:length(nexpFiles) % for each experiment.
   
    cd(datadir)
    cd(nexpFiles(ifol).name)
    
    %how many participants?
    ndatafiles = dir([pwd filesep 'pasrep*']);
     
%     for each participant
      for   ippant = 1:length(ndatafiles)
        
          
    cd(datadir)
    cd(nexpFiles(ifol).name)
    
    %Dlmread can't import files. Dlmread can't handle non-numeric
    %entries.
    % There are 33 columns, and to use textscan, we need to specify the
    % filetype in each. Though first row (headers) can be read as strings
    %%
    ncol=33;
    %assign file as readable
    ftmp= fopen(ndatafiles(ippant).name);
    % use string type '%s' for headers (first row only),
     disp(['reading in data for ' ndatafiles(ippant).name ])
    header=textscan(ftmp,[repmat('%s', 1, ncol)],1,'delimiter','\t');
    % note that one header name is incompatible with matlab (was 'Q1Thr.SD'):
    header{33} = {'Q1Thr_SD'};
    
    %now note the specific filetype order, mixture of floating points and strings:
    col_formatspec = ['%s%s',repmat('%f',1,11), '%s%f%s',repmat('%f', 1,17)];
    %check length (should be 33 specifiers);
  
    %now we can read the data,
    p_data=textscan(ftmp,col_formatspec,'delimiter','\t');
    
    fclose(ftmp);
    
%% each column now contains cell arrays. 
% easiest to save as table.
    p_table=table();

    for ivar=1:ncol
        
        %first convert cell to matrix, then add as new column to table        
        try dataIN=cell2mat(p_data{ivar}); % if string col.
        catch
            dataIN=p_data{ivar}; %in case not a string.
        end
   
        %add to table
    p_table = [p_table,table(dataIN)];
   
    p_table.Properties.VariableNames{1,ivar} = cell2mat(header{ivar});
    end

    %save based on actual number in the EEG.
      searchf= ['pasrep' num2str(ifol+1) '_%d'];
        str=ndatafiles(ippant).name;
        ppantnum = sscanf(str, searchf);
    
        disp([' saving EEG for ' str ' as participant number ' num2str(ppantnum)]);
    
    if ifol==1
        savename= ['Conf_Attn_participant' num2str(ppantnum)];
    else
        savename= ['PAS_Attn_participant' num2str(ppantnum)];
    end
    %%
        cd(datadir)
        cd ../
        % change into output folder, create if needed.
        try cd('Exp 2 and 3 mat files')
        catch
            mkdir('Exp 2 and 3 mat files')
        end
        %%
    save(savename, 'p_table');
      end
    
end