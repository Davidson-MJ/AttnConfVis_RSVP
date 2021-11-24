%mydirs
%% %%% %% set your home directory first! 

%this determines which user is running the script:
% Query type depends on Mac vs Windows log.

if isunix() % ID if Mac/Windows
    queryuser = 'USER'; % for MAC
else
    queryuser = 'username'; % for Windows
end

switch getenv(queryuser) % find username, set data directory.
    
    case 'mdavidson'    
        
        % My local data repository
        homedir='/Users/mdavidson/Desktop/Frontiers Project/Data';
       

    case 'matthewdavidson' % working off external HD
        homedir= '/Volumes/MattsBackup (2TB)/Frontiers Project/Data';
        
     
        %     case  'newuser' . 
% YOUR local data depository.
%         homedir = YOUR DATA DIRECTORY


end
datadir=[ homedir filesep 'Exp 2 and 3 mat files'];

cd(homedir)
cd ../Documents/Figures
figuredir=pwd;
cd(homedir);

cd(homedir);
cd('Exp 2 and 3 processed EEG');
eegdir = pwd;
%set up plot output.

xratings = {'Confidence' , 'Visibility'};
xlabsare = {'Conf', 'PAS'};
