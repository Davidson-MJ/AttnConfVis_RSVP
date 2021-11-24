%setpaths
if isunix() % ID if Mac/Windows
    queryuser = 'USER'; % for MAC
else
    queryuser = 'username'; % for Windows
end

switch getenv(queryuser) % find username, set data directory.
    
    case 'mdavidson'
        
        % My local data repository
        homedir='/Users/mdavidson/Desktop/Frontiers Project/Data';
        datadir=[ homedir filesep 'Exp 2 and 3 cnt files'];
        
        %add function support files.
        addpath([pwd filesep 'support funcs'])
        
        %     case  'newuser' .
        % YOUR local data depository.
        %         homedir = YOUR DATA DIRECTORY
        
end
figuredir = '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures/Alpha and BEH';

% set up directories
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
