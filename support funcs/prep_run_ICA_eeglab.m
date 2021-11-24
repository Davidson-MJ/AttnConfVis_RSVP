%called from Data_explore_B_EEG  
% function    prep_run_ICA_eeglab(homedir, ~)

 %following order:


% Imports minimally preprocessed, runs ICA for eyeblink detection.


dbstop if error
cd(homedir)
cd('Exp 2 and 3 processed EEG');

nfiles = dir([pwd filesep '*.mat']);
getelocs
%%
for ifile = 12%:length(nfiles)

    cd(homedir)
    cd('Exp 2 and 3 processed EEG');
    
    % load each file in turn.
    
    load(nfiles(ifile).name, 'pEEG_stim_detr_dsamp');

    %% convert to EEGlab structure, then run ICA.
    %initialize
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%%
data = pEEG_stim_detr_dsamp;

for i=1:size(data.trial,2)
  EEG.data(:,:,i) = single(data.trial{i});
end
%%
% load chanlocs.mat
getelocs;

EEG.setname    = data.cfg.dataset;
EEG.filename   = '';
EEG.filepath   = '';
EEG.subject    = '';
EEG.group      = '';
EEG.condition  = '';
EEG.session    = [];
EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};
EEG.ref        = []; %'common';
EEG.event      = [];
EEG.epoch      = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.saved      = 'no';
EEG.chanlocs = []; 
EEG = eeg_checkset( EEG );    

%match data to channels
EEG = pop_select(EEG, 'channel', [1:32]);
% %%
EEG.chanlocs = elocs32;
EEG = eeg_checkset( EEG );    

eeglab redraw

%%
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%% reduce to correct chan number.

EEG = eeg_checkset( EEG );    

%% run ICA.
%%
EEG = pop_runica(EEG, 'extended',0,'interupt','on');
EEG = eeg_checkset( EEG );    
%%
eeglab redraw
%%
pop_selectcomps(EEG, 1:32);

%% perform rejection, and save using the below:
EEG = eeg_checkset( EEG );    



%%

pEEG_stim_prepdICA= eeglab2fieldtrip(EEG, 'preprocessing');

% append 
save(nfiles(ifile).name, 'pEEG_stim_prepdICA', '-append');


end