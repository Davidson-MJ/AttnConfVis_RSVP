function calcAlphaPHASE_pertrial(homedir,~, ~)
% function to calculate the pre RSVP stream, alpha-power, per trial.
% This can then be averaged across different conditions, or used as an
% index of attention / subjective ratings.


%Will try to follow the methods described in MacDonald et al Frontiers 2011


% Focusing on the target-locked ERPs.
dbstop if error
%
% called from Data_explore_EEG
% MDavidson mjd070 dot gmail dot com
%%

%which experiment to analyze? Conf or PAS?

% for
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
%     for
for ippant = 1 :length(nfiles)

cd(eegdir)

rej_trials =[];
% load(nfiles(ippant).name, 'pEEG_stim_prepd', 'rej_trials');
load(nfiles(ippant).name, 'pEEG_stim_detr_dsamp');

%% also load corresponding behaviour this ppant:

searchf= [xlabis '_Attention_participant%d'];
str=nfiles(ippant).name;
ppantnum = sscanf(str, searchf);
cd(behdatadir)
allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
%
cd(behdatadir)

load(allf(1).name, 'p_table', 'SDTindex'); % load behavioural data in table form.
% adjust trial info:
p2= p_table;
p2(rej_trials,:)=[];
trial_table = p2;
%%


dataEEG = pEEG_stim_detr_dsamp;
%  dataEEG = pEEG_stim_prepd;

 time_secs = dataEEG.time{1};
%
%%

%reshape data to speed up filtering:
dforfilt =  cat(2, dataEEG.trial{:});


%unpack: Channels x samples x trials

nchans = length(dataEEG.label);
nsamps = length(dataEEG.time{1});
ntrials = length(dataEEG.trial);

EEGdatat= reshape(dforfilt, [nchans, nsamps,ntrials]);
EEGdatat=EEGdatat(1:32,:,:);
%apply highpass filter at 0.5 Hz
% filtd = eegfilt(dforfilt, 250, 0.5, 0);
%%
% ParietoOccChan = 24:32; % as per MacDonald et al.
% 

% peakPerparticipant = 10; % set peak to 10 Hz.

%Bandpass filter, 4 Hz around peak per participant:
% bp_rs_filt = eegfilt(filtd, 250, peakPerparticipant-2,peakPerparticipant+2);
%%
% EEGdata= reshape(bp_rs_filt, [nchans, nsamps,ntrials]);
% EEGdata= reshape(bp_rs_filt, [nchans, newnsamps,ntrials]);
%%
%where shall we take the phase 'at'
timewin_marks = dsearchn(time_secs', [-.3 .5 .7 1.7 ]'); % note that the 'get ready' presentation was at the -0.3s mark, RSVP onset at .7

%let's check the ITPC in the prestimulus window:
window(1,:) = [timewin_marks(1) , timewin_marks(3)]; %prestimulus. (get-ready to RSVP)
window(2,:) = [timewin_marks(3) , timewin_marks(4)]; %RSVP 
window(3,:) = [timewin_marks(2) , timewin_marks(3)]; %short prestimulus (200ms preRSVP)
%not that this will be updated:
window(4,:) = [timewin_marks(2) , timewin_marks(3)]; %short window.(200ms)

%%

% Now take hilbert envelope of  alpha filtered data
complexEEG = zeros(32,3,size(EEGdatat,3)); % channels,  windows, trials
params.tapers = [1,1];
params.Fs = 250;
params.fpass = [0 40];


%% for the 3 window sizes, precompute the FFT parameters for speed:
nfft=[];
fvecs=[];
tapersvec=[];
%%
for iwin=1:size(window,1)

    N=length(window(iwin,1):window(iwin,2));
% per window size, compute tapers, nfft, and frequency vector
nfft(iwin)=max(2^(nextpow2(N)),N);        

f=getfgrid(params.Fs,nfft(iwin),params.fpass);

fvecs(iwin).f=f;

tapers=dpsschk(params.tapers,N,params.Fs); % check tapers
tapersvec(iwin).t = tapers;
end
%%
%now feed these parameters in where needed, below:
for itrial = 1:size(EEGdatat,3)
    
    % actual target time on this trial =
    tloc = trial_table.TargChip(itrial);
    % onset in seconds:
    onset_s = .7 + (tloc-1)*.1;
    % specifiy new timewindow.
    timewintmp = dsearchn(time_secs', [onset_s-.1 onset_s+.1]'); % 100 ms short window.
    
    window(4,:) = [timewintmp(1), timewintmp(2)]; %short prestimulus.
    
    for iwin=1:size(window,1)
        
        %all channels, this trial and specific window:
        winEEG = squeeze(EEGdatat(:,window(iwin,1):window(iwin,2),itrial));
        
        %take FFT of this trial, operates on columns
        
        J=mtfftc(winEEG',tapersvec(iwin).t,nfft(iwin),params.Fs);
        
        %take complex value, closest to 10 Hz across channels:
        TENat = dsearchn(fvecs(iwin).f', 10');
        complexEEG(:, iwin, itrial) = J(TENat,1,:);
        
    end
end

%%
%sanity check ITPc should be greatest at occipital channels.
% geteloc;
% Ntrials = size(complexEEG,3);
% clf
titlesare={'1s preRSVP', '1s RSVP', '200ms periRSVP', '200ms periTarget'};
% for ip=1:4
%     
%     subplot(2,2,ip);
% TP1 = squeeze(complexEEG(1:32,ip,:)); % within RSVP stream.
%   ITPC=TP1./abs(TP1); %divide by amp to make unit length
%   ITPC=(sum(ITPC,2)); % sum angles in complex plane,
%   ITPC=squeeze(abs(ITPC)./Ntrials); %norm. average of these
% 
%   topoplot(ITPC, elocs32); colorbar; caxis([0 .5]); title(titlesare{ip});
% end

%% save data per ppant.

cd(eegdir)
cd([ xlabis ' alpha data'])
%%
% time_new = newtimevec;
phasedimsare={'1s preRSVP', '1s RSVP', '200ms preRSVP', '200ms periTarget'};
save(allf(1).name, 'complexEEG','phasedimsare', '-append');
disp(['saved phase results for ' allf(1).name])
%%
end
% %% 
% % clf

end

end



