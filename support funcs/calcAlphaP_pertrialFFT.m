% function calcAlphaP_pertrialFFT(homedir,~, ~)
% function to calculate the pre RSVP stream, alpha-power, per trial.
% This can then be averaged across different conditions, or used as an
% index of attention / subjective ratings.


% This script uses stationary FFT, to avoid contamination of stimulus
% evoked responses.


%
dbstop if error
%
% called from Data_explore_EEG
% MDavidson mjd070 dot gmail dot com
%%


% for
for idatatype=1%:2
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
% load(nfiles(ippant).name, 'pEEG_stim_detr_dsamp_ICACOMPSr');%, 'rej_trials');
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


dataEEG = pEEG_stim_detr_dsamp; %! use correct data!
%  dataEEG = pEEG_stim_prepd;
%  dataEEG=pEEG_stim_detr_dsamp_ICACOMPSr;

 time_secs = dataEEG.time{1};
%
%%

%reshape data to speed up filtering:
dforfilt =  cat(2, dataEEG.trial{:});


%unpack: Channels x samples x trials

nchans = length(dataEEG.label);
nsamps = length(dataEEG.time{1});
ntrials = length(dataEEG.trial);

if ntrials~= size(trial_table,1)
    error('check code')
end
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
% window(2,:) = [timewin_marks(3) , timewin_marks(4)]; %RSVP 
% window(3,:) = [timewin_marks(2) , timewin_marks(3)]; %short prestimulus (200ms preRSVP)
% %not that this will be updated:
% window(4,:) = [timewin_marks(2) , timewin_marks(3)]; %short window.(200ms)

%%

% Now take hilbert envelope of  alpha filtered data

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



% nfft(iwin)=max(2^(nextpow2(N)+2),N);        %2 includes zero padding.
nfft(iwin)=max(2^(nextpow2(N)),N);        %2 includes zero padding.




f=getfgrid(params.Fs,nfft(iwin),params.fpass);

fvecs(iwin).f=f;

tapers=dpsschk(params.tapers,N,params.Fs); % check tapers
tapersvec(iwin).t = tapers;
end

complexEEG = zeros(32,size(EEGdatat,3), length(f)); % channels,  trials, freqs.
%%
%now feed these parameters in where needed, below:
for itrial = 1:size(EEGdatat,3)
    
  
    for iwin=1:size(window,1)
        
        %all channels, this trial and specific window:
        winEEG = squeeze(EEGdatat(:,window(iwin,1):window(iwin,2),itrial));
        
        %take FFT of this trial, operates on columns
        
        J=mtfftc(winEEG',tapersvec(iwin).t,nfft(iwin),params.Fs); % multiplies by hann window.
        
       

complexEEG(:,itrial,:) = squeeze(J(1:length(f),1,:))';

 
    end
end

%%
%sanity check, power and ITPC should be greatest at occipital channels.
% geteloc;
% Ntrials = size(complexEEG,3);
% clf
% for ip=1:4
    
%     subplot(2,2,ip);
% TP1 = squeeze(complexEEG(1:32,ip,:)); % within RSVP stream.
%   ITPC=TP1./abs(TP1); %divide by amp to make unit length
%   ITPC=(sum(ITPC,2)); % sum angles in complex plane,
%   ITPC=squeeze(abs(ITPC)./Ntrials); %norm. average of these
% 
%   topoplot(ITPC, elocs32); colorbar; caxis([0 .5]); title(titlesare{ip});
% end

%% save data per ppant.

ampEEG = abs(complexEEG);

% plot the mean across channels.
% plot(f, squeeze(nanmean(ampEEG,2)));
% shg
%%

%8-12 hz is: 
keepf=dsearchn(f', [8 12]');

%find peak from Oz.
% tmpa = squeeze(mean(ampEEG(31,:,:),2));
% pkat = max(tmpa(keepf(1):keepf(2)));
% fmax = dsearchn(tmpa, pkat'); 
%ind adj range (adjust around individual peak freq)
% keepfP = dsearchn(f', [f(fmax)-2 f(fmax)+2]');





%%
%take relevant average
% AlphaMean_FFT= squeeze(nanmean(ampEEG(:,:,keepfP(1):keepfP(2)),3));
AlphaMean_FFT= squeeze(nanmean(ampEEG(:,:,keepf(1):keepf(2)),3));
ampEEG_faxis= f;
%%
% time_new = newtimevec;
cd(eegdir)
cd([ xlabis ' alpha data'])
save(allf(1).name, 'AlphaMean_FFT', 'trial_table','ampEEG','ampEEG_faxis', '-append');
disp(['saved alpha (fft) results for ' allf(1).name])
%%
end
% %% 
% % clf

end

% end



