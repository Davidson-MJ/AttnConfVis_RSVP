%plotExamplePretrialFeatures 
dbstop if error
plotPOWER =1;
plotPHASE =0;
plotSW    =0;
%%

RespColours = brewCOLOURS;

idatatype=1;
ippant = 1;
%%
xlabsare = {'Conf', 'PAS'};

%load channel data:
elocs=getelocs(2);
ParietoOccChan = 24:32; % as per MacDonald et al.
%%
xlabis = xlabsare{idatatype};

%% %set up directories
% connect to behavioural data:
cd(homedir)
cd('Exp 2 and 3 mat files')
behdatadir = pwd;

cd(homedir);
cd('Exp 2 and 3 processed EEG');
eegdir = pwd;


%collect participant index of interest
cd(eegdir)
nfiles = dir([pwd filesep xlabis '*.mat']);
%% load EEG
cd(eegdir)

if plotPOWER==1 || plotPHASE==1
load(nfiles(ippant).name, 'pEEG_stim_detr_dsamp');
dataEEG = pEEG_stim_detr_dsamp;
%% we also want the alpha power and complex phase angles per trial
cd([xlabis ' alpha data']);
nfiles2 = dir([pwd filesep xlabis '_Attn' '*.mat']);
load(nfiles2(ippant).name, 'AlphaMean', 'complexEEG');
Alpha_av = squeeze(mean(AlphaMean(28:32,:),1));
disp(['PAIRING EEG ' nfiles(ippant).name ' w BEH ' nfiles2(ippant).name ]);
elseif plotSW==1
    %need to use the less preprocessed ver.
    load(nfiles(ippant).name, 'pEEG_stim_prepdICA');
    dataEEG = pEEG_stim_detr_dsamp;
end 

 time_secs = dataEEG.time{1};
 
 
 %% now extract 2 example trials based on type.
 if plotPOWER==1
     ylimlab = {'Low', 'High'};
     pcol= 'k';
     %max and minimum?
     id1 =find(Alpha_av ==min(Alpha_av));
     id2 =find(Alpha_av ==max(Alpha_av));
     
     % quantiles.
qs= quantile(Alpha_av, 3);
% find nearest values to the quantiles.
val=[]; idx=[]
for iq= 1:length(qs)
[val(iq),idx(iq)]=min(abs(Alpha_av-qs(iq)));
end

showtrials = [id1, idx, id2]; % quantiles and  max alpha
     %% show result:
     dforfilt =  cat(2, dataEEG.trial{[showtrials]});
%feature sizes
nchans = length(dataEEG.label);
nsamps = length(dataEEG.time{1});
ntrials = length(dataEEG.trial);

     %high pass filt:
% filtd = eegfilt(dforfilt, 250, 0.5, 0);
%reshape and plot:

tmp= reshape(dforfilt, [nchans, nsamps,length(showtrials)]);
%take mean over chans:
plotDATA= squeeze(nanmean(tmp(ParietoOccChan,:,:),1));
 end
 
 if plotPHASE==1
%% find all phase angles at 10 Hz.
  usechan=29; % POz
   ylimlab = {'pi/2', 'pi'};
   pcol= 'b';
  %restrict range to trials with some alpha amplitude.
  cats = splitTrialsintoBins(AlphaMean, 5); % take top quarter.
  restrange = cats{6};
  tmp_Ad = squeeze(complexEEG(usechan,1,restrange));            

  % collect phase angles (in radians), of all trials:
            trial_angles = angle(tmp_Ad);
            
            %unwrap
            trial_angles = unwrap(trial_angles);

            %trial angles, to degrees
            degAng = rad2deg(trial_angles);
            
            %max and minimum?
     id1 =dsearchn(degAng,90);
     id2 =dsearchn(degAng, 180);
      % show result:
     dforfilt =  cat(2, dataEEG.trial{[restrange]});
%feature sizes
nchans = length(dataEEG.label);
nsamps = length(dataEEG.time{1});
ntrials = length(restrange);

tmp= reshape(dforfilt, [nchans, nsamps,length(restrange)]);
     
%high pass filt:
% filtd = eegfilt(dforfilt, 250, 0.5, 0);
%reshape and plot:

%take mean over chans:
plotDATA= squeeze(tmp(usechan,:,[id1,id2]));
 end
 
 
 
%% NOW plot
figure(1); clf; 
set(gcf, 'units', 'normalized', 'position', [0, .6, .27, .27])
scaleY = linspace(1, 250, length(showtrials));
alphaCols = flipud(squeeze(RespColours(4,:,:)));
ctmp = cbrewer('seq', 'Purples', 6);
for itrial = 1:size(plotDATA,2)

    
    plot(time_secs, plotDATA(:,itrial)+scaleY(itrial), 'linew', 2, 'col', ctmp(itrial+1,:));
%%
hold on

end
%
set(gca, 'Ytick', scaleY, 'YTickLabel',[], 'XTickLabel', [], 'fontsize', 20)
xlim([-.5 3])
% plot([0 0], ylim, [ 'k' ':'], 'linew', 1)
% plot([1 1], ylim, [ 'k' ':'], 'linew', 1)
% xlabel('Time from trial start [s]')

box off; 
axis off
% xlim([-.3 1])

c=colorbar('Ticks', [0,1], 'Ticklabels', [], 'location', 'WestOutside');
yl=ylabel(c, 'Alpha power', 'fontsize', 25);
% set(yl, 'Rotation', -90, 'VerticalAlignment', 'bottom');
ctmp = cbrewer('seq', 'Purples', 6);
% colormap(flipud(squeeze(RespColours(4,:,:))));
colormap(ctmp(2:6,:));
%%
figure(2); clf;
set(gcf, 'units', 'normalized', 'position', [0, .6, .27, .27])
c=colorbar('Ticks', [0,1], 'Ticklabels', [], 'location', 'WestOutside');
hold on;
axis off
ylabel(c, 'Attention', 'fontsize', 25)
colormap(flipud(squeeze(RespColours(3,:,:))));
set(gcf, 'color', 'w')

%%

