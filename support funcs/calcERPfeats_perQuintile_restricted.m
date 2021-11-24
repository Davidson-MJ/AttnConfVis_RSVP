function  calcERPfeats_perQuintile_restricted(homedir)
% uses alpha quintiles, with subselection of behavioural outcomes.
%for example, quintiles for HITS only.
% need to adjust each epoch to use target locked responses.

% Focusing on the target-locked ERPs.
dbstop if error
%
% called from Data_explore_EEG
% MDavidson mjd070 dot gmail dot com
%%

%which experiment to analyze? Conf or PAS?
showPeakssaved=0; %this slows down the pipeline, but enables to debug whether correct peaks are being saved.
useChans = 29:32; % parieto-occipital subset (for alpha)


%new addition.
MeasureProminence =0; % 0=baseline to peak, 1=measure trough to peak

ParietoOccChan = 29:32; % 24:32
for idatatype=1:2
    switch idatatype
        case 1
            xlabis = 'Conf';
        case 2
            xlabis = 'PAS';
    end
    
    
    dbstop if error
    
    %load channel data:
    elocs=getelocs(2); % 2nd study
    
    
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
        
        
        %% save data per ppant, with ALPHA, to make it easy.
        
        cd(eegdir)
        cd([ xlabis ' alpha data'])
        searchf= [xlabis '_Attention_participant%d'];
        str=nfiles(ippant).name;
        ppantnum = sscanf(str, searchf);
        
        allf = dir([pwd filesep xlabis '*pant' num2str(ppantnum) '.mat']);
        %
        load(allf(1).name,  'AlphaMean_FFT', 'trial_table'); % load alpha per trial.
        AlphaMean= squeeze(nanmean(AlphaMean_FFT(ParietoOccChan,:),1));
        
        
        %% first correct for pre Stream baseline.
        
        % note that stream onset is at the 0.7s mark.
        onset = 0.7; % ms
        
        cfg=[];
        cfg.demean        = 'yes';
        cfg.baselinewindow = [onset-.15 onset]; % in seconds
        cfg.hpfilter = 'yes';
        cfg.hpfreq= 1 ; % particularly high-cut-off, as we are interested in transient P1-N1 (as per Rajagovindan & Ding, Jneuro).
        cfg.hpfiltwintype = 'hann';
        cfg.lpfilter = 'yes';
        cfg.lpfreq= 25 ; % Hz
        cfg.lpfiltwintype = 'hann';
        
        br_data = ft_preprocessing(cfg, pEEG_stim_detr_dsamp);
        nchans = 32;
        
        %% sanity check, plot the mean ERP, per ppant, POchans
        ParietoOccChan = 24:32; % as per MacDonald et al.
        POchans = ft_channelselection(ParietoOccChan, br_data);
        %
        cfg=[];
        cfg.channel = POchans;
        avERP = ft_timelockanalysis(cfg, br_data);
        % needs avERP as input >
        timeax=avERP.time;
        %         showP1N1_SSVEPintrialERP;
        
        % continue with analysis.
        % will perform the same, after first averaging within alpha quintiles.
        
        %but by selecting quintiles within certain outcomes (H,M, FA,CR).
        [P1_bySDTandAlpha, P2_bySDTandAlpha,f1_bySDTandAlpha, Twof1_bySDTandAlpha]=  deal(nan(nchans, 4, 5));
        
        %
        for isub = 1:2%4 % H,M,CR,FA
            restrange= trial_table.Outcome ==isub;
            %now take quintiles of this range:
            
            tmp_Ad = AlphaMean(restrange);
            
            %also restrict EEG to this range
            cfg=[];
            cfg.trials= find(restrange);
            cfg.keeptrials=1;
            br_datatmp  = ft_timelockanalysis(cfg, br_data);
            
            
            %>Continue by taking zscore of remaining relant trials, for
            %alpha amplitudes.
            tmp_Ad_z = zscore(tmp_Ad);
            %% now split alpha into quintiles.
            
            RespCats = splitTrialsintoBins(tmp_Ad_z,4);
            
            
            
            %% only continue if there are enough trials,
            if length(tmp_Ad_z)>5 % atleast one trials each quintile
                %now calculate for each quintile.
                for iquin = 1:5
                    trialsvector = RespCats{iquin};
                    
                    %% take average ERP over these channels.
                    %Keep trials separate.
                    cfg=[];
                    cfg.trials= trialsvector;
                    cfg.keeptrials=1;
                    %sub select from the available trials
                    try    allD = ft_timelockanalysis(cfg, br_datatmp);
                    catch
                        checkcode
                    end
                    %%
                    
                    % now take P1,N1 per channel (POchans).
                    for ichan = 29:32 % don't waste time with frontal chans for vis response.
                        
                        % dive into a similar script as above...
                        
                        
                        tmpERP = squeeze(nanmean(allD.trial(:,ichan,:),1));
                        % find zero crossings:
                        % find zero crossings:
                        if MeasureProminence~=1
                            
                            measureBasetoPeak;
                            
                            %rename the output of the above script, for
                            %this current job.
                            P1_bySDTandAlpha(ichan,isub,iquin) = P1_byalpha(ichan,iquin);
                            P2_bySDTandAlpha(ichan,isub,iquin) = P2_byalpha(ichan,iquin);
                        elseif MeasureProminence==1
                            
                            %%
                            [pks,locs,~,prom]=findpeaks(tmpERP');
                            [troughs,troughloc]= findpeaks(tmpERP'.*-1);
                            %we want the prominence, of first two large peaks in
                            %% windows:
                            figure(10); clf
                            plot(timeax, tmpERP); hold on;
                            plot([.7 .7], ylim, ['k:'])
                            xlim([.6 1.5]);
                            % plot first 5 pks (after onset).
                            onsetat = dsearchn(timeax', .7);
                            
                            avail = find(locs>onsetat);
                            %add troughs
                            availtr = find(troughloc>onsetat);
                            %add pk and prom markers
                            for ip=1:3
                                %                         indextmp = avail(ip);
                                %                         pktmp = pks(indextmp); % height
                                %                         loctmp = locs(indextmp); %latency
                                %                         promtmp = prom(indextmp); %prominence - the key.
                                %                         hold on
                                %                         plot(timeax(loctmp) , pktmp, '*')
                                %                         plot([timeax(loctmp) timeax(loctmp)], [pktmp-promtmp, pktmp], '-')
                                %
                                %add trough
                                indextmp = availtr(ip);
                                loctmp = troughloc(indextmp); %latency
                                plot(timeax(loctmp), tmpERP(loctmp), 'bo')
                            end
                            
                            %%
                            %note that we want to use the same critical windows as
                            %above
                            %find positive peaks in different windows:
                            winds(1,:)=[0.08, 0.160]; % initial search windows (90 - 160 ms) (after onset t 0.7)
                            winds(2,:)=[0.161, 0.260]; % initial search windows (90 - 160 ms)
                            
                            %retain max peak from within these windows.
                            for iwind = 1:2
                                searcht = onset+winds(iwind,:); % find peaks in this window.
                                searchtax = dsearchn(timeax', searcht');
                                avail = find(locs>searchtax(1) & locs<searchtax(2));
                                
                                
                                P1tmp = pks(avail); % magnitude of peaks
                                P1= max(P1tmp); %retain larger
                                
                                
                                if ~isempty(P1)
                                    
                                    loct= find(pks==P1); % when was max peak?
                                    %in sec
                                    pk_s = timeax(locs(loct));
                                    tr_s = timeax(troughloc);
                                    
                                    %adjust height, by subtracting previous trough.
                                    prevt = find(tr_s<pk_s,1, 'last');
                                    
                                    %adjusted height
                                    troughM = troughs(prevt)*-1; % flip again (since we used peak finder)
                                    
                                    P1 = P1-troughM;
                                    
                                    if iwind==1
                                        P1_bySDTandAlpha(ichan,isub,iquin) = P1;
                                        P1loct = locs(loct);
                                    else
                                        P2_bySDTandAlpha(ichan,isub,iquin) = P1;
                                        P2loct = locs(loct);
                                    end
                                else  % no peaks found.
                                    if iwind==1
                                        P1_bySDTandAlpha(ichan,isub,iquin) = nan;
                                        P1loct = 301; % marks at onset if missed peak.
                                    else
                                        P2_bySDTandAlpha(ichan,isub,iquin) = nan;
                                        P1loct = 301;
                                    end
                                end
                            end
                            %                  %show the stored peak(s).
                            if showPeakssaved
                                %%
                                P1t=P1_bySDTandAlpha(ichan,iquin);
                                P2t=P2_bySDTandAlpha(ichan,iquin);
                                
                                figure(10); clf
                                plot(timeax, tmpERP); hold on;
                                %addsearch window
                                %     Xs =[P1win(1), P1win(1), P1win(2), P1win(2)];
                                %     Ys = [-4 4 4 -4];
                                %     ph=patch(Xs, Ys, 'b');
                                %     Xs =[P1win(1), P1win(1), P1win(2), P1win(2)];
                                %
                                %     ph=patch(Xs, Ys, 'b');
                                %     ph.FaceAlpha = .1;
                                plot([onset onset], [ylim ], 'k:')
                                if ~MeasureProminence
                                    % add p1 peak stored.
                                    plot([timeax(P1loct) timeax(P1loct)], [0 tmpERP(P1loct)], 'b-', 'linew', 3)
                                    % add p2 peak stored.
                                    plot([timeax(P2loct) timeax(P2loct)], [0 tmpERP(P2loct)], 'b-','linew', 3)
                                else % show prominence.
                                    % add p1 peak stored, and prominence.
                                    
                                    plot([timeax(P1loct) timeax(P1loct)], [tmpERP(P1loct)-P1t, tmpERP(P1loct)], 'b-', 'linew', 3)
                                    % add p2 peak stored.
                                    
                                    plot([timeax(P2loct) timeax(P2loct)], [tmpERP(P2loct)-P2t, tmpERP(P2loct)], 'b-', 'linew', 3)
                                end
                                
                                title({['ippant ' (num2str(ippant)) ', peaks found (blue)'];...
                                    ['chan ' num2str(ichan) ', quin ' num2str(iquin)]})
                                xlim([.5 1.5])
                                pause(.5);
                                
                            end
                            
                        end
                    end % ichans
                    
                    %% now also take the SSVEP amplitude.
                    %we want to project the SSVEP of each trial, to the mean SSVEP, to only
                    %look at phase-locked component (so unrelated to prestim alpha).
                    
                    %unpack per chan;
                    trialdata.data = permute(allD.trial(:,1:nchans,:), [2,3,1]); % to EEGlab order (chans, samps, trials).
                    trialdata.srate = br_data.fsample;
                    trialdata.xmin = timeax(1);
                    trialdata.xmax = timeax(end);
                    
                    % apply gabor wavelet, and take FFT, keep in complex plane.
                    data10 = eegF_Gabor(trialdata, 10, 0.5,'complex');
                    amp10 = data10.data;
                    data20 = eegF_Gabor(trialdata, 20, 0.5,'complex');
                    amp20 = data20.data;
                    
                    % get evoked amplitudes
                    [evamp10, evamp20] =deal([]);
                    
                    for c = 1:size(amp10,1) % channel
                        for s = 1:size(amp10,2) % sample
                         %first calculate mean complex value / vector
                         %magnitude of this value (normalizes per sample)
                            phase10 = mean(amp10(c,s,:)) / abs(mean(amp10(c,s,:)));
                            %divide all individual trials by this
                            %normalized, complex value , retain real
                            %component.
                            evamp10(c,s,:) = real(amp10(c,s,:)/phase10);
                            
                            phase20 = mean(amp20(c,s,:)) / abs(mean(amp20(c,s,:)));
                            evamp20(c,s,:) = real(amp20(c,s,:)/phase20);
                        end
                    end
                    
                    
                    
                    
                    %now store the mean amplitude, over the SSVEP window, for this
                    %phase-projected envelope.
                    
                    SSVEPt = dsearchn(timeax', [1.15 2.15]'); % sec
                    
                    %save mean amplitude over SSVEP window:
                    f1_bySDTandAlpha(:, isub,iquin) =  mean(squeeze(nanmean(evamp10(:, SSVEPt(1):SSVEPt(2), :),2)),2);
                    Twof1_bySDTandAlpha(:, isub, iquin) =  mean(squeeze(nanmean(evamp20(:, SSVEPt(1):SSVEPt(2), :),2)),2);
                    
                    
                    %% %sanity check: visualize result:
%                     clf
%                     colormap('jet');
%                     subplot(1,2,1); topoplot(f1_bySDTandAlpha(:,1,iquin), elocs32); colorbar;
%                     subplot(1,2,2); topoplot(Twof1_bySDTandAlpha(:,1,iquin), elocs32); colorbar;
%                     
%                     shg
%                     colormap('jet');
                end % for 5 quintiles.
            else
                
                disp(['>>>>>>> WARNING, CATEGORY ' num2str(isub) ' skipped '])
                disp(['>>>>>>> For ' allf(1).name ])
                disp(['>>>>>>> >>>>>>> >>>>>>> '])
            end % if length(quintiles) >5
            %%
        end % by SDT metric
        
        %% now reorient and save with alpha results:
        %san check.
        
        %%
        % time_new = newtimevec;
        save(allf(1).name, 'f1_bySDTandAlpha', 'Twof1_bySDTandAlpha',...
            'P1_bySDTandAlpha','P2_bySDTandAlpha' ,'-append');
        disp(['... saved ' allf(1).name  ' alpha x ERP feats (restricted by outcome)'])
        
        
    end%ppant
end % datatype
% function

%% also using this guy:
function[EEG]=eegF_Gabor(EEG,Frequency,HalfBandwidth,Output)
if nargin ~=3 && nargin ~=4, help(mfilename), return,end
if nargin<4, Output='amplitude'; end
%%
% Frequency = 5,
% HalfBandwidth = 1;
sigma=HalfBandwidth/sqrt(2*log(2)); %Amplitude
% sigma=HalfBandwidth/sqrt(log(2)); %Power
FWHMt=2*log(2)/(pi*HalfBandwidth);

disp(['GABOR-FILTER: Center Frequency: ' num2str(Frequency) ', FWHM(Frequency): +/-' num2str(HalfBandwidth) 'Hz, FWHM(Time): +/-' num2str(500*FWHMt,4) 'ms']);
%%
f=(0:(size(EEG.data,2)-1))/size(EEG.data,2)*EEG.srate; % keine neg.Frequenzen
H=single(2*exp(-0.5*((f-Frequency)./sigma).^2)); % Gabor Kernel direkt im Frequenzraum erzeugen
H=repmat(H,[size(EEG.data,1),1,size(EEG.data,3)]);
S=ifft(H .* fft(EEG.data,[],2),[],2);
switch lower(Output)
    case 'amplitude'
        EEG.data=abs(S);  %Betrag in EEG.Data zurückschreiben
    case 'complex'
        Phasenkorrektur=exp(-1i*Frequency*2*pi*(EEG.xmin:1/EEG.srate:EEG.xmax)); %Sollte so sein, dass bei t=0 Phase von cos(x)+i*sin(x)=0 ist
        S=S.*repmat(Phasenkorrektur,[size(EEG.data,1),1,size(EEG.data,3)]);
        % Anm.:
        % tatsächl frequenz>analysefrequenz: phasenzunahme --> beschleunigung
        % tatsächl frequenz<analysefrequenz: phasenabnahme --> verlangsamung
        EEG.data=S;  %Komplexen Wert in EEG.Data zurückschreiben
    otherwise
        error([mfilename ': Unknown value ''' Output ''' for ''output''.'])
end
