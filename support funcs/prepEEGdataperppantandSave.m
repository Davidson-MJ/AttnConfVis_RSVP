function     prepEEGdataperppantandSave(homedir, datadir, Mycfg);
%% cycles through participant EEG structure, and performs preprocessing in
% following order:


% rejects based on artefacts (automatic)
% references to average reference
% baseline subtraction (for ERP stimulus locked analyses).


dbstop if error
cd(homedir)
cd('Exp 2 and 3 processed EEG');

nfiles = dir([pwd filesep '*.mat']);
getelocs
%%
%step through and perform preprocessing jobsL

fixbaseline_subtraction=1;

for ifile=1%:length(nfiles)% 13:nfiles are the PAS subset
    
    clearvars  pEEG*
    close all
    
    load(nfiles(ifile).name); % stimulus aligned epochs
    
    
    % continue with ICA and art rej only if not previously completed!
    
    %% perform ICA, followed by ART rej?
%     if  ~exist('pEEG_stim_detr_dsamp_ICACOMPSr', 'var')
%         
%         
%         %need to perform ICA, and remove blinks etc.
%         
%         %filter and detrend to improve stationarity of EEG data:
%         cfg=[];
%         cfg.demean = 'yes';
%         cfg.detrend = 'yes';
%         cfg.baselinewindow = [-0.1 0];
%         cfg.hpfilter = 'yes';
%         cfg.hpfiltertype = 'firls';
%         cfg.hpfreq = 1; % just for ICA, 
%         
%         %new data for ICA:
%         dforICA = ft_preprocessing(cfg, pEEG_stim_detr_dsamp);
%         
%         cfgica=[];
%         cfg.method = 'runica';
%         ICAcomps = ft_componentanalysis(cfgica, dforICA);
%         
%        
%         %% save comps quickly for next time
%         save(nfiles(ifile).name,'ICAcomps', '-append');
%         
%         %%
%     
%         %see Chaumon et al., JNeuroMethods, 2015. for criteria used to reject blink and muscle
%         %components:        
%         
%          %% %% identify artifacts via ICA
%         
%          %Place topoplots on top right of screen:
%         figure(1); clf; 
%         set(gcf, 'units', 'normalized','position', [-.5 0 .5 1])
%         
%         cfg = [];
%         cfg.component = 1:25;       % specify the component(s) that should be plotted
%         cfg.layout = 'quickcap64.mat';
%         cfg.comment   = 'no';
%         %plot
%         ft_topoplotIC(cfg, ICAcomps)
%         %%
%         % also plot time-course of components      
%         cfg.viewmode = 'component';      
%         %plot
%         ft_databrowser(cfg, ICAcomps);
%         
%         %% also plot mean ERP for each component        
%         TL =ft_timelockanalysis(cfg, ICAcomps);        
%         figure(); clf; set(gcf, 'units', 'normalized','position', [.5 0 .5 1])
%         for ICAc=1:25
%             subplot(5,5,ICAc);
%             plot(TL.time, TL.avg(ICAc,:)); hold on; plot(xlim,[ 0 0], 'color', 'r', 'linew', 1);
%             ylim([-1.5 1.5])
%             title([num2str(ICAc) ' comp ERP'])
%         end
%         
%         
%         %% Now await input, and continue with ICA rejection after visual inspection.
%         disp('------- Press [enter] to continue ---------');
%         yesKeys = KbName('return');
%         prog=0; ListenChar(2)
%         % wait for enter key to continue.
%         while ~prog
%             [a,b,keyCode] = KbCheck;
%             if any(keyCode(yesKeys))
%                 prog=1;
%                 ListenChar(0)
%             end
%         end
%         close all
%         %%
%         
%         answer= cell2mat(inputdlg('Which components do you want to reject?:'));
%         rejcomps = str2num(answer);
%         %%
%         
%     end
    if fixbaseline_subtraction==1
        
    load(nfiles(ifile).name, 'rejcomps', 'pEEG_stim_detr_dsamp', 'ICAcomps'); % stimulus aligned epochs
    
        %remove the components:
        cfg = [];
        
        cfg.component = rejcomps; 
        % to be removed component(s), note that we are using the ICAs from
        % data with improved stationarity,.
        gd_data = ft_rejectcomponent(cfg, ICAcomps, pEEG_stim_detr_dsamp);
        
        %%
        %     %next, because our data contains large ERPs (SSVEPs), as well as motor
        %     %responses, we will reject more carefully, i.e. interactively.
        %     now demean, to interpret next stage (artf rejection).
        cfg=[];
        cfg.demean = 'yes';
        cfg.detrend = 'no';
        cfg.baselinewindow = [-0.5 -0.3]; % beware, the 'trial start' triggercode is actually screen flip offset (on screen for 0.3s).
        cfg.hpfilter = 'no';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 48; % as per MacDonald et al.,
        cfg.hpfiltertype = 'firws';
        cfg.lpfiltertype = 'firws';        
        pEEG_stim_detr_dsamp_ICACOMPSr = ft_preprocessing(cfg, gd_data);
        
       %% 
        %save as we go
        
        %sanity check, 
        
        
        save(nfiles(ifile).name, 'rejcomps', 'pEEG_stim_detr_dsamp_ICACOMPSr', '-append');
        close all
    end
    
    %%
    %perform trial and channel rejection if it hasn't been done
    if ~exist('pEEG_stim_prepd', 'var')
        cfg=[];
        cfg.method='summary';
        cfg.layout = 'quickcap64.mat';
        
        %don't show EOG/ MASTOIDS in rejection view.
        chansel = ft_channelselection({'all', '-VEOG', '-HEOG', '-M2'}, pEEG_stim_detr_dsamp_ICACOMPSr.label);        
        cfg.channel = chansel;
        
        % load GUI
        [pEEG_stim_prepd,rej_trials, rem_chans] =ft_rejectvisual(cfg,pEEG_stim_detr_dsamp_ICACOMPSr);
        
        if ~isempty(rem_chans)
        %% repair missing channels:
        %prepare neighbours, from all data (pre channel rejection)                
        cfg=[]; 
        cfg.method = 'triangulation';
        cfg.layout = 'quickcap64.mat';
        [neighbs]= ft_prepare_neighbours(cfg, pEEG_stim_detr_dsamp_ICACOMPSr);
        
        %% now interpolate the missing channels.
        cfg=[];
        cfg.method = 'average';
        cfg.neighbours = neighbs;
        badchannel = ft_channelselection(rem_chans, pEEG_stim_detr_dsamp_ICACOMPSr.label);
        cfg.badchannel = badchannel;
        cfg.layout = 'quickcap64.mat';
        pEEG_stim_prepd = ft_channelrepair(cfg, pEEG_stim_prepd);
        end
        %% note that we have already rejected trials in the above, but also store index of
        % rejected trials to match with Behaviour.
        
%         %sanity check before saving, quick plot of ERP at POz
%         tmpT= cat(2, pEEG_stim_prepd.trial{:});
%         EEGd= reshape(tmpT, [32, length(pEEG_stim_prepd.time{1}), length(pEEG_stim_prepd.trial)]);
%         clf; plot(pEEG_stim_prepd.time{1}, squeeze(mean(EEGd(29, :, :),3))); shg;
        %%
        %ssave per participant
        save(nfiles(ifile).name, 'rej_trials', 'rem_chans', 'pEEG_stim_prepd', '-append');
    end
    
end



end