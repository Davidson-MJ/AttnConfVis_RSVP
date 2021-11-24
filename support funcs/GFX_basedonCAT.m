
function GFX_basedonCAT(homedir)
%plots  EEG on group level, showing the results based on categories: Hit,
%Miss, FA, and CR. also sub options for comparing these.
%%

% % % type of data to plot:
%which type of data to plot? adjusts xlimits, and print name when outgoing.
useTargetLockedorwholetrial = 2;% 1 or 2


figdir= '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures';
getelocs
for dset = 1%:2
    figure(dset);
    
    switch dset
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


%% %%%% 
%%%%%% 
%%%%%% %%  PLOTTING
%%%%%% 
%%%%%% 


if useTargetLockedorwholetrial == 1
% cd([ xlabis ' VAN data lowpfilt'])
cd([ xlabis ' targlocked data lowpfilt'])
typeis  = 'target-locked ERPs';
usexlim=[-.1 1];
else
    cd([ xlabis ' whole epoch data lowpfilt'])
    typeis  = 'whole epoch ERPs';
    usexlim=[-.5 3];
end
   
%% load data
load(['GFX_' xlabis '_Attn_ERPs.mat']);
%%
if useTargetLockedorwholetrial == 1
    usetime = time_new;
    usexlabel = 'Time from target present [s]';
else
    usetime = time_secs;
    usetime = usetime+.3;
    usexlabel = 'Time from trial start [s]';
end


%plot at which channel?
plotchannel=[16,21,26,29,31];
%set colours:
[cmaps] = cbrewer('qual', 'Set1', 5);
figure(1); clf; set(gcf, 'units', 'normalized', 'position', [0 0 .4 .6])
% col='k';
for idata = 1:4
    switch idata
        case 1
            datatoplot = GFX_HIT_EEGd;
        case 2
            datatoplot = GFX_MISS_EEGd;
        case 3
            datatoplot = GFX_FA_EEGd;
        case 4
            datatoplot = GFX_CR_EEGd;
        
    end

    col = cmaps(idata,:);

% plot(time_new, squeeze(mean(datatoplot(:,plotchannel,:),1)), 'color', cmaps(idata,:), 'linew', 3);
dis = squeeze(mean(datatoplot(:,plotchannel,:),2));
stE= CousineauSEM(dis);

subplot(211);
sh=shadedErrorBar(usetime, squeeze(mean(dis,1)), stE, [], 1);
sh.mainLine.Color = col ;
sh.mainLine.LineWidth = 2;
sh.patch.FaceColor= col ;
sh.edge(1).Color= col ;
sh.edge(2).Color= col ;
hold on
% title([elocs32(plotchannel).labels])
shleg(idata) = sh.mainLine;
end
lg=legend([shleg], {'Hit', 'Miss', 'FA', 'CR'}, 'autoupdate', 'off');
set(gca, 'fontsize', 20);
% xlim(usexlim)
axis tight
ylim([-3 10])
set(gcf, 'color' ,'w')
hold on; plot([0 0], ylim, ['k:'], 'linew', 1)
if useTargetLockedorwholetrial==2
hold on; plot([1 1], ylim, ['k:'], 'linew', 1)
for im = 1:10
hold on; plot([1+im/10 1+im/10], [8 10], ['k:'], 'linew', 1)
end
 plot([1 2], [8 8], ['k:'], 'linew', 1)
 plot([2.2 2.2 ], [ylim], ['k:'], 'linew', 1)
hold on; plot(xlim, [0 0 ], ['k-'])
end
ylabel('uV');
xlabel(usexlabel)
hold on; 

%%
subplot(212)

% TIMINGS = nan(1, length(dis));
% getready = dsearchn(usetime', [.3 
% datatoplot = TIMINGS;

%%
datatoplot = GFX_HIT_EEGd - GFX_MISS_EEGd;
dis = squeeze(mean(datatoplot(:, plotchannel,:),2));
stE = CousineauSEM(dis);
col='k';
sh=shadedErrorBar(usetime, (squeeze(mean(dis,1))), stE, [], 1);
sh.mainLine.Color = col ;
sh.mainLine.LineWidth = 2;
sh.patch.FaceColor= col ;
sh.edge(1).Color= col ;
sh.edge(2).Color= col ;
hold on
% plot(time_new, squeeze(mean(datatoplot(:,plotchannel,:),1)), 'color', 'k');
hold on; plot([0 0], ylim, ['k-'])
hold on; plot(xlim, [0 0 ], ['k-']);
set(gca, 'fontsize', 25);
ylabel('\Delta uV');
xlabel('Time from target present [s]');
title([elocs32(plotchannel).labels])
xlim(usexlim)
% ylim([-2 4])
axis tight
pval = [];
for ip = 1:length(usetime)
    
    [~, pval(ip)] =ttest(dis(:,ip));
end
sigs = find(pval<.05);

ytx= get(gca, 'ylim');
for iplace = 1:length(sigs)
    pt= plot([usetime(sigs(iplace))], ytx(1)+.1, 'Marker', '*', 'MarkerSize', 15, 'Color', 'k');
end
%%
xlim 'auto'
legend(sh.mainLine, 'Hit-Miss')
cd(figdir)
cd('Target locked ERPs');
cd('GFX_both experiments')
print('-dpng', ['GFX_' xlabis ' exp all ' typeis ' CATs']) 
%% plot topography over cluster:
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                
                timewClust = [clusterSTandEND(maxClust,1):clusterSTandEND(maxClust,2)];
                
                eegtime = sigs(timewClust);
                topoClust = squeeze(mean(datatoplot(:,1:32, eegtime),3)); % mean over cluster.
                topoClustm = squeeze(mean(topoClust,1)); %mean over ppants;

                
                figure(2);
clf; set(gcf, 'units', 'normalized', 'position', [0 0 .35 .45])
                topoplot(topoClustm, elocs32);
                set(gcf, 'color', 'w')
                c=colorbar;
                ylabel(c, '\Delta uV');
                caxis([ -4 4])
                colormap('magma');
                set(gca, 'fontsize', 25);
                               
                print('-dpng', ['GFX_' xlabis ' exp all, ' typeis ' HIT-MISS cluster topography'])
% %%
% figure(dset+10); 
% 
% 
%     %times for topography
%     showt=[.1,.2];    
%     topoX=dsearchn(time_new', showt');
%     
%     
%     
%     for plotts = 1:2
%         plotspot = plotts + 2*(dset-1);
%         
%         subplot(3,4,plotspot);
%         
% 
%         mtdata=squeeze(mean(datatoplot(:,1:32,:),1));
%         topoplot(mtdata(:,topoX(plotts)), elocs32);
%         
%         c=colorbar;
%         caxis([-1 1])
%         title([num2str(showt(plotts)) 'ms'])
%         set(gca, 'fontsize', 15)
%         ylabel(c, 'uV')
%     end
%     %%
% %     plotspot = [5:6,9:10] + 2*(dset-1);
%      figure(1)
%     plot(time_new, mtdata, 'k')
%     set(gca, 'ydir', 'reverse')
%     hold on;
%     plot([showt(1) showt(1)], [ -2 2], 'r', 'linew', 4)
%     plot([showt(2) showt(2)], [ -2 2], 'r', 'linew', 4)
% 
% %also plot the peaks.
% t_VANwin = dsearchn(time_new', [.15, .25]');
% t_LPwin = dsearchn(time_new', [.3, .4]');
% 
% % %% average over these windows and show topography of each effect:
% % subplot(2,2,2);
% % %VAN:
% % topoplot(mean(mtdata(:, t_VANwin(1):t_VANwin(2)),2), elocs32); 
% % colorbar
% % %LP
% % subplot(2,2,4)
% % topoplot(mean(mtdata(:, t_LPwin(1):t_LPwin(2)),2), elocs32); 
% % colorbar



end
