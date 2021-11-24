
function PFX_basedonCAT(homedir)
%plots  EEG on group level, showing the results based on categories: Hit,
%Miss, FA, and CR. also sub options for comparing these.
%%
dbstop if error
figdir = '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures';
% % % type of data to plot:
getelocs

%which type of data to plot? adjusts xlimits, and print name when outgoing.
useTargetLockedorwholetrial = 2;% 1 or 2


for dset = 1:2
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


    
    % orient and load group data
cd(homedir);
cd('Exp 2 and 3 processed EEG');
eegdir = pwd;


if useTargetLockedorwholetrial == 1
cd([ xlabis ' VAN data lowpfilt'])
typeis  = 'target-locked ERPs';
usexlim=[-.1 1];
else
    cd([ xlabis ' whole epoch data lowpfilt'])
    typeis  = 'whole epoch ERPs';
    usexlim=[-.5 3];
end

eegdir = pwd;
nfiles = dir([pwd filesep '*' 'participant*']);
%%
for ippant = 1:length(nfiles)

    cd(eegdir);
    
    
    load(nfiles(ippant).name);
    
    if useTargetLockedorwholetrial == 1
        usetime = time_new;
    else
        usetime = time_secs;
    end
    
% 
%%
%plot at which channel?
plotchannel = [16,21,26,29,31];
%set colours:
[cmaps] = cbrewer('qual', 'Set1', 5);
figure(1); clf; set(gcf, 'units', 'normalized', 'position', [0 0 .4 .6])

for idata = 1:4
    switch idata
        case 1
            datatoplot = HIT_EEGd;
        case 2
            datatoplot = MISS_EEGd;
        case 3
            datatoplot = FA_EEGd;
        case 4
            datatoplot = CR_EEGd;
        
    end

    col = cmaps(idata,:);

% plot(time_new, squeeze(mean(datatoplot(:,plotchannel,:),1)), 'color', cmaps(idata,:), 'linew', 3);
dis = squeeze(mean(datatoplot(plotchannel,:,:),1))';
stE= CousineauSEM(dis);
%%
subplot(211);
try
    sh=shadedErrorBar(usetime, squeeze(mean(dis,1)), stE, [], 1);
    sh.mainLine.Color = col ;
    sh.mainLine.LineWidth = 2;
    sh.patch.FaceColor= col ;
    sh.edge(1).Color= col ;
    sh.edge(2).Color= col ;
    hold on
    title([elocs32(plotchannel).labels])
    shleg(idata) = sh.mainLine;
catch %in case only 1 trial!
    shleg(idata)= plot(usetime, dis', 'color', col);
end
end
%%
lg=legend([shleg], {'HIT', 'MISS', 'FA', 'CR'}, 'autoupdate', 'off');
set(gca, 'fontsize', 25);
xlim(usexlim)
ylim([-10 16])
set(gcf, 'color' ,'w')
hold on; plot([0 0], ylim, ['k-'])
hold on; plot(xlim, [0 0 ], ['k-'])
ylabel('uV');
xlabel('Time from target present [s]')
%%
subplot(212)
datatoplot = squeeze(mean(HIT_EEGd,3)) -squeeze(mean(MISS_EEGd,3));

dis = squeeze(mean(datatoplot(plotchannel,:),1));
plot(usetime, dis, 'color', 'k', 'linew', 3)


% sh=shadedErrorBar(time_new, (squeeze(mean(dis,1))), stE, [], 1);
% sh.mainLine.Color = col ;
% sh.mainLine.LineWidth = 2;
% sh.patch.FaceColor= col ;
% sh.edge(1).Color= col ;
% sh.edge(2).Color= col ;
hold on
% plot(time_new, squeeze(mean(datatoplot(:,plotchannel,:),1)), 'color', 'k');

set(gca, 'fontsize', 25);
ylabel('\Delta uV');
xlabel('Time from target present [s]');
title({['ppant ' num2str(ippant) ];[elocs32(plotchannel).labels]})
xlim(usexlim)
ylim([-5 10])
pval = [];
hold on; plot([0 0], ylim, ['k-'])
hold on; plot(xlim, [0 0 ], ['k-']);
%%
legend('HIT-MISS')


cd(figdir); cd('Target locked ERPs')
cd('PFX_splitbycategory');
print('-dpng', ['participant count ' num2str(ippant) ', ' typeis ' ' xlabis]);


end %ippant

%
end
