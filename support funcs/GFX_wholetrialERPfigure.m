
% function GFX_wholetrialERPfigure(homedir, figdir)
%plots  EEG on group level, showing the results based on  Experiment

% FOr results split by categories: Hit, Miss, FA, and CR. and 
% also sub options for comparing these. See: 'GFX_basedonCAT'
%%

% % % type of data to plot:
%which type of data to plot? adjusts xlimits, and print name when outgoing.

%plot at which channel?
  usexlim=[-.2 2.8]; 
  useylim=[-3 7]; 
  
ParietoOccChan = 29:32; % 
OccChan = 30:32; %  % for the P2 

CPchans= [11,16,21,15,17];
CustomChan=  [16,21,26,29,31]; % midline (same as Frontiers 2010).
plotchannel = CustomChan;

%set colour
RespColours= brewCOLOURS;
% create grey map.
grad= linspace(1,.9, 100);
greymap= [ grad',grad',grad'; flipud(grad'), flipud(grad'),flipud(grad')];


P1xlims= [0.9 1.2];
figure(1); clf; set(gcf, 'units', 'normalized', 'position', [0 0 1 1])
figure(2); clf; set(gcf, 'units', 'normalized', 'position', [0 0 .4 .9])
elocs= getelocs(2);

fs=30; % fontsize.
for dset = 1:2
    hold on
%     subplot(2,1,1);
    switch dset
        case 1
xlabis = 'Conf';
useCol= RespColours(1,1,:);
        case 2
xlabis = 'PAS';
useCol= RespColours(2,1,:);
    end


dbstop if error

%load channel data:
% getelocs


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
%load whole Epoch data:
    cd([ xlabis ' grandaverage whole epoch for plots'])
    typeis  = 'whole epoch ERPs';
     
    %% load data    
    load(['GFX_' xlabis ' grandERPs']);
%     GrandEEG = GFX_grandAverageERP_P1ver; % extra filtering.
    GrandEEG = GFX_grandAverageERP;
    
    % load TL data for insert panel:
        cd(eegdir)
        cd( [xlabis ' alpha data']);      
        %
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsTL_erp', 'time_TLerp');
 

% time axis
    usetime = timeax;
    usetime = usetime+.3;
    usexlabel = 'Time [s]';


% col='k';
%%
figure(1); hold on;
% shading first:
if dset==1
xpatch = [1.3 1.3 1.8 1.8];
ypatch =  [-8 10 10 -8];
colorpatch= [0 1 1 0 ]; % sets a gradient
ph= patch(xpatch, ypatch, colorpatch);
ph.FaceAlpha=.8;
ph.LineStyle='none';
colormap(greymap) % created above.

% same for P1
xpatch = [P1xlims(1) P1xlims(1) P1xlims(2) P1xlims(2)];
ph= patch(xpatch, ypatch, colorpatch);
ph.FaceAlpha=.8;
ph.LineStyle='none';
end
%%

dis = squeeze(mean(GrandEEG(:,plotchannel,:),2));
stE= CousineauSEM(dis);

% plot mean over ppants with SE
sh=shadedErrorBar(usetime, squeeze(mean(dis,1)), stE, [], 1);
sh.mainLine.Color = useCol ;
sh.mainLine.LineWidth = 4;
sh.patch.FaceColor= useCol ;
sh.edge(1).Color= useCol ;
sh.edge(2).Color= useCol ;
hold on
leg(dset)= sh.mainLine;



set(gca, 'fontsize', 20);
% xlim(usexlim)
axis tight

ylim([useylim])

set(gcf, 'color' ,'w')
hold on; 
rs=plot([1 1], [-8 10], ['k--'], 'linew', 2);
hold on; plot(xlim, [0 0], ['k-'])


% 
% % guide lines for image onset?
for im = 1:9
% hold on; plot([1+im/10 1+im/10], [18 20], ['k:'], 'linew', 1)
hold on; plot([1+im/10 1+im/10], [-8 -6], ['k:'], 'linew', 1)
end
 plot([2.2 2.2], [-8 -6], ['k:'], 'linew', 1)


ylabel('\muV');
xlabel(usexlabel)
hold on; 
box off
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 30)
%%
xlim([usexlim])
%%
figure(2);
subplot(221);
%plot P1

dis = squeeze(mean(GFX_grandAverageERP_P1ver(:,plotchannel,:),2));
stE= CousineauSEM(dis);
sh=shadedErrorBar(usetime, squeeze(mean(dis,1)), stE, [], 1);
sh.mainLine.Color = useCol ;
sh.mainLine.LineWidth = 2;
sh.patch.FaceColor= useCol ;
sh.edge(1).Color= useCol ;
sh.edge(2).Color= useCol ;

pleg(dset)= sh.mainLine;
hold on
hold on; 
rs=plot([1 1], ylim, ['k--'], 'linew', 2);
hold on; plot(xlim, [0 ,0], ['k-'], 'linew', 1)
xlim(P1xlims)
ylim([-2 5])
hold on;
ylabel('\muV')
set(gca, 'ytick', [-2 5], 'xtick', [0.9 1.2], 'fontsize', fs, 'color',[.95 .95 .95]);
xlabel('Time')
set(gcf, 'color', 'w')
% find peak at P1 and add line
searchlims= dsearchn(usetime', [1.05 1.15]');
Yd=sh.mainLine.YData;
mxY = max(Yd(searchlims(1):searchlims(2)));
mxAT = find(Yd==mxY);

pleg(1)= plot([usetime(mxAT) usetime(mxAT)], [0 mxY], 'color','r', 'linew', 3);



legend([pleg(1) rs],{'P1', 'RSVP onset'}, 'Location', 'NorthWest', 'fontsize', 15, 'color', 'w');
box on
%% topo under:
subplot(223); 
tmp=zeros(1,32);
% tmp(plotchannel)=1;
topoplot(tmp, elocs(1:32), 'emarker2', {plotchannel, 'o','r', 10,2})
topoplot(tmp, elocs(1:32), 'emarker2', {plotchannel, 'o','r', 10,2})
% colormap(flipud(greymap))
colormap(greymap)
%% and TL component:
subplot(222); hold on;
%across across quints and chans:
dis = squeeze(mean(mean(GFX_alphavsTL_erp(:,CPchans,:,:),2),4));
stE= CousineauSEM(dis);
sh=shadedErrorBar(time_TLerp, squeeze(mean(dis,1)), stE, [], 1);
sh.mainLine.Color = useCol ;
sh.mainLine.LineWidth = 2;
sh.patch.FaceColor= useCol ;
sh.edge(1).Color= useCol ;
sh.edge(2).Color= useCol ;
xlim([-.1 .6])
ylim([-2 12])
ylabel('\muV');
xlabel('Time from target onset')
hold on
pl=plot([0 0], ylim, ['k',':'], 'linew', 2);
plot(xlim, [0 0], ['k','-']);
box on
set(gca,'ytick', [-2 12], 'xtick', [-.1 .6], 'fontsize', fs, 'color', [.95 .95 .95])

% find peak at CPP and add line
% searchlims= dsearchn(time_TLerp', [0.4 0.6]');
Yd=sh.mainLine.YData;
% mxY = max(Yd(searchlims(1):searchlims(2)));
% mxAT = find(Yd==mxY);
% pleg(1)= plot([time_TLerp(mxAT) time_TLerp(mxAT)], [0 mxY], 'color','r', 'linew', 3);
% 

% try patch?
% plot filled shading:
xts= dsearchn(time_TLerp', [.25, .55]');
xvec= [time_TLerp(xts(1)-1:xts(2)+1)]; 
yvec = [0, Yd(xts(1):xts(2)), 0];
ph = patch(xvec, yvec, 'r', 'LineStyle', 'none', 'FaceAlpha', .2);
pleg(1)= ph;

legend([pleg(1) pl],{'CPP', 'Target onset'}, 'Location', 'NorthWest', 'fontsize', 15, 'color', 'w');
box on
ylim([-2 14])



%% topo under:
subplot(224); 
tmp=zeros(1,32);
% tmp(plotchannel)=1;
topoplot(tmp, elocs(1:32), 'emarker2', {CPchans, 'o','r', 10,2})
% colormap(flipud(greymap))
colormap(greymap)
set(gcf,'color', 'w')
box on
end
%%
%%
figure(1);
legend([leg(1) leg(2) rs] ,{'Exp 1', 'Exp 2', 'RSVP onset'}, 'Location', 'NorthWest', 'fontsize', 25);
