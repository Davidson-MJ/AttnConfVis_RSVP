function plotAlpha_mainresults(homedir)
% stripped back script to just plot alpha x P1, and alpha x CPP
% see plotAlpha_vs_ERPfeats for the details. or to plot
% ppant level effects, fits across pp's,


% plot P1, and CPP, per exp, then together.

cd(homedir)
eegdir = [homedir filesep 'Exp 2 and 3 processed EEG'];
%% House keeping:
runplotGoF = 1; % for plotting final fit:
normON=0;


% for different features, we should focus on separate channels:
OccChan = [30,31,32];%[29:32]; %  % for the P2

% CPchans= [11,15:17,21];
CPchans= [15:17,20:22];%, 25:27];

elocs=getelocs(2); % chan loc data for topoplots.

% checkwins= [.25 .55]; % in sec, for CPP calc
checkwins= [.35 .65]; % in sec, for CPP calc
% load colours.
RespColours = brewCOLOURS;
alphaCols = flipud(squeeze(RespColours(4,:,:))); % purples in ascending contrast.

%grey map:
grad= linspace(1,.9, 100);
greymap= [ grad',grad',grad'; flipud(grad'), flipud(grad'),flipud(grad')];


xratings = {'Confidence' , 'Visibility'};
xlabsare = {'Conf', 'PAS'};
quinsare = {'Alpha'};
usequin=1; % separate by alpha quintiles
featsare = {'P1 [\muV]', 'CPP [\muV]'};

figure(3); clf
set(gcf, 'units', 'normalized', 'position', [0 0 .5 1], 'color', 'w'); shg
figure(4); clf
set(gcf, 'units', 'normalized', 'position', [0 0 .5 1], 'color', 'w'); shg

fs=20; % fontsize
%% Concat GFX across ppants.

fcount=1; % plot count


for idatatype=1:2; % each exp.
    xlabis = xlabsare{idatatype};
    cd(eegdir);
    cd( [xlabis ' alpha data']);
    
    nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
    loaddir = pwd;    
   
    GFX_alphavsP1=[]; % we use all trials for P1,
    GFX_alphavsTL_erp = []; % we use only targ present trials for CPP.


    for ippant =1:length(nfiles);%[1:5,7:9]%
        %%
        cd(loaddir);
        %plot attention and state ratings. The process is:
        load(nfiles(ippant).name, 'TL_bySDTandAlpha_erp', 'P1_byalpha', 'time_TLerp');%
        
        %store together
        TRIALrestriction= 1; % Hits only.
        GFX_alphavsTL_erp(ippant,:,:,:)= squeeze(TL_bySDTandAlpha_erp(1:32,:,TRIALrestriction,:));
        
        TRIALrestriction= 3; % all trials.
        GFX_alphavsP1(ippant,:,:) = P1_byalpha; % all trial version
        
        
    end
    
    
    
    % plot P1
    figure(3);
        
    usechans = OccChan;    
    %prep data, mean over chans.
    barD = squeeze(nanmean(GFX_alphavsP1(:, usechans,:),2));    
    %^^^^^^^^^^^^^^
    %Plot! 
    subplot(2,2,fcount);
    my_barplot(barD, runplotGoF, normON);
    
    ylabel(featsare{1})    
    legend off       
    title(xlabis) % show which exp.
%     ylim([4 7]);
    fcount= fcount+1;
   
    %% ^^^^^^^^^^^^^
    % Now plot the CPP bar 
    usechans=CPchans;
    
    figure(3);
    dataIN=GFX_alphavsTL_erp;
    tav = dsearchn(time_TLerp', [checkwins(1) checkwins(2)]');
    tmpB = squeeze(nanmean(dataIN(:, usechans,:,:),2));
    
    tmpBar = squeeze(nanmean(tmpB(:, tav(1):tav(2),:),2));
    
    subplot(2,2,fcount);
    my_barplot(tmpBar, runplotGoF,normON);    
    ylabel([featsare{2}])
  
    title(xlabis) % wshow which exp.
%     ylim([3 8]);
    fcount=fcount+1;
     %%
    
    % now prep the CPP (ERP)
    %^^^^^^^^^^^^^^
    usechans= CPchans;

    figure(5)
    subplot(1,2,idatatype)
    % plot the patch we average over:
    xv = [checkwins(1), checkwins(1), checkwins(2), checkwins(2)];
    yv = [-2, -1.7, -1.7, -2];
    ph=patch(xv, yv, [.9 .9 .9]);
    ph.FaceAlpha = .5;
    %prep data:
    
    for iquin = 1:size(dataIN,4)
        
        tmpD=squeeze(nanmean(dataIN(:,usechans,:,iquin),2));
        
        mP = squeeze(nanmean(tmpD,1));
        stE = CousineauSEM(tmpD);
        
        sh=shadedErrorBar(time_TLerp, mP, stE, [], 1);
        sh.mainLine.Color = alphaCols(iquin,:);
        sh.mainLine.LineWidth =2;
        sh.patch.FaceColor= alphaCols(iquin,:);
        sh.edge(1).Color  = alphaCols(iquin,:);
        sh.edge(2).Color  = alphaCols(iquin,:);
        sh.patch.FaceAlpha=0.2;
        hold on
        lg(iquin)= sh.mainLine;
    end
    
    title(xlabis) % wshow which exp.
    xlabel('Time from target onset');
    binlab = quinsare{usequin}(1:3);
    legend(lg, {[binlab ' 1'],[binlab ' 2'], [binlab ' 3'],[binlab ' 4'],[binlab ' 5']});%, 'autoupdate', 'off')
    ylabel( '\muV');
    shg
    
    %% rename to concat across (for both EXP results).
    if idatatype==1
        p1_exp1 = barD;
        CPP_exp1= tmpBar;
    else       
        p1_exp2 = barD;
        CPP_exp2= tmpBar;
    end
    
    %% topoBetaweights=nan(32,2);

for ichan=1:size(dataIN,2)
    tmpD = squeeze(nanmean(dataIN(:, ichan, tav(1):tav(2), :),3));
    %fits
    mD=mean(tmpD,1);
    [F1,G1] = fit([1:size(mD,2)]', mD', 'poly1');
    [F2,G2] = fit([1:size(mD,2)]', mD', 'poly2');
    topoBetaweights(ichan,1)=F1.p1;
    topoBetaweights(ichan,2)=F2.p1;
end

figure(10);% clf

%compute difference,
topoBetaweights(:,3)= abs(topoBetaweights(:,2))-topoBetaweights(:,1);
titlesare = {'beta linear', 'beta quad', 'abs(q)-lin)'};
for it = 1:3
    ploc= it+(3*(idatatype -1));
    subplot(2,3,ploc)
    topoplot(topoBetaweights(:,it), elocs(1:32));
    c=colorbar;
    title({[titlesare{it}];['time = ' num2str(checkwins(1,:))]})
    caxis([-.1 .1]);
    
end

    
    
end % idatatype (Exp)


%% ^^^^^^^^^^
% Now concat across exps, and plot combined:
% see command window for stats output.
%^^^^^^^^^^^^
%P1:
P1tot = [p1_exp1; p1_exp2];

figure(4); clf
subplot(221);
my_barplot(P1tot, runplotGoF,normON);
ylabel(featsare{1})
ylim([1 3]);

subplot(223); % plot chans used:
tmpD=squeeze(nanmean(nanmean(GFX_alphavsP1,1),3));

topoplot(tmpD, elocs(1:32), 'emarker2', {OccChan, '*', 'w', 5,4});
topoplot(tmpD, elocs(1:32), 'emarker2', {OccChan, 'o', 'k',8,3});
c=colorbar;
caxis([0 2]);
set(gca, 'fontsize', 15)
colormap('inferno');
%% 
%^^^^^^^^^^^^
%CPP:
CPPtot = [CPP_exp1; CPP_exp2];

figure(4)
subplot(222);
my_barplot(CPPtot, runplotGoF,normON);
ylim([5 10]);
ylabel([featsare{2} ]);%', ' num2str(checkwins(1)*100) ':' num2str(checkwins(2)*100) ' ms']);
%
subplot(224); % plot chans used:
% tmpD=zeros(1,32);
% topoplot(tmpD, elocs(1:32), 'emarker2', {CPchans, 'o', 'r'});
% colormap(greymap)

%or plot CPP, mean over window:
plotTopo = squeeze(mean(mean(mean(dataIN(:,:,tav(1):tav(2),:),3),4),1));
topoplot(plotTopo, elocs(1:32), 'emarker2', {CPchans, '*', 'w', 5,4});
topoplot(plotTopo, elocs(1:32), 'emarker2', {CPchans, 'o', 'k',8,3});

colormap('inferno')
c=colorbar;
ylabel(c, '\muV', 'fontsize', 15)
set(gca, 'fontsize', 15)
caxis([0 8])
set(gcf, 'color', 'w')

%% show beta topoweights:
%calculate beta per channel, plot beta weights.





function h= my_barplot(plotdata, runplotGoF,normON)

    
    %plot the same bar chart each time.

if normON
   
%    pM = squeeze(nanmean(plotdata,2));
   pM = max(plotdata,[],2);
   pDiv = repmat(pM, 1, size(plotdata,2));
   newd = plotdata ./pDiv;
   plotdata = newd;
    
end


barM= squeeze(nanmean(plotdata,1));
%w/in ppant STE
stEtmp = CousineauSEM(plotdata);

%^^^^^^^^^^^^^^

%^^^^^^^^^^^^^^
% plot separately, to add pretty colours to bar.
v=1:5;
nansfill =  nan(size(plotdata));
for ib=1:5
    % for all but this column, plot nan.
    tmp= nansfill;
    tmp(:,ib) = plotdata(:,ib);
    bh=bar(nanmean(tmp)); hold on;
    bh.FaceColor = alphaCols(ib,:);
    
end
eb=errorbar(1:size(plotdata,2), barM, stEtmp, 'k','LineStyle', 'none');%

% add best fit (linear or quadratic).
if runplotGoF==1
    statsOUT=plotGoF(plotdata);

    % print relevant params into command window:
    
    
end
xlabel({['Alpha quintiles'];['(low-high)']});
legend off
set(gca, 'fontsize', 20)

end
end