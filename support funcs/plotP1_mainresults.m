function plotP1_mainresults(homedir)
% stripped back script to just plot  P1 x Beh ratings (Hits only)


% plot P1, 

cd(homedir)
eegdir = [homedir filesep 'Exp 2 and 3 processed EEG'];
%% House keeping:
runplotGoF = 1; % for plotting final fit:
normON=0;


% for different features, we should focus on separate channels:
OccChan = [30,31,32];%[29:32]; %  % for the P2

elocs=getelocs(2); % chan loc data for topoplots.

% load colours.
RespColours = brewCOLOURS;


xratings = {'Confidence' , 'Visibility'};
xlabsare = {'Conf', 'PAS'};
quinsare = {'Alpha'};
usequin=1; % separate by alpha quintiles
featsare = {'P1 [\muV]'};

figure(3); clf
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
   
    GFX_BEHvsP1=[]; % we use all trials for P1,
    
    

    for ippant =1:length(nfiles);%[1:5,7:9]%
        %%
        cd(loaddir);
        %plot attention and state ratings. The process is:
        load(nfiles(ippant).name, 'P1_byAttention', 'P1_byXrating');%
                
        GFX_BEHvsP1(ippant,:,:,1) = P1_byXrating; % Hits only
        GFX_BEHvsP1(ippant,:,:,2) = P1_byAttention; % Hits only ( see calcERPfeats_perQuintile).
        
        
    end
    
    
    
    
    figure(3);
    % plot P1 (xrating)
    %prep data, mean over chans.
    barD = squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,1),2));    
    %^^^^^^^^^^^^^^
    %Plot! 
    subplot(2,2,fcount);
   
    if idatatype==1
        useCols = cbrewer('seq', 'Blues', 5);
    else
        useCols = cbrewer('seq', 'Greens', 5);
    end
    
    my_barplot(barD, runplotGoF, normON);
    
    ylabel(featsare{1})    
    legend off       
    title(xratings{idatatype}) % show which exp.
    fcount= fcount+1;
   
        
    % plot P1 (Attn)
    %prep data, mean over chans.
    
    barD = squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,2),2));    
    %^^^^^^^^^^^^^^
    %Plot! 
    subplot(2,2,fcount);
    
    useCols = cbrewer('seq', 'Reds', 5);
    my_barplot(barD, runplotGoF, normON);
    
    ylabel(featsare{1})    
    legend off       
    title('Attention') % show which exp.
    fcount= fcount+1;
   
    
    %% rename to concat across (for both EXP results).
    if idatatype==1
        p1_exp1_X =  squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,1),2));    
        p1_exp1_Y=  squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,2),2));    
    else       
        p1_exp2_X =   squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,1),2));    
        p1_exp2_Y=   squeeze(nanmean(GFX_BEHvsP1(:, OccChan,:,2),2));    
    end
    
    
    
    
end % idatatype (Exp)


%% ^^^^^^^^^^
% Now concat across exps, and plot combined:
% see command window for stats output.
%^^^^^^^^^^^^
%P1:
P1tot = [p1_exp1_X; p1_exp2_X];

figure(4); clf
subplot(221);
useCols = cbrewer('seq', 'GnBu', 5);
my_barplot(P1tot, runplotGoF,normON);
ylabel(featsare{1})
xlim([0 6])
title('Conf and vis') % show which exp.
%% 
%^^^^^^^^^^^^
%CPP:
P1tot = [p1_exp1_Y; p1_exp2_Y];

figure(4)
subplot(222);
useCols = cbrewer('seq', 'Reds', 5);
my_barplot(P1tot, runplotGoF,normON);
% ylim([5 10]);
ylabel([featsare{1} ]);%', ' num2str(checkwins(1)*100) ':' num2str(checkwins(2)*100) ' ms']);
title('Attention')
xlim([0 6])
%%


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
    bh.FaceColor = useCols(ib,:);
    
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