% plotBEH_vs_ERPfeats
%script called from Data_explore_B_EEG.

%In this, script we have already calculated the mean ERPs (Target locked),
%to divisions based on x axis (Conf/Vis) and y-axis (attention).
%Here we plot the results of this criteria on target locked ERPs.


% Participant level analysis
%one at a time(!)

job.concatGFX           =1;

normon=0; % change to normalize per ppant., see script for different methods.

%use feats = 1,2
%1 for Hit only.
%2 for miss only,
%3 for no restriction (all trials), 
TRIALrestriction = 1; %
TRIALstype = {'H only', 'M only', 'All Targets'};

job.plotPpantlevel      =0;
job.fitFeatsPpantlevel  =0; % fits linear and quad to features, per channel.

job.plotGFX             =0; % plots stationary effects.
job.plotGFX_ERP         =1; % plots stationary effects.

runplotGoF = 0; % for plotting final fit:




usequin = 1; % 1 for Xratings, 2 for Attention.

fontsize = 15;

%% for different features, we should focus on separate channels:
ParietoOccChan = 29:32; % 
OccChan = 29:32; %  % for the P2 
% CPchans = [10:12, 15:17, 20:22];% for Target locked.
CPchans = [15:17, 20:22];% for Target locked.
elocs= getelocs(2);
% for different plots of the data (Target perceived-as ? classes):


% will convert the above to AUC.
useAUC=0;
AUCtrue = 1; % what are the truth values for calculating AUC.
%1 for target present, 0 for target absent.


%also plots group GFX at end
figuredir = '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures/Alpha and ERP features';

%set up directories
% connect to behavioural data:
cd(homedir)
cd('Exp 2 and 3 mat files')
behdatadir = pwd;

cd(homedir);
cd('Exp 2 and 3 processed EEG');
eegdir = pwd;
%set up plot output.

xratings = {'Confidence' , 'Visibility'};
xlabsare = {'Conf', 'PAS'};
%%
%for each experiment
featsare = {'P2(.15 .25)' , 'P3(.35 .45)'};

for idatatype=1%:2
    xlabis = xlabsare{idatatype};
    
    quinsare = {xratings{idatatype}, 'Attention'};    
    TargClass = TRIALstype{TRIALrestriction};
    % concat across ppants.
    if job.concatGFX           ==1
        %% for each participant,
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        %%
        GFX_QUINvsERPfeats = zeros(length(nfiles),length(featsare), 32, 5) ; % splitting into quintiles x chans
        %dims are Alpha, N1, P1, SSVEPf1, SSVEPf2
        GFX_QUINvsTL_erp = zeros(length(nfiles), 32, 426, 5);
        for ippant =1:length(nfiles)
            %%
            cd(loaddir);
            %plot attention and state ratings. The process is:
            load(nfiles(ippant).name);
            
            %store together
            plotOutput = zeros(2,32,5); % features, chans, by feature.
            
            if usequin==1
                %%mean  amp
%                 error('this data not completed - see calcTarglocked...wrestrctions ');
                plotOutput(1,:,:) = squeeze(TL_P2_bySDTandXresp(1:32,TRIALrestriction,:));
                plotOutput(2,:,:) = squeeze(TL_P3_bySDTandXresp(1:32,TRIALrestriction,:));
                TL_erp = squeeze(TL_bySDTandXresp_erp(1:32,:,TRIALrestriction, :));
                barcol='b';
            else
                plotOutput(1,:,:) = squeeze(TL_P2_bySDTandAttresp(1:32,TRIALrestriction,:));
                plotOutput(2,:,:) = squeeze(TL_P3_bySDTandAttresp(1:32,TRIALrestriction,:));
                TL_erp = squeeze(TL_bySDTandAttresp_erp(1:32,:,TRIALrestriction, :));
                barcol='r';
            end
            
            if normon==1
                normalize_ERPfeats
%                 normdfeats = output
                plotOutput =normdfeats;
            end

            
            GFX_QUINvsERPfeats(ippant,:,:,:) = plotOutput;
            
            GFX_QUINvsTL_erp(ippant,:,:,:)= TL_erp;
        end
        
        save(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_QUINvsERPfeats','GFX_QUINvsTL_erp','time_TLerp','-append');
    end
    
    cmapB=cbrewer('qual', 'Set2', 5);
    if     job.plotPpantlevel ==1
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsERPfeats');
        
        for ippant =1:size(GFX_QUINvsERPfeats,1)
            %%
            %ppant level data.
            plotOutput = squeeze(GFX_QUINvsERPfeats(ippant,:,:,:));
            
            
            %% prepare for plotting:
            cd(figuredir)
            TargClass = 'All targets';
            
            figure(1); clf;
            set(gcf, 'units', 'normalized', 'position', [.5 1 .5 1], 'color', 'w'); shg
            for ifeat = 1:length(featsare)
               usechans = CPchans;
                
                subplot(2,2,ifeat);
                %restrict to POchans.
                tmp = squeeze(mean(plotOutput(ifeat, usechans, :),2));
                bh= bar(tmp);
                bh.FaceColor = cmapB(ifeat,:);
                %             bh.FaceColor = cmap(ifeat,:);
                set(gca, 'fontsize', 15)
                
                %% will add best fit (linear or quadratic).
                %               if runplotGoF==1
                %               plotGoF(plotOutput(ifeat,:));
                %               end
                % fix labels (after poly plot)
                %             ptitle = ['normalized (' featsare{ifeat} ')'];
                ptitle = [ featsare{ifeat}  ];
                xlabel({[quinsare{usequin} ' quintiles'];['(low-high)']});
                title(ptitle)
                ylabel('zscored');
%                 ylim([0 2])
                curax = get(gca);
                limits = max( abs(curax.YLim) );  % take the larger of the two "nice" endpoints
                
                 subplot(2,2,ifeat+2);
                 %which channels are retained?
                 tmpT = squeeze(plotOutput(ifeat, :, :));
                 topoplot(mean(tmpT,2), elocs, 'emarker2', {usechans 'o' 'w', 10, 3})
            end
            %%
            sgtitle(['Participant ord ' num2str(ippant) ', ' TargClass ]);
          
            
            %%
            
            %%
            print('-dpng', ['Participant ord ' num2str(ippant) ', ' quinsare{usequin} ' x TargLocked ERPfeats ' xlabis '-Exp, ' TargClass 'zscored']);
            
        end
    end % ppant plot job.
    
    if job.fitFeatsPpantlevel==1
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_QUINvsERPfeats');
        [nppants, nfeats, nchans, nquins] = size(GFX_QUINvsERPfeats);
        
        allFeatFits= zeros(2, nppants, nfeats, nchans); % first dim linear vs quad.
        %%
        clf
        for ifeat = 1:size(plotOutput,1)
            for ippant =1:size(GFX_QUINvsERPfeats,1)
                %%
                %ppant level data.
                plotOutput = squeeze(GFX_QUINvsERPfeats(ippant,:,:,:));
                
                %test linear and quad fit, for alpha, at each channel.
                for ichan = 1:size(plotOutput,2)
                    tmpD = squeeze(plotOutput(ifeat,ichan,:));
                    %% will add best fit (linear or quadratic).
                    try 
                    [F1,G1] = fit([1:5]',tmpD, 'poly1');
                    [F2,G2] = fit([1:5]', tmpD, 'poly2');
                    %%
                    %Goodness of fits:
                    goodnessfitsare = [G1.adjrsquare, G2.adjrsquare];
                    %Output of Fits
                    Fitsare  = {F1, F2};%
                    allFeatFits(1, ippant, ifeat, ichan) = F1.p1;
                    allFeatFits(2, ippant, ifeat, ichan) = F2.p1;
                    catch
                        
                    end
                end
            end
            
            %% plot mean fits across ppants.
            
            subplot(211);
            mF = squeeze(nanmean(allFeatFits(1,:,ifeat, :),2));
            topoplot(mF, elocs), title(['lin fit for ' featsare{ifeat}]);
            set(gca, 'fontsize', 20); c=colorbar;
            ylabel(c,'coeff.')
            %%
            subplot(212);
            mF = squeeze(nanmean(allFeatFits(2,:,ifeat, :),2));
            allF = squeeze(allFeatFits(2,:,ifeat, :));
            pvals= [];
            for ichan = 1:size(allF,2)
                [~,pvals(ichan)] = ttest(allF(:,ichan),0);
            end            
            %mask nonsig:
            allel = ones(1,length(pvals));
            %nonsig
            ns = find(pvals<.05);
            allel(ns)=0;

            topoplot(mF, elocs, 'pmask',allel), title(['quad fit for ' featsare{ifeat}]);
            set(gcf, 'color', 'w');  
            ylabel(c,'coeff.')
            c=colorbar;
             %%
        cd(figuredir)
        set(gcf, 'color', 'w');
        set(gca, 'fontsize', 20);
        %%
        print('-dpng', ['GFX, alpha x ERP features ' xlabis '-Exp, Feature fits:' featsare{ifeat}]);
        end
    end
        
    %%
    
    if    job.plotGFX             ==1
        TargClass = 'All targets';
        TRIALstype{TRIALrestriction};
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        %%
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_QUINvsERPfeats');
        %% after all ppants, plot GFX:
        
        figure(1); clf
        nfeats=  length(featsare);
        for ifeat = 1:nfeats
            
            usechans =CPchans;
            %
%             usechans= [16,21,26,29,31];
            
            
               %% ANOVA topo of effects by quintile:

%                         ANOVAERPfeats;
            %%
            figure(1); 
            set(gcf, 'units', 'normalized', 'position', [0 1 1 1], 'color', 'w'); shg
            subplot(2,nfeats,ifeat);
            %mean topoplot of the effect first.
            tpPxCh = squeeze(nanmean(GFX_QUINvsERPfeats(:, ifeat, :, :),4));
            %topoplot with subchans.
            topoplot(squeeze(nanmean(tpPxCh,1)), elocs, 'emarker2', {usechans '*' 'w', 10, 3});
            c=colorbar;
            ylabel(c, 'normalized amplitude');
            caxis([0 1])
            title(featsare{ifeat})
            set(gca, 'fontsize', 20);
            
            %% plot the bar.
            subplot(2,nfeats,ifeat+nfeats);
            barD = squeeze(nanmean(GFX_QUINvsERPfeats(:, ifeat, usechans,:),3));
            barMean = squeeze(nanmean(barD,1));
            stE = CousineauSEM(barD);
            
            bh=bar(barMean); hold on;
            bh.FaceColor = cmapB(ifeat,:);
            eb=errorbar(1:size(barD,2), barMean, stE, 'k','LineStyle', 'none');%
            %                 bh.FaceColor = cmap(ifeat,:);
            set(gca, 'fontsize', 15)
            
            %% will add best fit (linear or quadratic).
            if runplotGoF==1
                plotGoF(barMean);
            end
            % fix labels (after poly plot)
            %                 ptitle = ['norm(' featsare{ifeat} ')'];
            ptitle = ['normd '  featsare{ifeat}];
            xlabel({[quinsare{usequin} ' quintiles'];['(low-high)']});
            ylabel(ptitle)
            
            curax = get(gca);
            limits = max( abs(curax.YLim) );  % take the larger of the two "nice" endpoints
%             ylim( [0, 1] );
        end
        %%
        sgtitle(['Group, ' TargClass]);
        set(gca, 'fontsize', 15)
        
        cd(figuredir)
        set(gcf, 'color', 'w')
        %%
        print('-dpng', ['GFX, ' quinsare{usequin} ' x ERP features ' xlabis '-Exp, ' TargClass]);
    end
    
    
    
    %%
    if job.plotGFX_ERP         ==1 % plots erps effects.
       
            %%
        cd(eegdir)
        cd( [xlabis ' alpha data']);        
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_QUINvsTL_erp', 'time_TLerp');
    %%
        plotERPbyquintilesplits;
       %% 
        print('-dpng', ['GFX,  ' quinsare{usequin} ' x Targlocked ERP by quintile detrend ' num2str(detrON) ' ' xlabis '-Exp, ' TargClass]);
    
    end
end

