% plotAlpha_vs_ERPfeatswrestricted
%script called from Data_explore_B_EEG.

% Having already calculated alpha per trial (previous job). Here we will
% normalize (zscore) per participant. Then show relationship to subjective
% ERPfeatures, such as P1 amplitude, and CPP strength.


% Participant level analysis
%one at a time(!)


%use feats = 1,2

TRIALrestriction = 1; %
%1 for Hit only.
%2 for miss only,
%3 for no restriction (all trials), (use this for P1).

job.concatGFX           =1; % update after TRIALrestriction changes.

TRIALstype = {'H only', 'M only', 'All Targets'};
normon=0; % change to normalize per ppant., see script for different methods.


job.plotPpantlevel      =0;
job.fitFeatsPpantlevel  =0; % fits linear ansd quad to features, per channel.

job.plotGFX             =0; % plots effects of quintiles as boxplots.
job.plotGFX_ERP         =1; % plots (target locked) ERPs by quintiles..

runplotGoF = 1; % for plotting final fit:

%% for different features, we should focus on separate channels:
ParietoOccChan = 29:32; % 
OccChan = 29:32; %  % for the P2 
CPchans = [10:12, 15:17, 20:22];% for Target locked.

% getelocs;
% for different plots of the data (Target perceived-as ? classes):

usingPhase=0;
% will convert the above to AUC.
useAUC=0;
AUCtrue = 1; % what are the truth values for calculating AUC.
%1 for target present, 0 for target absent.


%also plots group GFX at end
figuredir= '/Volumes/MattsBackup (2TB)/Frontiers Project/Documents/Figures';

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
quinsare = {'Alpha'}; 
usequin=1;
barcol = [.5 .5 .5];

RespColours = brewCOLOURS; % load colours.
alphaCols = flipud(squeeze(RespColours(4,:,:))); % purples in ascending contrast.
%%
%for each experiment
featsare = {'P1 amplitude [\muV]', 'P2 amp', 'SSVEP amplitude [\muV]', '2f1 phase-amp', 'Targlocked (.35 .45)'};

for idatatype=1%1:2
    xlabis = xlabsare{idatatype};
    
    % concat across ppants.
    if job.concatGFX           ==1
        %% for each participant,
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        
        GFX_alphavsERPfeats = [];%zeros(length(nfiles),length(featsare), 32, 5) ; % splitting into quintiles x chans
        %dims are Alpha, N1, P1, SSVEPf1, SSVEPf2
        
        GFX_alphavsTL_erp = [];%zeros(length(nfiles), 32, 426, 5);
        for ippant =1:length(nfiles);%[1:5,7:9]%
            %%
            cd(loaddir);
            %plot attention and state ratings. The process is:
            load(nfiles(ippant).name);%
            
            %store together
            plotOutput = [];%zeros(5,32,5); % features, chans, by feature.
            
            % exctract relevant data per participants, f
            if TRIALrestriction <3                
            plotOutput(1,:,:) = squeeze(P1_bySDTandAlpha(:,TRIALrestriction,:));
            plotOutput(2,:,:) = squeeze(P2_bySDTandAlpha(:,TRIALrestriction,:));
            plotOutput(3,:,:) = squeeze(f1_bySDTandAlpha(:,TRIALrestriction,:));
            plotOutput(4,:,:) = squeeze(Twof1_bySDTandAlpha(:,TRIALrestriction,:));
%             plotOutput(5,:,:) = TL_P3_byalpha;
            plotOutput(5,:,:) = squeeze(TL_P3_bySDTandAlpha(1:32,TRIALrestriction,:));
            
            GFX_alphavsTL_erp(ippant,:,:,:)= squeeze(TL_bySDTandAlpha_erp(1:32,:,TRIALrestriction,:));
            elseif TRIALrestriction==3 % use all trial version
            
           plotOutput(1,:,:) = P1_byalpha;
            plotOutput(2,:,:) = P2_byalpha;
            plotOutput(3,:,:) = f1_byalpha;
            plotOutput(4,:,:) = Twof1_byalpha;
            plotOutput(5,:,:) = TL_P3_byalpha(1:32,:);
%             plotOutput(5,:,:) = TL_byalpha;
            
            GFX_alphavsTL_erp(ippant,:,:,:)= TL_P3_byalpha_erp;
            end
            
            
            %add a mean value for TL?
            
            % normalize output?
            if normon==1
                normalize_ERPfeats
%                 normdfeats = output
                plotOutput =normdfeats;
            end

            
            GFX_alphavsERPfeats(ippant,:,:,:) = plotOutput;
            


        end
        
        save(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsERPfeats','GFX_alphavsTL_erp','time_TLerp','-append');
    end
    
    cmapB=cbrewer('qual', 'Set2', 5);
    if     job.plotPpantlevel ==1
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsERPfeats');
        
        for ippant =1:size(GFX_alphavsERPfeats,1)
            %%
            %ppant level data.
            plotOutput = squeeze(GFX_alphavsERPfeats(ippant,:,:,:));
            
            
            %% prepare for plotting:
            cd(figuredir)
            TargClass = TRIALstype{TRIALrestriction};
            
            figure(1); clf;
            set(gcf, 'units', 'normalized', 'position', [0 0 1 .5], 'color', 'w'); shg
            for ifeat = [1,5]%
                
                if ifeat<3
                    usechans = 29:32;
                    
                elseif ifeat>2 && ifeat < 5
                    usechans =ParietoOccChan;
                else
                    usechans =CPchans;
                end
                
                subplot(2,5,ifeat);
                %restrict to POchans.
                tmp = squeeze(nanmean(plotOutput(ifeat, usechans, :),2));
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
                xlabel({['alpha quintiles'];['(low-high)']});
                title(ptitle)
%                 ylabel('zscored');
%                 ylim([0 2])
                curax = get(gca);
                limits = max( abs(curax.YLim) );  % take the larger of the two "nice" endpoints
                
                 subplot(2,5,ifeat+5);
                 %which channels are retained?
                 tmpT = squeeze(plotOutput(ifeat, :, :));
                 try topoplot(mean(tmpT,2), elocs32, 'emarker2', {usechans 'o' 'w', 10, 3})
                 catch
                     disp(['no topoplot info for feature ' num2str(ifeat)])
                 end
            end
            %%
            sgtitle(['Participant ord ' num2str(ippant) ', ' TargClass ', using prominence']);
           
            
            %%
            
            %%
            print('-dpng', ['Participant ord ' num2str(ippant) ', alpha x ERPfeats ' xlabis '-Exp, ' TargClass ', using prominence']);
            if ippant==12
                pause
            end
        end
    end % ppant plot job.
    
    if job.fitFeatsPpantlevel==1
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsERPfeats');
        [nppants, nfeats, nchans, nquins] = size(GFX_alphavsERPfeats);
        
        allFeatFits= zeros(2, nppants, nfeats, nchans); % first dim linear vs quad.
        %%
        clf
        for ifeat = 1:size(plotOutput,1)
            for ippant =1:size(GFX_alphavsERPfeats,1)
                %%
                %ppant level data.
                plotOutput = squeeze(GFX_alphavsERPfeats(ippant,:,:,:));
                
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
            topoplot(mF, elocs32), title(['lin fit for ' featsare{ifeat}]);
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

            topoplot(mF, elocs32, 'pmask',allel), title(['quad fit for ' featsare{ifeat}]);
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
        
        runplotGoF=1;
        normOUTPUT=0;
        TargClass = TRIALstype{TRIALrestriction};
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        %%
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsERPfeats');
        %% after all ppants, plot GFX:
        
        figure(1); clf
        nfeats=  length(featsare);
        fcount=1;
        for ifeat = [1,5]
            
         
            if ifeat<3
               
                usechans = [29:32]; % Occ chans.
%                 usechans=29%30:32;
            elseif ifeat>2 && ifeat < 5
                usechans =ParietoOccChan;
%              
            else
                usechans =CPchans;
%              
            end
            
            
            
               %% ANOVA topo of effects by quintile:

%                         ANOVAERPfeats;
            %%
            figure(1); clf
            set(gcf, 'units', 'normalized', 'position', [0 0 .5 .5], 'color', 'w'); shg
            subplot(1,2,fcount);
            %mean topoplot of the effect first.
            tpPxCh = squeeze(nanmean(GFX_alphavsERPfeats(:, ifeat, :, :),4));
            %topoplot with subchans.
%             maskat = zeros(1,32);
%             maskat(usechans)=1;
%             tp=topoplot(squeeze(nanmean(tpPxCh,1)), elocs32, 'emarker2', {usechans '*' 'w', 10, 3}, 'pmask', maskat);
%             c=colorbar;
%             ylabel(c, 'normalized amplitude');
%             caxis([0 1])
            title(featsare{ifeat})
            
            % plot the bar.
            
            barD = squeeze(nanmean(GFX_alphavsERPfeats(:, ifeat, usechans,:),3));
            %%
            if normOUTPUT==1
                %normalize, per ppant, by dividing be mean across bins.
                tmpBar = barD;
                %         pmax = max(tmpBar,[],2);
                pmax = nanmean(tmpBar,2);
                tmpbarn = tmpBar ./ (repmat(pmax, 1, size(tmpBar,2)));
                tmpBar=tmpbarn;
                barD=tmpBar;
            end
            %
            
            barMean = squeeze(nanmean(barD,1));
            stE = CousineauSEM(barD);
            % add colours to bar.
            vec=1:5;
            nansfill =  nan(size(barD));
            
            
            
            for ib=1:5
            
                % for all but this column, plot nan.
                tmp= nansfill;
                tmp(:,ib) = barD(:,ib);
            bh=bar(nanmean(tmp)); hold on;
            bh.FaceColor = alphaCols(ib,:);
            
            end
            eb=errorbar(1:size(barD,2), barMean, stE, 'k','LineStyle', 'none');%
            %                 bh.FaceColor = cmap(ifeat,:);
            set(gca, 'fontsize', 15)
            
            %% will add best fit (linear or quadratic).
            if runplotGoF==1
                statsOUT=plotGoF(barMean);
            end
            % fix labels (after poly plot)
            %                 ptitle = ['norm(' featsare{ifeat} ')'];
            ptitle = [ featsare{ifeat}];
            xlabel({['Alpha quintiles'];['(low-high)']});
            ylabel(ptitle)
            legend off
            curax = get(gca);
            limits = max( abs(curax.YLim) );  % take the larger of the two "nice" endpoints
    
            fcount= fcount+1;
        end
        %%
%         sgtitle(['Group, ' TargClass]);
        set(gca, 'fontsize', 15)
        
        cd(figuredir)
        set(gcf, 'color', 'w')
        %%
%         print('-dpng', ['GFX, alpha x ERP features ' xlabis '-Exp, ' TargClass]);
    end
    
    
    
    
    if job.plotGFX_ERP         ==1 % plots erps effects.
  %%
        TargClass = TRIALstype{TRIALrestriction};
        
        
        cd(eegdir)
        cd( [xlabis ' alpha data']);
      
        %%
        load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_alphavsTL_erp', 'time_TLerp');
        %%
        GFX_QUINvsTL_erp = GFX_alphavsTL_erp;
        plotERPbyquintilesplits; % helper function, same one is called to plot ERPs x Behaviour.
%%
         print('-dpng', ['GFX, alpha x Targlocked ERP by quintile detrend ' num2str(detrON) ' ' xlabis '-Exp, ' TargClass]);
    %%
    end
% end
end

