% comparePhase_andBEH

%script called from Data_explore_B_EEG.

% Having already calculated alpha per trial (power and complex values) here
% we will test for differences in phase between conditions.

% 
% % Participant level analysis
% %one at a time(!)
% job1.plotIndividualAlpha_byAttention=0;% plot alpha bins by attention ratings.
% 
% job1.plotIndividualAlpha_byXrating=0; % plot alpha bins vs xaxis ratings (visibility or confidence).
% 
% job1.plotIndividualAlpha_byDetectionAccuracy =1;

% DETECTclass=1; %; output as accuracy (1), HRr (2), FAr (3), Crit (4), or d' (5)

%note that detection can be subdivided into objectively present or absent
%classes.
TARGtype = 1; %use targets when objectively present, absent, or all targets, 

topPower=1; % use all trials (0), or top half based on median split of power (1).

addGoF = 0; % for plotting final fit:

normON=0;

getelocs;
    %% prep directories
    
    %also plots group GFX at end
    figuredir = '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures/Alpha and BEH';
    
    %set up directories
    % connect to behavioural data:
    cd(homedir)
    cd('Exp 2 and 3 mat files')
    behdatadir = pwd;
    
    cd(homedir);
    cd('Exp 2 and 3 processed EEG');
    eegdir = pwd;
    %set up plot output.
    getelocs;
    xratings = {'Confidence' , 'Visibility'};
    xlabsare = {'Conf', 'PAS'};
    %%
for usewin = 1
%phasedimsare = 1s pre RSVP, 1s RSVP, 200ms preRSVP 200mspreTarget};

% for different plots of the data (Target perceived-as ? classes):

%(Perceived Present, Absent, BOTH). (1,2,3)
% H,M,CR, FA % (4,5,6,7)

    
    %for each experiment
    for idatatype=1%:2
        xlabis = xlabsare{idatatype};
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        
        %% for each participant,
        
%         GFX_alphavsBEH = nan(length(nfiles), 2,nquintiles) ; % splitting into quintiles.
%         GFX_alphaTOPO = nan(32, length(nfiles));

for ippant =1:length(nfiles)
            %%
            cd(loaddir);
            clearvars 'complex*'
            %plot attention and state ratings. The process is:
            load(nfiles(ippant).name, 'complexEEG', 'trial_table' ,'ppantnum', 'phasedimsare', 'time_secs', 'AlphaMean');
            
            %% before continuing, restrict to just a subset of trials based on high power, if required.
            if topPower==1
                %median split based on power.
                cutoff =median(AlphaMean);
                powlimitd= find(AlphaMean>=cutoff);
                
                %restrict all data to only those relevant trials:
                complexEEGin = complexEEG(:,:,powlimitd);
                trial_tablein = trial_table(powlimitd,:);
            else 
                %no restriction based on power.
                complexEEGin = complexEEG;
                trial_tablein = trial_table;
                
                
            end
            
            %%
            % compare the phase angles for all hits and misses.
            
            trialOutcomes =  trial_tablein.Outcome;
            
            allH = find(trialOutcomes==1);
            allM = find(trialOutcomes==2);
            allCR = find(trialOutcomes==3);
            allFA = find(trialOutcomes==4);
           
            useTrials(1).c = allH;
            useTrials(2).c = allM;
            useTrials(3).c = allCR;
            useTrials(4).c = allFA;
            
            %plot the distribution of phase angles for each:            
            usechan=29; % POz
            tmp_Ad = squeeze(complexEEGin(usechan,usewin,:));
            
            %% 
            clf;
            for iOutcome=1:4
            % collect phase angles (in radians), of all trials, from complex plane:
            
            trial_angles = angle(tmp_Ad(useTrials(iOutcome).c));
            
            %convert to degrees,
%             trial_angles = rad2deg(trial_angles);            
            %fit to circle,             
%             trial_angles = wrapTo360(trial_angles);
                  
            
            subplot(2,2,iOutcome)
%             [t,r]=rose(trial_angles);
%             pp= polar(t,r);
%             pp.Color = usecol;
%             pp.LineWidth = 3;
            
                circ_plot(trial_angles,'hist',[],20,true,true,'linewidth',2,'color','r')

            end
            shg
            %%
%            end
           for ichan = 1:32
             tmptrials = squeeze(complexEEGin(ichan, usewin, :));
                
                        ITPCt=tmptrials./abs(tmptrials); %divide by amp to make unit length
                        ITPCt=(sum(ITPCt,1)); % sum angles in complex plane,
                        ITPC=squeeze(abs(ITPCt)./N); %norm. average of these                                                  
             GFX_alphaTOPO(ichan, ippant) = ITPC;
           end
            
           
           
            trialOutcomes =  trial_tablein.Outcome;
            
            %what beh data will we compare for plotting?
            if job1.plotIndividualAlpha_byAttention==1
                tmp_Behdata= trial_tablein.AttResp;
                %-100 to 100 for both datasets. to avoid problems, recode to
                %all positive (for attention):
                tmp_Behdata = tmp_Behdata+100;
                
                compwas = 'Attention';
                usecol =[1,0,0];
            elseif job1.plotIndividualAlpha_byXrating==1
                tmp_Behdata= trial_tablein.Resp;
                compwas = xratings{idatatype};
                usecol = [0,0 1]; % for plotting.
                
                
            elseif job1.plotIndividualAlpha_byDetectionAccuracy==1
                % we want to plot the detectiona accuracy. so use as input
                % data, the response (H,M,FA,CR).
                tmp_Behdata = trialOutcomes; 
                compwas= 'Accuracy';
                usecol = [.6 .6 .6];
            end
            
            
            
            
            % we might want to restrict our analysis only to cases when targets
            % were perceived as present, absent, or both.
            % True target classes:
            
            TargPresent= trial_tablein.TargPresent;
              %% note that there are subclasses for detection type.
                switch  TARGtype
                    case 1 % present targets only.
                    restrange = find(trial_tablein.TargPresent==1);
                    TargClass = 'Objectively present only';
                    
                    case 2 % objectivel absent only.
                    restrange = find(trial_tablein.TargPresent==0);
                    TargClass = 'Objectively absent only';
                    
                    case 3
                        restrange = 1:length(TargPresent);
                        TargClass = 'All targets';
                                       
                end
            
            %%
            % Restrict range as appropriate.
            % for TargPresent, Alpha, and ratings, all used in later stages:
            
            tmpBeh_adapted = tmp_Behdata(restrange);
            TargPresent_adapted = TargPresent(restrange);
            
            %>Continue by sorting all trials into phase bins.
            %Start tby looking at just one channel:
            usechan=29; % POz
            tmp_Ad = squeeze(complexEEGin(usechan,usewin,restrange));
            
            % collect phase angles (in radians), of all trials, from complex plane:
            trial_angles = angle(tmp_Ad);
            
            %convert to degrees,
            trial_angles = rad2deg(trial_angles);            
            %fit to circle,             
            trial_angles = wrapTo360(trial_angles);
            
            %% now split phase angles into bins.
            Qp = quantile(trial_angles, nquintiles);
            cattrials=[];
            % collect index of the trials in each quintile.
            cattrials(1).t = find(trial_angles<=Qp(1));
            
            for icats = 1:(nquintiles-1)
                
                ind1= find(trial_angles>Qp(icats));
                ind2= find(trial_angles<=Qp(icats+1));
                %find members of both:
                catn=intersect(ind1,ind2);
                cattrials(icats+1).t = catn;
            end
            
            cattrials(nquintiles).t = find(trial_angles>Qp(length(Qp)));
            %%
            
            
            %% prepare data for plots
            
            if job1.plotIndividualAlpha_byDetectionAccuracy ~=1;
                % either calcuulate AUC for each quintile, or take the average
                % zscored subjective rating. Otherwise, calculate detection
                % accuracy.
                % we will look at the mean per category.
                    plotOutput = zeros(2,nquintiles);
                    
                    %% take z scores of all behavioural ratings also
                    tmp_Behdata_z = zscore(tmpBeh_adapted);
                    
                    for icat = 1:length(cattrials)                        
                        plotOutput(1,icat) = nanmean(tmp_Behdata_z(cattrials(icat).t));
%                         plotOutput(1,icat) = nanmean(tmpBeh_adapted(cattrials(icat).t));

                    end
                    
                    ptitle = ['z (' compwas ')'];
                    %add AUC and overwrite
                if useAUC==1
                    storeAUC=zeros(2,length(cattrials));
                    
                    
                    
                        if idatatype==1 && restricttoTargPres_Abs_All==2 % absent
                            %change behavioural data (conf values), to increase with absent cases.
                            tmpBeh_adapted = abs(tmpBeh_adapted);
                        end
                        
                    for icat= 1:length(cattrials)                        
                        Datanow = tmpBeh_adapted(cattrials(icat).t,:);
                        Targnow = TargPresent_adapted(cattrials(icat).t,:);
                        try
                            [x,y, ~, AUC]= perfcurve(Targnow, Datanow, AUCtrue);
                            
                        catch
                            % if no neg classes:
                            if sum(Targnow)==length(Targnow)
                                AUC=1; % perfect classification.
                            else
                                error('check code');
                            end
                        end
                        
                        %store AUC                        
                         %store AUC, 
                            plotOutput(1,icat) = AUC;
                              
                    end
                    ptitle = ['AUC'];
                else 
                    
                end
                
            else % we need to calculate detection rates for quintiles
                %%
                plotOutput = zeros(2,nquintiles);                
                for icat = 1:length(cattrials)
                    tmptrials = cattrials(icat).t;
                    %%
                    %calculate detection, (H+CR)/ (H+M+CR+FA);
                    TP = length(find(tmpBeh_adapted(tmptrials)==1)); %(Hits= true positive)
                    Mc = length(find(tmpBeh_adapted(tmptrials)==2));
                    TN = length(find(tmpBeh_adapted(tmptrials)==3)); % (CR =true negative)
                    FAc = length(find(tmpBeh_adapted(tmptrials)==4));
                    %%
                     %% also HRr and FAr for SDT metrics.
                    %HRr
                       HRr= TP/(TP +Mc) ;
                    %FAr
                        FAr= FAc/(FAc+TN);
                        
                    %true class (all targ present)
                    P= sum(TargPresent_adapted(tmptrials));
                    %true negative class (all targ absent)
                    N=length(find(TargPresent_adapted(tmptrials)==0));
                    
                    AccuracyTmp= (TP+TN)/(P+N);
                    if DETECTclass==1 % all trials, accuracy
                    %calculate accuracy, (H+CR)/ (H+M+CR+FA);                    
                    %%
                    %true class (all targ present)
                    P= sum(TargPresent_adapted(tmptrials));
                    %true negative class (all targ absent)
                    N=length(find(TargPresent_adapted(tmptrials)==0));
                    
                    AccuracyTmp= (TP+TN)/(P+N);
                    ptitle = ['accuracy'];
                    
                    elseif DETECTclass==2 % calculate HR
                        AccuracyTmp= HRr;
                        ptitle = ['Hit rate'];
                        
                    elseif DETECTclass==3 % calculate FAr
                        %compute
                        AccuracyTmp= FAr;
                        ptitle = ['FA Rate'];
                        
                    elseif DETECTclass==4 % calculate criterion.
                        
                       
                        AccuracyTmp= -0.5*(norminv(HRr)+ norminv(FAr));
                         
                        ptitle = ['criterion'];
                        
                    elseif DETECTclass==5 % calculate dprime
                        
                        AccuracyTmp= norminv(HRr)-norminv(FAr);
                         
                        ptitle = ['d-prime'];
                    end
                    
                    plotOutput(1,icat) = AccuracyTmp;
                    %                     plotOutput(2,icat) = nanmean(tmp_Ad_z(tmptrials)); % alpha for plot
                     %%
              
                end
            end
            
            %replace any inf with NAN.
            rem= isinf(plotOutput);
            plotOutput(rem) = nan;
            GFX_alphavsBEH(ippant,:,:) = plotOutput;
            
            
            
            %% now before plotting, arrange to preferred phase in second row.
%                 i.e.
                %because we are working with circular stats, centre on the preferred phase:
                % align to individual preferred phase for detection.
                if AligntoPref_orDynAlign==2 % dynamic alignment
                prefP = max(plotOutput(1,:));
                centreat = dsearchn([plotOutput(1,:)]', prefP');
                
                %% store for other plots:
                %actual preferred phase
                phasedim = find(plotOutput(1,:)==prefP,1,'first');                
                
                PP_preferredPhase(ippant)= Qp(phasedim);
                
                else % use the alignment best for HITS.
%                     centreat =  PP_preferredPhaseHITS(ippant);
                    centreat =  PP_preferredPhasePRESENT(ippant);
                    prefP = plotOutput(1,centreat);
                end
                
                
                %easiest way is to repeat the array.
                reparr = repmat(plotOutput(1,:), [1,3]);
                newc = find(reparr==prefP);
                %now rearrange
                shoulder  = floor(nquintiles/2); % min values before centre value.
                st = find(newc>shoulder,1, 'first');
                useID=newc(st);
                try
                plotOutput(2,1:shoulder) = reparr(useID-shoulder:useID-1);
                plotOutput(2,shoulder+1) = reparr(useID); % centre value
                plotOutput(2,shoulder+2:length(Qp)) = reparr(useID+1:useID+shoulder);
                catch 
                end
            
                
                
            %% prepare for plotting:
            cd(figuredir)
            
            
            figure(1);clf; 
            for isub=1:2 % phase bins then aligned...
                subplot(2,2,isub);
                bh=[];
                bh=bar(plotOutput(isub,:));
                title({[phasedimsare{usewin} ' , ' TargClass]})
                legend(compwas);
                bh.FaceColor = usecol;
                if restricttoTargPres_Abs_All ==2
                    bh.FaceAlpha = .2;
                end
                set(gca, 'fontsize', 15, 'XTickLabel', round(Qp))
                
                xlabel({['alpha phase angle [deg]'];['']});
                
                ylabel(ptitle)
                hold on
                
                sgtitle(['Participant ' num2str(ppantnum) ', ' TargClass]);

                
                if isub==2
                    legend('aligned alpha phase ')
                    xlabel('Aligned alpha phase')
                    set(gca, 'XTickLabel', {'-pi', '','','','','','','','','','pi'})
                   bh.FaceColor='flat';
                    bh.CData(6,:) = [1,1,1];
                end
                set(gca, 'fontsize', 15)
                ylabel(ptitle)
%                  ylim([.4 1])
            end
            %%
            %             if useAUC==1 && strcmp(compwas, 'Attention')
            %                 ylim([.45 .6])
            %
            
            set(gca, 'fontsize', 15)
            if addGoF==1
                %% test linear or quad to fit to data.
                plotGoF(plotOutput(1,:));
            end
            % fix labels (after poly plot)
            xlabel({['alpha phase bins']}); ylabel(ptitle)
            
            
            set(gcf,'color', 'w', 'units', 'normalized', 'position', [0 .4 .25 .65]); shg
            %%
            %finish with topoplot
             subplot(2,2,3);
            topoplot(GFX_alphaTOPO(:,ippant), elocs32, 'emarker2', {ParietoOccChan, '*' , 'w', 15, 3});
            c= colorbar; ylabel(c, 'Alpha ITPC')
            caxis([0 .3])
            colormap('magma');
            set(gcf, 'color' ,'w')
            
            
            set(gca, 'fontsize', 15)
            shg
            %%
             
       
        
        
        %%    
            
            
            
            print('-dpng', ['Participant ' num2str(ppantnum) ', alpha phase x ' ptitle  ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass]);
            
            %replace any inf with NAN.
            rem= isinf(plotOutput);
            plotOutput(rem) = nan;
            GFX_alphavsBEH(ippant,:,:) = plotOutput;
        end
        
        
        %% after all ppants, plot GFX:
        
        figure(2); clf;
        set(gcf,'color', 'w', 'units', 'normalized', 'position', [0 .4 .25 .65]); shg
        %
        for isub=1:2
            %
            subplot(2,2,isub);
            set(gca, 'fontsize', 15)
            
            barMean = squeeze(nanmean(GFX_alphavsBEH(:,isub,:),1));
            stE = CousineauSEM(squeeze(GFX_alphavsBEH(:,isub,:)));
            
            bh=bar(barMean); hold on;
            eb=errorbar(1:length(Qp), barMean, stE, 'k','LineStyle', 'none');%
            title({[phasedimsare{usewin} ' , ' TargClass]})
            legend(compwas, 'autoupdate', 'off');
            
            bh.FaceColor = usecol;
            if restricttoTargPres_Abs_All ==2
                bh.FaceAlpha = .2;
            end
            xlabel({['alpha phase bins'];[]}); ylabel(ptitle)
            hold on
            
            shg
            if isub==2
                xlabel('Aligned alpha phase')
%                 set(gca, 'XTickLabel', {'-pi', '','','','','','','','','','pi'})                
                bh.FaceColor=usecol; 
                bh.CData(6,:) =[.9,.9,.9];
                bh.FaceColor='flat';
                 bh.CData(1:5,:)= repmat(usecol,5,1);
                 bh.CData(7:11,:)= repmat(usecol,5,1);
                 
                
            end
            axis tight
            set(gca, 'fontsize', 15)
% axis tight
        end
        %
        if useAUC==1 && strcmp(compwas, 'Attention')
            ylim([.45 .9])
        elseif useAUC==1
            ylim([ .5 .8])
        end
        set(gca, 'fontsize', 15)
        %
        axis tight
        if addGoF==1
            %% test linear or quad to fit to data.
           plotGoF(barMean);            
            % fix labels (after poly plot)
            xlabel({['alpha quintiles'];['(low-high)']}); ylabel(ptitle)
        end
         %
        subplot(2,2,3);
        topoplot(nanmean(GFX_alphaTOPO,2), elocs32, 'emarker2', {usechan, '*' , 'w', 15, 3});
        c= colorbar; ylabel(c, 'Alpha ITPC')
        caxis([0 .3])
        colormap('magma');
        set(gcf, 'color' ,'w')
        
        
        set(gca, 'fontsize', 15)
        shg
       %
        subplot(224);
        [t,r]=rose(PP_preferredPhase);
       pp= polar(t,r);
       pp.Color = usecol;
       pp.LineWidth = 3;
% polarhistogram(
        
        title('Preferred phase across participants')
        
        [pval, z] = circ_rtest(PP_preferredPhase', repmat(28,[1,length(nfiles)]));
        shg
        
      
        %%
        print('-dpng', ['GFX, alpha phase x ' ptitle ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass]);
        
        %% final cheeky plot:
        figure(3); set(gcf, 'position', [0 .7 .25 .4])
        tmp1=squeeze(GFX_alphavsBEH(:,2,:));
%         tmp1(:,shoulder+1)= nan;
        
        %norm per pp
        if normON==1
            pM = squeeze(nanmean(tmp1,2));
            pMD = repmat(pM, [1, size(tmp1,2)]);
            tmp1adj = tmp1./pMD;
            tmp1=tmp1adj;
        end
        
        
        stE= CousineauSEM(tmp1);
        eh=errorbar(1:size(tmp1,2), nanmean(tmp1,1), stE, 'LineWidth', 3, 'LineStyle', [':']);
%         eh.Color = 'm';
        eh.Color = usecol;
        eh.MarkerFaceColor = 'k';
        
        xlabel('Aligned alpha phase')
                set(gca, 'XTickLabel', {'-pi','' ,'','','pi'})
                set(gca, 'fontsize', 20);
                ylabel(ptitle)
        xlim([0 length(Qp)+1]);
        if job1.plotIndividualAlpha_byDetectionAccuracy==1            
        ylim([.82 .96])
        ylim([.78 .96])
        else
            ylim([-.5 .5])
        end
        if useAUC==1
            ylim([.6 1])
        end
        title({[phasedimsare{usewin} ' , ' TargClass]});
        shg
%         axis tight
        %%
        print('-dpng', ['GFX, alpha phase x ' ptitle ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass ', clean']);
        %%
        %save for stats outside of matlab.
        %convert to table: 
%         tmp1=tmp1([1:8,10:12],:);
        %remove centre column of nans:
        testdata = tmp1(~isnan(tmp1));
        testd=reshape(testdata, [size(tmp1,1), size(tmp1,2)-1, 1]);
        
% output as matrix:
TableforExport= splitvars(table(testd));
%name columns for stats:
cyclePoint = [ones(1,5),ones(1,5)*2]; 
phaseDist = [fliplr([1:5]),[1:5]];

for icol = 1:size(TableforExport,2)
    varN = ['PhaseDist' num2str(phaseDist(icol)) '_cyclep' num2str(cyclePoint(icol))];
    
         TableforExport.Properties.VariableNames(icol)={varN};

end
%%       
%alternateively create columns for other outuput types:

%         %prep for output in columns:
%         dataCol = reshape(testd, [length(testd(:)),1]);
%         subsCol = repmat([1:size(testd,1)]', [size(testd,2),1]);
%         phaseDist = [fliplr([1:5]),[1:5]];
%         PhaseDisttmp= repmat(phaseDist, [size(testd,1),1]);        
%         phaseDistCol =PhaseDisttmp(:); % unwrap into column.
%         %pre or post  aligned phase, test for interaction
%         cyclePoint = [ones(1,5),ones(1,5)*2]; 
%         cycletmp = repmat(cyclePoint, [size(testd,1),1]);        
%         cycleCol = cycletmp(:); % unwrap into column.
%         %create table:
%         tabOUT = [dataCol, subsCol, phaseDistCol, cycleCol];
%         TableforExport= splitvars(table(tabOUT));
        %%
%         TableforExport.Properties.VariableNames(1)={'Data'};
%         TableforExport.Properties.VariableNames(2)={'Ppant'};
%         TableforExport.Properties.VariableNames(3)={'PhaseDist'};
%         TableforExport.Properties.VariableNames(4)={'CyclePos'};
        
        %%
       cd('/Users/mdavidson/Desktop')    
       cd('JASP stats')
    writetable(TableforExport, ['GFX, alpha phase x ' ptitle ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass ', clean.csv'])
%         save(['GFX, alpha phase x ' ptitle ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass ', clean'], 'tmp1')


% %%
%         cd(eegdir)
%         cd( [xlabis ' alpha data']);
%         PP_preferredPhasePRESENT = PP_preferredPhase;
%         save('PPpreferredphase' , 'PP_preferredPhasePRESENT', '-append');
    end

end %usewin
  %%
        plotSINEphasebyBEH;