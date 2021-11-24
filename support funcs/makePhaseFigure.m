% makePhaseFigure;
clear all
setpaths;
close all

%piggy backing of calcAlphaITC_vs_Behaviour, this simply combines the data
%in one plot for the paper.

% calcAlphaITC_vs_Behaviour
%script called from Data_explore_B_EEG.

% Having already calculated alpha per trial (power and complex values) (previous job). Here we will
% normalize (zscore) per participant. Then show relationship to subjective
% report (attention , confidence ,visibility ratings, and MCOG).


% Participant level analysis
%one at a time(!)

classesare = 1:2;%1:2; % 1:3' attention, xratings, objective performance. plots both in a subplot.

%if useClass==3, which type of performance to plot?
DETECTclass=2; %; output as accuracy (1), HRr (2), FAr (3), Crit (4), or d' (5)

job1.Crunchacrossppants = 1; %after changing the above, crunch across ppants. All chans, all phase bins.

job1.plotIndividualFX=0;% plot ppant level.
job1.plotGFX = 1;

exportforJasp=1;
%note that detection can be subdivided into objectively present or absent
%classes. We will actually overlay these types.

% TARGtype = 1:3; %use targets when objectively present, absent, or all
% targets,  (now specified in for loop below).


%other specs.
topPower=0; % use all trials (0), or top half based on median split of power (1).

addGoF = 0; % for plotting final fit:

normON=0;
% what type of phase alignment
AligntoPref_orDynAlign = 2; % 1 for aligning to preferred phase for 'HITS'; 2 for dynamic (best result).



nquintiles = 11;% 11
%%
% will convert the above to AUC.
useAUC=0;
AUCtrue = 1; % what are the truth values for calculating AUC.
%1 for target present, 0 for target absent.
ParietoOccChan = 24:32;
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
usewin = 1; % time window for analysis.
%phasedimsare = 1s pre RSVP, 1s RSVP, 200ms preRSVP 200mspreTarget};

% for different plots of the data (Target perceived-as ? classes):

%(Perceived Present, Absent, BOTH). (1,2,3)
% H,M,CR, FA % (4,5,6,7)

pc=1;
%for each experiment
for idatatype=1:2
    for  useClass= classesare % attention and conf ratings.
        xlabis = xlabsare{idatatype};
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);
        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        
        if AligntoPref_orDynAlign == 2 % 1 for aligning to preferred phase for 'HITS'; 2 for dynamic (best result).
            PP_preferredPhase= zeros(length(nfiles),2); % saves per ppant, and target type.
        else
            load('PPpreferredphase.mat'); % loads previous
        end
        
        
        %         load(['GFX_' xlabis '_ratingsandalpha']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if job1.Crunchacrossppants == 1 %after changing the above, crunch across ppants. All chans, all phase bins.
            
            
            GFX_alphavsBEH = nan(length(nfiles), 2,nquintiles, 3) ; %
            
            %for reference:
            %             [nppants, alignment, phasebins, targclass] = size(GFX_alphavsBEH)
            
            %             GFX_alphaTOPO = nan(32, length(nfiles));
            %% for each participant,
            for TARGtype = 1:2
                
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
                    
                    
                    
                    %% store GFX of alpha for topoplot
                    prestim = dsearchn([time_secs]', [-.3 .7]');
                    
                    % calculate ITPC over all trials.
                    N= size(complexEEGin,3);
                    for ichan = 1:32
                        tmptrials = squeeze(complexEEGin(ichan, usewin, :));
                        
                        ITPCt=tmptrials./abs(tmptrials); %divide by amp to make unit length
                        ITPCt=(sum(ITPCt,1)); % sum angles in complex plane,
                        ITPC=squeeze(abs(ITPCt)./N); %norm. average of these
                        GFX_alphaTOPO(ichan, ippant) = ITPC;
                    end
                    
                    trialOutcomes =  trial_tablein.Outcome;
                    
                    %what beh data will we compare for plotting?
                    switch useClass
                        case 1
                            
                            tmp_Behdata= trial_tablein.AttResp;
                            %-100 to 100 for both datasets. to avoid problems, recode to
                            %all positive (for attention):
                            tmp_Behdata = tmp_Behdata+100;
                            
                            compwas = 'Attention';
                            usecol =[1,0,0; .99,.75, .79 ]; % red, pale red
                        case 2
                            tmp_Behdata= trial_tablein.Resp;
                            compwas = xratings{idatatype};
                            usecol = [0,0 1; .65, .74, .85]; % blue, pale blue
                            
                            
                        case 3 % accuracy.
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
                            TargClass = 'Target present only';
                            
                        case 2 % objectivel absent only.
                            restrange = find(trial_tablein.TargPresent==0);
                            TargClass = 'Targer absent only';
                            
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
                    usechan=31; % POz
                    tmp_Ad = squeeze(complexEEGin(usechan,usewin,restrange));
                    
                    % collect phase angles (in radians), of all trials, from complex plane:
                    trial_angles = angle(tmp_Ad);
                    
                    
                    %unwrap to smooth discontinuous phase angles (jumps larger than
                    %360deg)
                    trial_angles = unwrap(trial_angles);
                    
                    %convert to degrees,
                    %             trial_angles = rad2deg(trial_angles);
                    %fit to circle,
                    %             trial_angles = wrapTo360(trial_angles);
                    
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
                    
                    %% prepare data for plots
                    
                    if useClass<3 %i.e. not objective performance.
                        disp('calculating for z-scores of subjective response categories');
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
                                
                            elseif DETECTclass==3 % calculate FArs
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
                    %                 GFX_alphavsBEH(ippant,:,1, TARGtype) = plotOutput;
                    
                    
                    
                    %% now before plotting, arrange to preferred phase in second row.
                    %                 i.e.
                    %because we are working with circular stats, centre on the preferred phase:
                    % align to individual preferred phase for detection.
                    if AligntoPref_orDynAlign==2 % dynamic alignment
                        
                        prefP = max(plotOutput(1,:));
                        %or worst
%                         prefP = min(plotOutput(1,:));
                        
                        centreat = dsearchn([plotOutput(1,:)]', prefP');
                        
                        %% store for other plots:
                        %actual preferred phase
                        phasedim = find(plotOutput(1,:)==prefP,1,'first');
                        
                        PP_preferredPhase(ippant,TARGtype)= Qp(phasedim);
                        
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
                    
                    %% store across ppants.
                    GFX_alphavsBEH(ippant,:,:, TARGtype) = plotOutput;
                    
                    
                end % nppants
            end % targtype
        end % crunch job.
        
        if job1.plotIndividualFX==1% plot ppant level.
            
            %% prepare for plotting:
            cd(figuredir)
            
            
            figure(1);clf;
            for isub=1:2 % phase bins then aligned...
                subplot(2,2,isub);
                bh=[];
                bh=bar(plotOutput(isub,:));
                title({[phasedimsare{usewin} ' , ' TargClass]})
                legend(compwas);
                bh.FaceColor = usecol(1,:);
                
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
            
            
            print('-dpng', ['Participant ' num2str(ppantnum) ', alpha phase x ' ptitle  ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass]);
            
        end % PFX plot
        
        
        %% after all ppants, plot GFX:
        
        Targsare = {'Target present', 'Target absent','All targets'};
%         usecolsplot = [present,absent];

        if job1.plotGFX == 1
            normON=0;
            for TARGtype=2%:2
                
                figure(1); %if TARGtype ==1; clf; end
                %naming
                TargClass = Targsare{TARGtype};
                
                set(gcf,'color', 'w', 'units', 'normalized', 'position', [0 .8 .6 .65]); shg
                
                isub=2; % only plot phase aligned.
%                 subplot(2,2, TARGtype + (idatatype-1)*2);
                subplot(2,2, pc);
                set(gca, 'fontsize', 15)
                
                tmp1 = squeeze(GFX_alphavsBEH(:,isub,:, TARGtype));
                %                 norm per pp ?
                if normON==1
                    %                         pM = squeeze(nanmean(tmp1,2));;
                    pM= max(tmp1, [],2);
                    pMD = repmat(pM, [1, size(tmp1,2)]);
                    tmp1adj = tmp1./pMD;
                    tmp1=tmp1adj;
                end
                % first plot as grey, then add with centre removed.
                
                
                %                     tmp1(:,6)=nan;
                stE= CousineauSEM(tmp1);
                eh=errorbar(1:size(tmp1,2), nanmean(tmp1,1), stE, 'LineWidth', 3, 'LineStyle', [':']);
                eh.Color = [.9 .9 .9];
                hold on
                %adjust y value to include occluded range.
                ylim([ min(nanmean(tmp1,1))-.2 max(nanmean(tmp1,1))+.3])
                
                
                %occlude centre column.
                tmp1(:,6)=nan;
                eh=errorbar(1:size(tmp1,2), nanmean(tmp1,1), stE, 'LineWidth', 3, 'LineStyle', [':']);
                eh.Color = usecol(TARGtype,:);
                eh.MarkerFaceColor = 'k';
                %                     if TARGtype==2;
                %                         eh.
                
                xlabel('Aligned alpha phase')
                ylabel(ptitle)
                %                             set(gca, 'fontsize', 20);
                %                             ylabel(ptitle)
                
                set(gca, 'Xtick', 1:11,'XTickLabel', {'-pi','' ,'','','','','','','','','pi'})
                
%                 title({[phasedimsare{usewin} ' , ' TargClass]});
                title({['Experiment ' num2str(idatatype) ' ' compwas ' when ' TargClass]});
                
                set(gca, 'fontsize', 15)
                %                     axis tight
                xlim([0 length(Qp)+7]);
                xlim([0 length(Qp)+1]);
                set(gca, 'fontsize', 15)
                
              
                
%                 %% add rose plot of preferred phase per ppant.
                figure(13);  set(gcf, 'color', 'w')
%                 
                subplot(2,2, pc);
%                 %                     hold on ; %clf
% 
                nbins=20;
%                 hold on
%                 
             h=polarhistogram(PP_preferredPhase(:,TARGtype), nbins);
             h.FaceColor = usecol(TARGtype,:);
% %                     subp
  pc=pc+1;
%%
if exportforJasp==1;
testd=tmp1;
testd(:,shoulder+1)= [];

% output as matrix:
TableforExport= splitvars(table(testd));
%name columns for stats:
cyclePoint = [ones(1,(nquintiles-1)/2),ones(1,(nquintiles-1)/2)*2];
phaseDist = [fliplr([1:(nquintiles-1)/2]),[1:(nquintiles-1)/2]];

for icol = 1:size(TableforExport,2)
    varN = ['PhaseDist' num2str(phaseDist(icol)) '_cyclep' num2str(cyclePoint(icol))];
    
    TableforExport.Properties.VariableNames(icol)={varN};
    
end


%%
cd('/Users/mdavidson/Desktop')
cd('JASP stats')
writetable(TableforExport, ['GFX, alpha phase x ' ptitle ' ' phasedimsare{usewin} ' ' xlabis '-Exp, ' TargClass ', clean.csv'])
end

                
            end% TARGtype
           
            
            %%
           
            
        end % plot GFX
        
        tmpD= wrapTo360(PP_preferredPhase);
        if useClass==1
            PP_preferredphase_Attn_PresAbs = tmpD;
        else
            PP_preferredphase_Xrating_PresAbs = tmpD;
        end
        
    end % subj ratings.
    %         plotSINEphasebyBEH;
    
    %%
    %
%     save('PPpreferredphase.mat', 'PP_preferredphase_Attn_PresAbs', 'PP_preferredphase_Xrating_PresAbs');
end % datatype (3xp_

