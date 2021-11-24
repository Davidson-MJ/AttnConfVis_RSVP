% calcAlpha_vs_Behaviour_table
%script called from Data_explore_B_EEG.
%based on calcAlpha_vs_Behaviour.

% This script calculates various behavioural indices when split by alpha
% quintiles, and saves in a table per ppant.

%ready for plotting next jobs.


% Having already calculated alpha per trial (previous job). Here we will
% normalize (zscore) per participant. Then show relationship to subjective
% report (attention , confidence ,visibility ratings, MCOG, accuracy etc).

% which features to calc?
% job1.plotIndividualAlpha_byAttention=0;                 % plots attention ratings by alpha bins.
% job1.plotIndividualAlpha_byXrating=0;                   % plot alpha bins vs xaxis ratings (visibility or confidence).
% job1.plotIndividualAlpha_byDetectionAccuracy =1;        % select output type below:

%note that there are different options for detection output, spec below:
% DETECTclass = 1;                                        %; output as accuracy (1), HRr (2), FAr (3), Crit (4), or d' (5)
% useAUC=0;                                               % will convert the above to AUC.


%what to do with the output (on plots)?
% addGoF = 1;                                              % for plotting final fit:
% normOUTPUT = 0;                                          % normalize per ppant, to the max  across alpha bins
% useZSCORE_orRank =1;                                    % to determine alpha quintiles.

% which plots to show?
% job1.plotPPantFX = 0;                                    % plot ppant level.
% job1.plotGFX  =1;                                        % plot GFX
% job1.exporttoJasp =0;                                    % export table to JASP for further analysis.


job.calcTableperppant=1;
job.concatTable_GFX=1;


statsOUTPUT=[]; % this collects as we go,
% and outputs to command window upon completion



%%

% setpaths;
mydirs;
%
nquintiles = 5;                                          % how many to use?


%% BEGIN
cd(homedir)

%% load channe data and specify channels to plot.
% getelocs;

ParietoOccChan = 24:32; % conf,v
% ParietoOccChan = [29:32]; % attn, acc

%for which experiment
for idatatype=1:2
    
    %% load data
    
    xlabis = xlabsare{idatatype};
    
    cd(eegdir)
    cd( [xlabis ' alpha data']);
    nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
    loaddir = pwd;
    
    GFX_table = table();
    GFX_alphaTOPO=[]; % topoplot at alpha freqs.
    GFX_OccSpec =[]; %spectrum at PO chans.
    %% load ppant data, store for GFX
    for ippant =1:length(nfiles); %[1:5,7:9]%
    
        %%
        if job.calcTableperppant % calculate all behaviour per ppant?
        
        cd(loaddir);
        clearvars 'Alpha*'
        %load alpha per trial
        load(nfiles(ippant).name, 'AlphaMean_FFT', 'AlphaMean', 'trial_table' ,'ppantnum', 'ampEEG', 'ampEEG_faxis', 'time_secs');
        % alphamean using hilbEEG
        AlphaMean_chanav= squeeze(nanmean(AlphaMean_FFT(ParietoOccChan,:),1));
        GFX_alphaTOPO(:, ippant) = squeeze(nanmean(AlphaMean_FFT(1:32,:), 2));
        
        GFX_OccSpec(ippant,:) = squeeze(mean(nanmean(ampEEG(ParietoOccChan,:,:), 1),2)); % chans, trials, hz
        %alpha mean using fft
%         AlphaMean= squeeze(nanmean(AlphaMean_FFT(ParietoOccChan,:),1));
        
        dataTable = table(); % for all output
        %% store GFX of alpha for topoplot
                prestim = dsearchn([time_secs]', [-.3 .7]');
                
        
        %%  we might want to restrict our analysis only to cases when targets
        % were present, absent, or use both.
        % True target classes:
        
        %note that detection can be subdivided into objectively present or absent
        % target classes:
        % TargClasses ={'Target-present only', 'Target-absent only', 'All trials'};
        TargClasses ={'trgPres', 'trgAbs', 'trgAll'};
        TargPresent= trial_table.TargPresent;
        
        
        %% store data split by attention in a big table.
        %what beh data will we compare for plotting?
        
        for isubjectivedata=1:3 % attention, xaxis, H,M,FA,CR, accuracy, type2.
            switch isubjectivedata
                
                case 1 % yaxis rating
                    
                    tmp_Behdata= trial_table.AttResp;
                    %-100 to 100 for both datasets. to avoid problems, recode to
                    %all positive (for attention):
                    tmp_Behdata = tmp_Behdata+100;
                    
                    compwas = 'Attention';
                    
                    
                case {2 3}   % xaxis rating, xaxis using AUC
                    tmp_Behdata= trial_table.Resp;
                    compwas = xratings{idatatype};
                    
                    
                    
            end
            
            
            for TARGtype = [1:3]  %use targets when objectively present, absent, or all targets,
                %% note that there are subclasses for detection type.
                switch  TARGtype
                    case 1 % present targets only.
                        restrange = find(trial_table.TargPresent==1);
                    case 2 % objectively absent only.
                        restrange = find(trial_table.TargPresent==0);
                    case 3 % all trials
                        restrange = 1:length(TargPresent);
                end
                TargClass = TargClasses{TARGtype};
                %% restrict to relevant trials, then split by quintiles.
                % Restrict range as appropriate.
                % for TargPresent, Alpha, and ratings, all used in later stages:
                tmp_Ad = AlphaMean_chanav(restrange);
                tmpBeh_adapted = tmp_Behdata(restrange);
                TargPresent_adapted = TargPresent(restrange);
                
                %>Continue by taking zscore of remaining relevnt trials, for
                %alpha amplitudes.
                tmp_Ad_z = zscore(tmp_Ad);
                %% now split alpha into quintiles.
                RespCats = splitTrialsintoBins(tmp_Ad_z,nquintiles-1);
                
                %% prepare data for storage
                
                % either calcuulate AUC for each quintile, or take the average
                % zscored subjective rating. Otherwise, calculate detection
                % accuracy.
                
                % we will look at the mean per category.
                plotOutput = zeros(2,length(RespCats));
                
                %% take z scores of all behavioural ratings also
                
                tmp_Behdata_z = zscore(tmpBeh_adapted);
                
                % store per quintile.
                for icat = 1:length(RespCats)
                    plotOutput(1,icat) = nanmean(tmp_Ad_z(RespCats{icat})); % first dim alpha
                    plotOutput(2,icat) = nanmean(tmp_Behdata_z(RespCats{icat}));
                    
                    
                    % store alpha (will overwrite each loop)
                    colname = ['Alpha_q' num2str(icat)];
                    dataTable.(colname) = plotOutput(1,icat);
                end
                
                ptitle = [compwas '_' TargClass];
                
                % if use AUC, overwrite.
                if isubjectivedata ==3 && TARGtype==3 % need all targs for AUC calcs.
                    % AUC of xratings.
                    storeAUC=zeros(2,length(RespCats));
                    % also reclass conf varibles for AUC calc.
                    
                    for icat= 1:length(RespCats)
                        AUCtrue=1;
                        Datanow = tmpBeh_adapted(RespCats{icat},:);
                        Targnow = TargPresent_adapted(RespCats{icat},:);
                        try
                            [x,y, ~, AUC]= perfcurve(Targnow, Datanow, AUCtrue);
                            
                        catch
                            % if no neg classes:
                            if sum(Targnow)==length(Targnow)
                                AUC=1; % perfect classification.
                            else
                                % only absent targs this quintile,
                                % can't calculate 
                                error('check code');
                            end
                        end
                        
                        %store AUC, rearrange
                        plotOutput(1,icat) = plotOutput(2,icat);
                        plotOutput(2,icat) = AUC;
                        
                        
                    end
                    
                    ptitle= [ptitle '_AUC'];
                    
                end
                
                %add to table (each quintile).
                for icat= 1:length(RespCats);
                    colname = [ptitle 'q' num2str(icat)];
                    dataTable.(colname) = plotOutput(2,icat);
                    
                    
                end
                
            end % 3 Targ present classes
        end % % all subjective data
       
        
        %% Now calculate and  store objective/SDT cata
        %we want to plot the detectiona accuracy. so use as input
        % data, the response (H,M,FA,CR).
        tmp_Behdata = trial_table.Outcome;
        
        for iobjdata= 1:5 %(Acc, H,M,FA,CR)
            
            %    we need to calculate detection rates for quintiles
            plotOutput = zeros(2,length(RespCats));
            
            for icat = 1:length(RespCats) % per  alpha quintile.
                
                tmptrials = RespCats{icat};
                %%
                TP = length(find(tmp_Behdata(tmptrials)==1)); %(Hits= true positive)
                Mc = length(find(tmp_Behdata(tmptrials)==2));
                TN = length(find(tmp_Behdata(tmptrials)==3)); % (CR =true negative)
                FAc = length(find(tmp_Behdata(tmptrials)==4));
                
                %% also HRr and FAr for SDT metrics.
                %HRr
                HRr= TP/(TP +Mc) ;
                %FAr
                FAr= FAc/(FAc+TN);
                %%
                if iobjdata==1 %Accuracy
                    %calculate accuracy, (H+CR)/ (H+M+CR+FA);
                    %%
                    %true class (all targ present)
                    P= sum(TargPresent(tmptrials));
                    %true negative class (all targ absent)
                    N=length(find(TargPresent(tmptrials)==0));
                    
                    datatmp= (TP+TN)/(P+N);
                    ptitle = ['accuracy'];
                    
                elseif iobjdata==2 % calculate HR
                    datatmp= HRr;
                    ptitle = ['Hitrate'];
                    
                elseif iobjdata==3 % calculate FAr
                    %compute
                    datatmp= FAr;
                    ptitle = ['FARate'];
                    
                elseif iobjdata==4 % calculate criterion.
                    datatmp= -0.5*(norminv(HRr)+ norminv(FAr));
                    
                    ptitle = ['crit'];
                    
                elseif iobjdata==5 % calculate dprime
                    
                    datatmp= norminv(HRr)-norminv(FAr);
                    
                    ptitle = ['dprime'];
                end
                
                colname = [ptitle '_q' num2str(icat)];
                
                % add to table.
                dataTable.(colname) = datatmp;
            
            end % each quintile.
        end % each datatype.
        
    %save per ppant.
    behdata_byalphaTable = dataTable;
    save(nfiles(ippant).name, 'behdata_byalphaTable', '-append');
    disp(['Fin table for ppant ' num2str(ippant) ' Exp ' num2str(idatatype) ]);
    
        end % job 1 (calc per ppant).
        
        if job.concatTable_GFX
             % load within ppant loop to save across subjs             
             load(nfiles(ippant).name, 'behdata_byalphaTable');
             
             %append to table, each row a new ppant.
             t=array2table(ippant, 'VariableNames', {'Ppant'});
             GFX_table = [GFX_table; [t, behdata_byalphaTable]];
             
             
        end
end % perppant

%save GFX table
GFX_behdata_byalphaTable= GFX_table;
f_OccSpec= ampEEG_faxis;
save(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_behdata_byalphaTable', ...
    'GFX_alphaTOPO', ...
    'GFX_OccSpec', 'f_OccSpec','-append');
end % datat type (EXP)
