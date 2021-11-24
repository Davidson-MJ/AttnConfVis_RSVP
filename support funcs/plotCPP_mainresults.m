% plotCPP_mainresults

% support script to plot the precomputed data for CPP x behaviour ratings.
% grouped and arranged for Manuscript.

% (builds off the analysis in 'plotBEH_vs_ERPfeats')

%%
% load the data for Exp 1 and Exp 2.
mydirs;
fontsize = 15;

%1 for Hit only.
%2 for miss only,
%3 for no restriction (all trials), 
TRIALrestriction = 1; %
TRIALstype = {'H only', 'M only', 'All Targets'};

CPchans = [15:17, 20:22];% for Target locked.
elocs= getelocs(2);


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
runplotGoF=1;
usequin=1; % xratings, or attention

figure(1); clf
set(gcf, 'units', 'normalized', 'position', [0, 0, 1,1])
pcount =1; % subplot counter.
RESPcolours = brewCOLOURS;
% checkwins= [.2 .55];
checkwins= [.25 .55];


job.concatGFX =0; % pre sort and save for plot
job.plotGFX =1; % plot the above.
%%
for idatatype=1:2
    xlabis = xlabsare{idatatype};
    quinsare = {xratings{idatatype}, 'Attention'};    %for labelling.
    cd(eegdir)
    
    % concat across ppants.
    if job.concatGFX           ==1
        %% for each participant,
        %load data
        cd(eegdir)
        cd( [xlabis ' alpha data']);        
        nfiles = dir([pwd filesep xlabis '_Attn_' '*.mat']);
        loaddir = pwd;
        
        [GFX_Yaxis_vsTL_erp,GFX_Xaxis_vsTL_erp] = deal(zeros(length(nfiles), 32, 426, 5));
        
        for ippant =1:length(nfiles)
            %%
            cd(loaddir);
            %plot attention and state ratings. The process is:
            load(nfiles(ippant).name);
            
            
                GFX_Xaxis_vsTL_erp(ippant,:,:,:) = squeeze(TL_bySDTandXresp_erp(1:32,:,TRIALrestriction, :));
                GFX_Yaxis_vsTL_erp(ippant,:,:,:) = squeeze(TL_bySDTandAttresp_erp(1:32,:,TRIALrestriction, :));
                
        end
        
        save(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_Xaxis_vsTL_erp','GFX_Yaxis_vsTL_erp','-append');
    
    end % concat job (if first time)
    
    %% 
    % ^^^^^^^^^^^^^^^^^^^
    
    if job.plotGFX
        cd(eegdir)
    cd( [xlabis ' alpha data']);
    load(['GFX_' xlabis '_ratingsandalpha.mat'], 'GFX_Xaxis_vsTL_erp', 'GFX_Yaxis_vsTL_erp', 'time_TLerp');
    
    % colour to use?
        

    for usequin=1:2
        
        if usequin==1 % Xratings
            
            dataIN= GFX_Xaxis_vsTL_erp;
            
            if idatatype==1
            cmapB= cbrewer('seq', 'Blues', 6);
            else
             cmapB= cbrewer('seq', 'Greens', 6);
            end
        else
            %attention
            dataIN= GFX_Yaxis_vsTL_erp;
            
            cmapB= cbrewer('seq', 'Reds', 6);
        end
        cmapB =cmapB(2:6,:); % drop the first ( too light to see).
        %
        figure(1)
        %% plot ERP, overlaying for each quintile.
        subplot(2,4,pcount)
        for iquin = 1:size(dataIN,4)
            
            tmpD=squeeze(nanmean(dataIN(:,CPchans,:,iquin),2));
            
            mP = squeeze(nanmean(tmpD,1));
            stE = CousineauSEM(tmpD);
            
            sh=shadedErrorBar(time_TLerp, mP, stE, [], 1);
            sh.mainLine.Color = cmapB(iquin,:);
            sh.mainLine.LineWidth =2;
            sh.patch.FaceColor= cmapB(iquin,:);
            sh.edge(1).Color  = cmapB(iquin,:);
            sh.edge(2).Color  = cmapB(iquin,:);
            sh.patch.FaceAlpha=0.2;
            hold on
            lg(iquin)= sh.mainLine;
            
            
            if iquin==1 % plot patch?
                Yd    =sh.edge(1).YData;
                xts= dsearchn(time_TLerp', [.35, .65]');
                xvec= [time_TLerp(xts(1)-1:xts(2)+1)];
                yvec = [0, Yd(xts(1):xts(2)), 0];
                ph = patch(xvec, yvec, [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', .2);
                pleg(1)= ph;
            end
        end
        xlabel('Time from target onset');                
        ylabel( '\muV');
        binlab = quinsare{usequin}(1:3);
%         legend(lg, {[binlab ' 1'],[binlab ' 2'], [binlab ' 3'],[binlab ' 4'],[binlab ' 5']});%, 'autoupdate', 'off')
        
        hold on
        xlim([-.1 0.7]);
        hold on;
        
        
        ylim([-1 14])
        plot([0 0], ylim, ['k:'], 'linew', 2)
        plot(xlim, [0 0], ['k-'], 'linew', .5)
        set(gca, 'fontsize', 20);
        %% Add BAR:
        
        %
        pcount=pcount+1;
        for iwin = 1:size(checkwins,1)
            figure(1);
            subplot(2,4,pcount);
            xv = [checkwins(iwin,1), checkwins(iwin,1), checkwins(iwin,2), checkwins(iwin,2)];
            yv = [-2, -1.7, -1.7, -2];
            ph=patch(xv, yv, [.9 .9 .9]);
            ph.FaceAlpha = .5;
            
            %% add bar charts
            
            
            subplot(2,4,pcount);
            tav = dsearchn(time_TLerp', [checkwins(iwin,:)]');
            tmpB = squeeze(nanmean(dataIN(:, CPchans,:,:),2));
            plotdata = squeeze(nanmean(tmpB(:, tav(1):tav(2),:),2));
            
            
            
            mBar = nanmean(plotdata,1);
            nB= length(mBar);
            
            %plot separately, to add pretty colours to bar.
            v=1:5;
            nansfill =  nan(size(plotdata));
            for ib=1:5
                % for all but this column, plot nan.
                tmp= nansfill;
                tmp(:,ib) = plotdata(:,ib);
                bh=bar(nanmean(tmp)); hold on;
                bh.FaceColor = cmapB(ib,:);
                
            end
            
            stE= CousineauSEM(plotdata);
            hold on;
            errorbar(1:5, mBar, stE, 'k', 'LineStyle', 'none');
            
            if runplotGoF == 1 % for plotting final fit:
                
%                 statsOUT=plotGoF(mBar);
                statsOUT=plotGoF(plotdata);
            end
            %         ylim([0 8])
            xlabel({[quinsare{usequin} ' quintiles'];['(low-high)']});
            ylabel({['CPP [\muV] ']});%;['[' num2str(checkwins(iwin,1)) '-' num2str(checkwins(iwin,2)) '] s']})
            set(gca, 'fontsize', 20);
            legend off
            axis tight;
            xlim([0 6])
            ylim([0 11])
            
            
        end % window
        pcount=pcount+1;
    end % axis split (quintiles to plot).
    end
end % datatype

    
    