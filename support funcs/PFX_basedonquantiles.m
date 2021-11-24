function PFX_basedonquantiles(homedir)

%concat EEG on participant level, across all participants for saving.
%%
dbstop if error
uselowpfilt=0; % change to 1 for using lowpfiltered data.
plotTercileorQuartile = 1; % 1 = tercile split based on subj range. 2 = hard quartile boundary.

%which type of data to plot? adjusts xlimits, and print name when outgoing.
useTargetLockedorwholetrial = 2;% 1 or 2

% % % type of data to plot:
figdir= '/Users/mdavidson/Desktop/Frontiers Project/Documents/Figures';
for dset =1:2
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
    
    % orient and load group data
    cd(homedir);
    cd('Exp 2 and 3 processed EEG');
    
    
    
    if useTargetLockedorwholetrial == 1
        finddir = [ xlabis ' VAN data'];
        typeis  = 'target-locked ERPs';
        usexlim=[-.1 1];
    else
        finddir= [ xlabis ' whole epoch data'];
        typeis  = 'whole epoch ERPs';
        usexlim=[-.5 3];
    end
    if uselowpfilt ==1;
        finddir= [finddir ' lowpfilt'];
    end

    
    if plotTercileorQuartile ==1
        cd([finddir ' tercilesplit']);
    else
        cd([ finddir ' quartilesplit'])
    end
    
    eegdir= pwd;
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
        
        %plot at which channel?, replicating JMAC frontiers: Cz, CPZ,PZ,POZ,OZ
        plotchannel = [16,21,26,29,31];
        %set colours:
        [cmaps] = cbrewer('qual', 'Set1', 5);
        %%
        for iAXIStocompare = 1:2
            
            figure(iAXIStocompare); clf; set(gcf, 'units', 'normalized', 'position', [0 .4 .4 .45])
            
            if iAXIStocompare==1
                tis= 'Attention';
                
                % subselect data.
                DATAtoplot  = outgoingEEG_Attnsplit;
                
            elseif dset==1
                tis = 'Confidence';
                % subselect data.
                DATAtoplot  = outgoingEEG_XAXISsplit;
            elseif  dset==2
                tis = 'Visibility';
                DATAtoplot  = outgoingEEG_XAXISsplit;
            end
            %%
            
            for idata = 1%:4
                
                switch idata
                    case 1
                        subc = 'HITs';
                        
                    case 2
                        subc = 'MISS';
                        
                    case 3
                        subc = 'FAs';
                        
                    case 4
                        subc = 'CRs';
                        
                end
                %%
                datatoplot = squeeze(DATAtoplot(idata,:,:,:));
                
                clf
                
                for iquantile = 1:size(datatoplot,1)
                    % plot(time_new, squeeze(mean(datatoplot(:,plotchannel,:),1)), 'color', cmaps(idata,:), 'linew', 3);
                    dis = squeeze(mean(datatoplot(iquantile,plotchannel,:),2));
                    
                    col = cmaps(iquantile,:);
                    sh= plot(usetime, dis);
                    
                    sh.Color = col ;
                    sh.LineWidth = 1;%
                    
                    
                    shleg(iquantile) = sh;
                     hold on
                end
                
                title([subc ' split by ' tis ' at ' elocs32(plotchannel).labels ' (' xlabis ')'])
                try
                    if plotTercileorQuartile==1
                        lg=legend([shleg], {'Lowest', 'Medium', 'Highest'}, 'autoupdate', 'off');
                    else
                        lg=legend([shleg], {'Q1', 'Q2', 'Q3', 'Q4'}, 'autoupdate', 'off');
                    end
                catch
                    legend off
                end
                set(gca, 'fontsize', 25);
                xlim([usexlim])
                ylim([-15 25])
                set(gcf, 'color' ,'w')
                hold on; plot([0 0], ylim, ['k-'])
                hold on; plot(xlim, [0 0 ], ['k-'])
                ylabel('uV');
                xlabel('Time from target present [s]')
                
                
                %% print to figure folder
                cd(figdir)
                cd('Target locked ERPs');
                cd('PFX_splitbycategory')
            
                print('-dpng', ['participant count ' num2str(ippant) ', ' typeis ', '  subc ' split by ' tis ' ' xlabis]);
            end % each data type (Hit,M,FA,CR)
        end
            
            
        end
        % %%
        
        
    end
