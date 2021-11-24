function extractMouseClicksperppantandSave(homedir, datadir)
% steps through the matlab files per participant, extracting the scatter
% plot representation of their mouse clicks across all trials. Sorted by
% 4x SDT measures (SDTm; Hits, FA, Miss, CorR).

%saves:
% dimsare: SDTm - Hits, FA, Miss, CorR
% XY_RESP_data: (ntrials x XYdim) used for scatterplot
% SDTindex: vector containing logical index assigning SDTm to each point
% above.
% SDTmetrics_on_Yaxis: Cumulative total datapoints X attention (collapsing
% onto y axis)
% SDTmetrics_on_Xaxis: Cumulative total datapoints X Confidence or PAS 
% (depending on experiment, collapsing onto X axis).

%

%%%%%% MD June 2019

cd(homedir)
cd(datadir)

for iExp=1:2
    
    if iExp==1
        Axis_X = -100:100; % confidence measure.
        xlabis='Conf';
        loadnppants = [1:12];
    else
        Axis_X = 0:200; % PAS measure.
        xlabis='PAS';
        loadnppants = [2:10];
        
    end
        
    
%     dimsare = ['nppants', 'X_Attn', 'SDTm(H,M,FA,CR)', 'axis'];
    for ippant=loadnppants
        
          % open p_table
          loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
          load(loadvarname);
             
        
        %%  Columns in data set:
        %         DV               Var1
        %     _______________________    ____
        %
        %     'DateTime'                   1
        %     'Experiment'                 2
        %     'Refreshrate'                3
        %     'Session'                    4
        %     'Block'                      5
        %     'Trial'                      6
        %     'BlockCond'                  7
        %     'Condition'                  8
        %     'TargPresent'                9
        %     'TargOrientation'           10
        %     'TargChip'                  11
        %     'RespTrigger'               12
        %     'BinResp'                   13
        %     'RespID'                    14
        %     'Outcome'                   15
        %     'OutcomeID'                 16
        %     'RT'                        17
        %     'Resp'                      18
        %     'AttResp'                   19
        %     'StimStartError'            20
        %     'StimEndError'              21
        %     'DesiredDuration'           22
        %     'Duration1'                 23
        %     'Duration2'                 24
        %     'EstFramesMissed'           25
        %     'FramesMissed'              26
        %     'GreyFramesMissed'          27
        %     'FixTimeExceeded'           28
        %     'MakeTexturesTimeTaken'     29
        %     'Quest'                     30
        %     'TargContrast'              31
        %     'Q1Threshold'               32
        %     'Q1Thr_SD'                  33
        
        
        %% Columns of interest are the confidence response (±100), and attention 0-200.
        % these are located in columns 18,19 (RTs in 17)
        %SDT metrics: Hits, misses, FA, CORej . which were which rows?
        
        XY_RESP_data = table2array(p_table(:,18:19));
        
        XY_RT_data = table2array(p_table(:,17));
        
        SDTindex= zeros(4,size(XY_RESP_data,1));
        
        dimsare={'Hit', 'Miss', 'FA', 'CorR'};
        %collect which rows satisfy each condition:
        SDTindex(1,:)= ismember(p_table.OutcomeID,{'HIT'});
        SDTindex(2,:)= ismember(p_table.OutcomeID,{'MIS'});
        SDTindex(3,:)= ismember(p_table.OutcomeID,{'FAA'});
        SDTindex(4,:)= ismember(p_table.OutcomeID,{'COR'});
        
       % append to saved participant data.
        save([loadvarname], 'SDTindex', 'XY_RESP_data', 'XY_RT_data', '-append');
        
        %% now calc density across axes (for plotting distributions).
        %first create vector like cdf
        %confidence is on the x-axis, 2nd column, so collapse across attention.
        %we can sort the values. then plot by bin.
        
        for iDIM=1:2          
            %note that the range changes between experiments (0:200,
            %-100:100)           
            %default
            this_X=-100:100;           
            
            if iDIM==1 && iExp==2 % x-axis for  PAS scale.
                this_X=0:200;
            end
            
            
            outgoing_tmp=zeros(4,201);
            for iSDT=1:4
                %use relevant rows to sort data.
                plot_tmp = sort(XY_RESP_data(logical(SDTindex(iSDT,:)),iDIM),'ascend');
                
                % record instances per bin. 
                plot_dens=zeros(1,length(this_X));                
                for i = 1:length(this_X)
                    xpoint=this_X(i);
                    plot_dens(i) = sum(plot_tmp==this_X(i));
                end
                
                %store for saving per ppant:
               outgoing_tmp(iSDT,:)=plot_dens;
            end
            
              switch iDIM
                  case 1
                      SDTmetrics_on_Xaxis= outgoing_tmp;
                  case 2
                      SDTmetrics_on_Yaxis= outgoing_tmp;
              end
                      
        end
        
        %note that we can also compute the correlation value between x and
        %y axes
        if iExp==1 % attention x confidence
            %then x axis is -100 to 100. we can take the absolute value of
            %x axis
            corrXY=[abs(XY_RESP_data(:,1)), XY_RESP_data(:,2)];
        else
            % x values range from 0 :200, 
            % we are only interested in the non zero numbers of the xaxis.
                idx= XY_RESP_data(:,1)>0;
                corrXY = [(XY_RESP_data(idx,1)), (XY_RESP_data(idx,2))];
        end
        
        [rho, corr_pval]= corr(corrXY, 'type', 'Spearman', 'tail', 'both');
        
        
        save(loadvarname, 'SDTmetrics_on_Xaxis', 'SDTmetrics_on_Yaxis', 'dimsare', 'rho', 'corr_pval', 'corrXY','-append')
    end %nppants
    
end %iExp type.
