function  calc_ppantAttnCalib_objectiveAccuracy(homedir,datadir)
cd(homedir)
cd(datadir)
dbstop if error

nquintiles = 5;

for iExp=1:2
    
    if iExp==1
        
        xlabis='Conf';
         loadnppants = [1:12];
    else
        
        xlabis='PAS';
        loadnppants = [2:10];
    end
    cd(datadir)
    %how many data types/ppants this folder.
    nppants = length(dir([pwd filesep xlabis '_Attn_participant*']));
    
    %%     dimsare = ['nppants', 'X_Attn', 'SDTm(H,M,FA,CR)', 'axis'];
    for ippant=loadnppants
        cd(datadir)
        % open p_table
        loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
        load(loadvarname, 'p_table', 'SDTindex');
        
        
        %
        %% type 1 AUC analyses.
        
        % to perform this analysis, we need to compare a response category with a
        % correct class category.
        %
        %Cycle through different response categories, defined by attention
        %subset.
        
        
        %in each case the 'true' class labels, are target present vs absent:
        %true labels (present vs absent) in table:
        
        TargPresent = table2array(p_table(:,9)); % 9th dimension in participant table.
        PosClasslabel = 1; % what is the correct case in this vector
        
        % we will use the H, M, FA,COR for these data below.
        
        
        %collect and display ROC data per resp category.
        % note that we are interested in Attention data also
        AttentionData = table2array(p_table(:,19)); % the YAXIS response (in clicks)
        
        % collect data
        
        RespCats = splitTrialsintoBins(AttentionData,nquintiles-1);
        
        
        %%
        
        [storeAccuracy,storeDetection,...
            storeFA,storeDprime,...            
            storeCrit]= deal(zeros(1,length(RespCats)));
            
            
            for icat = 1:length(RespCats) % for each comparison 
                % note that we are using Xaxis data indexed by y-axis  
                % (attention) bins

                
                %CALCULATING ACCURACY instead of AUC:
                 %Accuracy = (TP+TN)/(P+N)
                 
                 %
                  %collect SDT data for all trials in these attention bins:  
                  SDTNow= SDTindex(:,RespCats{icat});
                  % collect true-label class for same trials
                  Targnow = TargPresent(RespCats{icat},:);
                  
                  %SDT rows are H,M,FA,CR:
                  % true positive (Hits)                 
                  TP=sum(SDTNow(1,:)); 
                  % true negative (CR)
                  TN= sum(SDTNow(4,:));
                  %true class (all targ present)
                  P= sum(Targnow);
                  %true negative class (all targ absent) 
                  N=length(find(Targnow==0));                  
                  
                  FA=sum(SDTNow(3,:));
                  
                  HR = TP/P;
                  FAr = FA/N;
                 %% store various metrics.
                 storeAccuracy(icat)= (TP+TN)/(P+N);
                 
                 storeDetection(icat) = HR;
                 storeFA(icat) = FAr;
                 
                 %avoid floor/ceiling effects.
                 if FAr==0
                     FAr=.001;
                 end
                 if HR==0
                     HR=.001;
                 end
                     
                 storeDprime(icat) =   norminv(HR)-norminv(FAr);           
                 if isinf(storeDprime(icat))
                     break
                 end
                 storeCrit(icat) = -0.5*(norminv(HR)+ norminv(FAr));
                             
            end
            
            %rename to save.
               RespxAttn_Multisplit_Accuracy = storeAccuracy;                                           
               %
               
        
        %save per ppant.
        
            save(loadvarname, 'RespxAttn_Multisplit_Accuracy',...
                'storeDetection', 'storeFA', 'storeDprime', 'storeCrit','-append');
        
        %%
        disp(['fin ippant ' num2str(ippant)])
    end % per ppant
end % per experiment.
%