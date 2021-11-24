% function  calc_ppantAcc_ROC(datadir, homedir);
cd(homedir)
cd(datadir)
dbstop if error

storeAccuracy = nan(2,12);
storeDetection = nan(2,12);
storeFA = nan(2,12);
storeDprime = nan(2,12); % dprime (sensitivity), for both experiments.
storeCrit = nan(2,12);

storeAUCxy = [];
storeAUC_Xax = nan(2,12,3); % last dim, include tota, present and absent case.
storeAUC_Yax = nan(2,12,3); % last dim, include tota, present and absent case.
storeConf = nan(2,12,4); % total absolute confidence in each type (SDT)
storeAttn = nan(2,12,4); % total absolute confidence in each type (SDT)
%colormap
cmapB= cbrewer('qual', 'Paired', 12);    
    col1 = cmapB(8,:);
    col2 = cmapB(10,:);
    %%
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
        
        load(loadvarname, 'p_table')
        % we want the accuracy and AUC for all detection (present absent).
        
        %Accuracy first.
        
        %in each case the 'true' class labels, are target present vs absent:
        %true labels (present vs absent) in table:
        
        TargPresent = table2array(p_table(:,9)); % 9th dimension in participant table.
        PosClasslabel = 1; % what is the correct case in this vector
        
        % collect data                       
                %CALCULATING ACCURACY instead of AUC:
                 %Accuracy = (TP+TN)/(P+N)                 
                  % true positive (Hits)                 
                  TP=sum((p_table.Outcome==1)); 
                  % true negative (CR)
                  TN=sum((p_table.Outcome==3));
                  %true class (all targ present)
                  P= sum(TargPresent);
                  %true negative class (all targ absent) 
                  N=length(find(TargPresent==0));
                  
                  AccuracyTmp= (TP+TN)/(P+N);

                 storeAccuracy(iExp, ippant)= AccuracyTmp;
              
                 %% detection rate = H when present                 
                 storeDetection(iExp, ippant) = TP/P;
                 
                 FA=sum((p_table.Outcome==4));
                 %% store FA rate
                 storeFA(iExp,ippant) = FA/N; 
                 
                 %% store d' c:
                 HR = TP/P;
                 FAR = FA/N;
                 
                 % d prime = z(h)-z(fA)
                 dp = norminv(HR)-norminv(FAR);
                 
                 % c = -0.5*[z(h)+z(fA)]
                 c = -0.5*(norminv(HR)+ norminv(FAR));
                 
                 storeDprime(iExp,ippant)=dp;
                 storeCrit(iExp,ippant) = c;
                 
                 RespData= table2array(p_table(:,18));
                 %% store AUC for all trials.
                [x,y, ~, AUC]= perfcurve(TargPresent, RespData, 1);
                storeAUC_Xax(iExp, ippant, 3) = AUC;
                 
                %also store AUC based on x axis descriptors.
                RespData= table2array(p_table(:,19));
                [~,~, ~, AUC]= perfcurve(TargPresent, RespData, 1);
                
                storeAUC_Yax(iExp,ippant,3) = AUC;
                
                usetrials=[];
                 
                 %% also grab AUC for present vs absent.
                 %x axis responses will be used:
                 usetrials(1).tr= find(RespData>2); % perceived present
                 usetrials(2).tr= find(RespData<=1);% perceived absent
                 
                % for both classes:
                classlabels = [1,0]; 
                cols= [col1;col2]; 
                %also store for subset of perceived present vs absent.
                
                for iPresorAbs =1:2
                   
                Datanow = RespData(usetrials(iPresorAbs).tr); 
                
                 if iPresorAbs==2
                     Datanow = Datanow*-1; % reward lower values.
                 end
                Targnow = TargPresent(usetrials(iPresorAbs).tr);
                
                [x,y, ~, AUC]= perfcurve(Targnow, Datanow, classlabels(iPresorAbs));
                    
                %%
                storeAUCxy(iExp, ippant,iPresorAbs).xdata=x';
                storeAUCxy(iExp, ippant,iPresorAbs).ydata=y';
                storeAUC(iExp, ippant, iPresorAbs)=AUC;
%                  hold on
%                 plot(x,y, 'color', cols(iPresorAbs,:))
                end
                
                
                allconf = p_table.Resp;
                allAttn = p_table.AttResp;
                %store total confidence per type
                for iSDT=1:4
                    outc = find(p_table.Outcome==iSDT);
                storeConf(iExp, ippant, iSDT) = nanmean(allconf(outc));
                storeAttn (iExp, ippant, iSDT)= nanmean(allAttn(outc));
                end
    
    end % ppant
  
      
      
    end % per exp
    %% plot results
    %First plot boxplots . with scattered points for accuracy
    
    plotBOXES_AccandROC;
    %% now plot the individual ROC curves.
    plotPPANT_ROCcurves;
    
    % end %function
    %