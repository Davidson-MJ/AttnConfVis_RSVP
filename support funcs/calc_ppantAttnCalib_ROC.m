function    calc_ppantAttnCalib_ROC(homedir, datadir)


cd(homedir)
cd(datadir)


nquintiles=5; 
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
%     nppants = length(dir([pwd filesep xlabis '_Attn_participant*']));
    
    %%     dimsare = ['nppants', 'X_Attn', 'SDTm(H,M,FA,CR)', 'axis'];
    for ippant=loadnppants
        cd(datadir)
        % open p_table
        loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
        load(loadvarname);
        
        
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
        
        %We will divide the response (confidence judgements), by attention bins
        % First: the vector of confidence (negative and positive) is
        
        RespData= table2array(p_table(:,18));
        % i.e. confidence response when present/absent
        
        
        %collect and display ROC data per resp category.
        % note that we are interested in Attention data also
        AttentionData = table2array(p_table(:,19)); % on y axis, from -100 to 100.
        
        % collect data
        
        % collect data
        
        RespCats = splitTrialsintoBins(AttentionData,nquintiles-1);
        
            %%
            storeAUC=zeros(1,length(RespCats));
            storeAUCxy=[];
            
            
            for icat = 1:length(RespCats) % for each comparison
                % note that we are using Confidence data, but indexed by attention
                % category
                
                Datanow = RespData(RespCats{icat},:);
                Targnow = TargPresent(RespCats{icat},:);
                
                try
                    [x,y, ~, AUC]= perfcurve(Targnow, Datanow, PosClasslabel);
                catch
                    if sum(Targnow)==length(Targnow)
                        AUC=1;
                        x=[1,1];
                        y=[1,1];
                    else
                        error('checkcode')
                    end
                end
                    %%
                     % confirm AUC with trapz function
                     
                 
                    
                %%
                storeAUCxy(icat).xdata=x';
                storeAUCxy(icat).ydata=y';
                storeAUC(icat)=AUC;
                
                
                %sanity check;
                %         plot(x,y);
                
                
                
            end
          
                    RespxAttn_Multisplit_AUC = storeAUC;
                   
           
        
        %save per ppant.
        
            save(loadvarname, 'RespxAttn_mediansplit_AUC', ...
                'RespxAttn_mediansplit_AUCxy',...
                'RespxAttn_Qtlsplit_AUC',...
                'RespxAttn_Qtlsplit_AUCxy',...
                'RespxAttn_Multisplit_AUC',...
                'Qp', 'attnMedian', '-append');
        
        %%
    end % per ppant
end % per experiment.
%



%% calculate d prime
%
%           %calculate frequencies per SDTm:
%           H2= length(find(SDTindex(1,:)));
%           M2= length(find(SDTindex(2,:)));
%           FA2=length(find(SDTindex(3,:)));
%           CR2=length(find(SDTindex(4,:)));
%
%
%           if sum([H2,M2, FA2, CR2]) ~=936
%               error('incorrect trial count!');
%           end
%
%           h = sum(H2) / (sum(H2)+sum(M2));
%           fA= sum(FA2) / (sum(FA2)+sum(CR2));
%
%
%           % d prime = z(h)-z(fA)
%           dprime = norminv(h)-norminv(fA);
%
%           % c = -0.5*[z(h)+z(fA)]
%           dcrit = -0.5*(norminv(h)+ norminv(fA));