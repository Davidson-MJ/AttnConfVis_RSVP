function combineAttnCalibacrossppantandSave(homedir,datadir)
% combines AUC data across participants
% MD June 2019


cd(homedir )
cd(datadir)
xlabsare={'Conf', 'PAS'};
%%
for iExp=1:2  % for both experiments
    
    xlabis = xlabsare{iExp};
    
    cd(datadir)
    %how many data types/ppants this folder.
    if iExp==1
        loadnppants = [1:12];
    else
        loadnppants = [2:10];
    end
    
    
    nppants = length(loadnppants);
    
    [all_AUCbymedian,...
        all_AUCbymedian_PRver,...
        all_AUCbymedian_nback,...
        all_Accuracybymedian] = deal(zeros(nppants,2));
    
    [all_AUCbyquantile,...
        all_AUCbyquantile_PRver,...
        all_AUCbyquantile_nback,...
        all_Accuracybyquantile]= deal(zeros(nppants, 3));
    
    [all_AUCbymultiforpred,...        
        all_Accuracybymultiforpred, ...
        all_HRbymulti,all_FArbymulti,all_Dprimebymulti,...
        all_Critbymulti]= deal(zeros(nppants, 5));
    
    
    [all_AUCbyquantile_PRver,all_AUCbyquantile_nback]= deal(zeros(nppants, 4));
    
    
    icounter=1;
    for ippant=loadnppants
        cd(datadir)
        % open p_table
        loadvarname= [ xlabis '_Attn_participant' num2str(ippant) '.mat'];
        load(loadvarname, 'RespxAttn*', 'store*', 'nshift');
        
        %add objective accuracy per participant:
        
%         all_Accuracybymedian(icounter,:) = RespxAttn_mediansplit_Accuracy;
%         all_Accuracybyquantile(icounter,:) = RespxAttn_Qtlsplit_Accuracy;        
        all_Accuracybymultiforpred(icounter,:) = RespxAttn_Multisplit_Accuracy;
        all_HRbymulti(icounter,:)= storeDetection;
        all_FArbymulti(icounter,:)=storeFA;
        all_Dprimebymulti(icounter,:)= storeDprime;
        all_Critbymulti(icounter,:)=storeCrit;
        
    % add AUC measure (confidence calibration)
%         all_AUCbymedian(icounter, :)=RespxAttn_mediansplit_AUC;
%         all_AUCbyquantile(icounter, :)=RespxAttn_Qtlsplit_AUC;
        all_AUCbymultiforpred(icounter,:) = RespxAttn_Multisplit_AUC;
        
        %add precision recall versions (for plot comparisons).        
        all_AUCbymedian_PRver(icounter, :)=RespxAttn_mediansplit_AUC_PRver;
        all_AUCbyquantile_PRver(icounter,:)=RespxAttn_Qtlsplit_AUC_PRver;
        
        
        %also save n-back versions
         
        all_AUCbymedian_nback(icounter, :)=RespxAttn_mediansplit_AUC_nback;
        all_AUCbyquantile_nback(icounter, :)=RespxAttn_Qtlsplit_AUC_nback;
        
        
        
        %add precision recall versions (for plot comparisons).
%         all_AUCbymedian_PRver_nback(ippant, :)=RespxAttn_mediansplit_AUC_PRver_nback;
%         all_AUCbyquantile_PRver_nback(ippant,:)=RespxAttn_Qtlsplit_AUC_PRver_nback;
        
        

       icounter=icounter+1; 
    end
%     save appropriately.
            savename= [xlabis '_Attn_Allpp'];
            save(savename, ...
                'all_Accuracybymultiforpred',... % accuracy by attention quintiles
                'all_AUCbymultiforpred',... % AUC
                'all_HRbymulti',... % % HRate
                'all_FArbymulti',...% FArate
                'all_Dprimebymulti',... % dprime
                'all_Critbymulti', '-append') % criterion.
end
end