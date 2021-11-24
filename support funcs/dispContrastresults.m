function     dispContrastresults(homedir)
%%
%across ppants. load contrast values, disp in command window the M and SD.

%called from Data_Expore.m
cd(homedir)
cd('Exp 2 and 3 mat files');

%%
Exps = {'Conf', 'PAS'};
allCV = nan(2,12);
for iexp=1:2
    clearvars -except allCV  Exps iexp
    
    allpp = dir([pwd filesep Exps{iexp} '*_participant*']);
    
    for ippant = 1:length(allpp)
        load(allpp(ippant).name , 'p_table');
        
        CVp = unique(p_table.TargContrast);
        
        allCV(iexp, ippant) = mean(CVp);
    end
end
%%
disp(['--------------------------'])
disp(['--------------------------'])
disp(['Mean Targ Congtast values:'])
disp(['Exp 1 = ' sprintf('%.2f', nanmean(allCV(1,:))), ' (sd=' sprintf('%.2f', nanstd(allCV(1,:),0,2)) ')'])
disp(['Exp 2 = ' sprintf('%.2f', nanmean(allCV(2,:))), ' (sd=' sprintf('%.2f', nanstd(allCV(2,:),0,2)) ')'])
%%
[h,p,CI,stats]=ttest2(allCV(1,:), allCV(2,:))
end