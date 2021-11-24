function     dispTargPosResults_onAccuracy(homedir)
%%
%across ppants. load contrast values, disp in command window the M and SD.

%called from Data_Expore.m
cd(homedir)
cd('Exp 2 and 3 mat files');

%%
Exps = {'Conf', 'PAS'};
TChip = nan(2,12,6,2); % 2 exps, 12 ppmax, 6 pos, 2 outcomes)
for iexp=1:2
    clearvars -except Exps iexp TChip
    
    allpp = dir([pwd filesep Exps{iexp} '*_participant*']);
    %%
    for ippant = 1:length(allpp)
        load(allpp(ippant).name , 'p_table');
        
        tableTChip(1,1) = table(ippant);
        % for each targ position, calculate the %correct.
        prestrials = find(p_table.TargPresent==1);
        
        for itarg = 3:8
            
            %position was
            reltrials = find(p_table.TargChip==itarg);
            
            %targ present (not absent trial).
            uset = intersect(prestrials, reltrials);
            
            resps = p_table.Outcome(uset);
            TChip(iexp,ippant,itarg-2, 1) =length(find(resps==1));% Hits.
            TChip(iexp,ippant,itarg-2, 2) =length(find(resps==2));% Misses.
            
        end
        
    end

%
disp(['--------------------------'])
disp(['Mean Targ Chip results:'])
%H by pos.
figure(1); clf; set(gcf, 'color','w');
%plot H, M,Acc by pos.
cols = {'b', 'r'};
end
%%
for iexp=1:2
Hpos = squeeze(nanmean(TChip(iexp,:,:,1),2));
Mpos = squeeze(nanmean(TChip(iexp,:,:,2),2));

subplot(2,2,iexp)
bh=bar([Hpos, Mpos]);

%accuracy (% H)
Ht =squeeze(TChip(iexp,:,:,1));
Mt =squeeze(TChip(iexp,:,:,2));
Gt = Ht+Mt;

Acc = Ht./Gt;
subplot(2,2,3); hold on
npp = length(find(~isnan(Ht(:,1))));
eh=errorbar(1:6,nanmean(Acc,1), nanstd(Acc,0,1)./sqrt(npp));
eh.LineWidth=2;
eh.LineStyle=':';
eh.Color = cols{iexp};
ylabel('Accuracy');
xlabel('Target position')
xlim([.5 6.5])
ylim([.6 .9])
set(gca, 'Xticklabel', {'3','4',' 5', '6', '7', '8'})


%% perform rm anova:
%convert to table.
dataT=Acc(1:npp,:);
dataIN = dataT(:);
subjectsIN = repmat([1:npp]', size(dataT,2),1);
conds= repmat([1:6], size(dataT,1),1);
condsIN=conds(:);

rmret = rmanova(dataIN, condsIN, subjectsIN);

rmret.ANOVAtable
end
legend('Exp1', 'Exp2')
%%
end