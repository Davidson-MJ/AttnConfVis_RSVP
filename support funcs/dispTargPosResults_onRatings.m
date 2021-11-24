function     dispTargPosResults_onRatings(homedir)
%%
%across ppants. load contrast values, disp in command window the M and SD.

%called from Data_Expore.m
cd(homedir)
cd('Exp 2 and 3 mat files');

%%
Exps = {'Conf', 'PAS'};
xvals= {'Confidence', 'Visibility'};
TChip = nan(2,12,6,2); % 2 exps, 12 ppmax, 6 pos, 2 outcomes)
for iexp=1:2
    clearvars -except Exps iexp TChip xvals
    
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
            TChip(iexp,ippant,itarg-2, 1) =mean(abs(p_table.Resp(uset))); %x axis resp
            TChip(iexp,ippant,itarg-2, 2) =mean(p_table.AttResp(uset) + 100); % Att resp (adjust by 100)
            
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
Xpos = squeeze(nanmean(TChip(iexp,:,:,1),2));
Ypos = squeeze(nanmean(TChip(iexp,:,:,2),2));
if iexp==1
    Xpos = Xpos./100; % normalize to 1 on x axis
else
    Xpos = Xpos./200; % normalize to 1 on x axis
end
Ypos = Ypos./200; % normalize to 1 on y axis
subplot(2,2,iexp)
bh=bar([Xpos, Ypos]);
legend(xvals{iexp}, 'Attention')
ylim([0 1])
%%
%Plot (% errorbars)
if iexp==1
XResp =squeeze(TChip(iexp,:,:,1))./100;
else
    XResp =squeeze(TChip(iexp,:,:,1))./200;
end
YResp =squeeze(TChip(iexp,:,:,2))./ 200;
%%
subplot(2,2,2+iexp); hold on
for id=1:2
    switch id
        case 1
            plotme = XResp;
            plotc= 'b';
            
        case 2
            plotme=YResp;
            plotc='r';
    end
    npp = length(find(~isnan(plotme(:,1))));
    eh=errorbar(1:6,nanmean(plotme,1), nanstd(plotme,0,1)./sqrt(npp));
    eh.LineWidth=2;
    eh.LineStyle=':';
    eh.Color = plotc;

%     ylabel('Accuracy');
    
    xlabel('Target position')
    xlim([.5 6.5])
%     ylim([.6 .9])
    set(gca, 'Xticklabel', {'3','4',' 5', '6', '7', '8'})
    

    %% perform rm anova:
    %convert to table.
    dataT=plotme(1:npp,:);
    dataIN = dataT(:);
    subjectsIN = repmat([1:npp]', size(dataT,2),1);
    conds= repmat([1:6], size(dataT,1),1);
    condsIN=conds(:);

    rmret = rmanova(dataIN, condsIN, subjectsIN);

    rmret.ANOVAtable
end % both axes
legend(xvals{iexp}, 'Attention')
axis tight

end % both exps.
legend('Exp1', 'Exp2')
%%
end