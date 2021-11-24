% plotBOXES_AccandROC

% called from: calc_ppantAcc_ROC.m
clf
fontsize = 10;
SCcol = {'b', 'r'}; % colours.

RESPColours = brewCOLOURS;


figsize = [ 0 0 1 .4];

%have already done the stats:


figure(1);    clf;
set(gcf, 'Units', 'normalized','Position',figsize, 'color', 'w');
%%
for iBOXtype =3;%[1,2,3]%1:3
    switch iBOXtype
        case 1
            alldataIN = storeAccuracy;
            ylabmark = 'Accuracy';
%           
            params.scatCol = [squeeze(RespColours(1,1,:))'; squeeze(RespColours(2,1,:))'];
            
            subplot(1,4,1)
            ylimsare=[.4 1.1];
        case 2
            alldataIN = squeeze(storeAUC_Xax(:,:,3));
            alldataINt = squeeze(storeAUC_Yax(:,:,3));
            alldataIN = cat(1, alldataIN, alldataINt);
            ylabmark = 'Type-2 (AUROC2)';
            
            %max colours, conf, vis, attention:
            Cmax = squeeze(RESPColours(1,1,:));
            Vmax = squeeze(RESPColours(2,1,:));
            Amax = squeeze(RESPColours(3,1,:));
%             SCcol = {Cmax, Vmax, Amax, Amax}; % colours.
            params.scatCol = [Cmax'; Vmax'; Amax'; Amax'];
            ylimsare=[.4 1.1];
            lgP= {'Confidence rating', 'Visibility rating', 'Attention rating'};
           subplot(1,4,3:4);
        case 3
            alldataIN = storeDprime;
            ylabmark = 'd''';
            subplot(1,4,2)
            ylimsare = [.51 4.5];
            params.scatCol = [squeeze(RespColours(1,1,:))'; squeeze(RespColours(2,1,:))'];
        case 4
            
            alldataIN = storeCrit;
            ylabmark = 'crit.';
        case 5
            
            alldataIN = storeDetection;
            ylabmark = 'HR.';
        case 6
            
            alldataIN = storeFA;
            ylabmark = 'FArate.';
      
    end
    %%
    params.cols = [0,0,0]; %cols for boxes.
    
    params.plotScatter=1;
    params.scatSize = 40;
    %params.scatCols% defined above.
            
    box_and_scatter(alldataIN', params);
  
 
 set(gca, 'Xticklabels', {'Exp. 1' 'Exp. 2'})

ylim([ylimsare]);
hold on;
plot(xlim, [0.5 .5], ['k--'])
set(gca, 'fontsize', fontsize*1.5)
ylabel(ylabmark)

if iBOXtype==2
% lg=legend(lgP, 'Location', 'NorthEastOutside')
end
%% display result summary in command window (tailored for type-2 sensitivity)
icontr = [1,2; 3,4; 1,3; 2,4]; % condition comparisons (refers to pos in box plot).
exps  =[1,2; 1,2; 1,1; 2,2]; % which experiment the data came from
ttesttype = [1,1,2,2];% independent or paired samples.
sigbar = [.99, .65, 1.03, .93];

for icomp = 1:3%size(icontr,1)
    condsare = icontr(icomp,:);
    expsare = exps(icomp,:);
    ttesttypetmp = ttesttype(icomp); 
disp(['-----------------------------------------------------------'])
disp([ ylabmark ' comparison ' num2str(icomp)])
disp([ 'Mean exp ' num2str(expsare(1)) ':' sprintf('%.2f',nanmean(alldataIN(condsare(1),:))) ', SD  = ' sprintf('%.2f',nanstd(alldataIN(condsare(1),:),0,2))])
disp(['Mean exp ' num2str(expsare(2)) ':' sprintf('%.2f',nanmean(alldataIN(condsare(2),:))) ', SD  = ' sprintf('%.2f',nanstd(alldataIN(condsare(2),:),0,2))])
%%
if ttesttypetmp==1
[h,p,CI,stat] = ttest2(alldataIN(condsare(1),:), alldataIN(condsare(2),:));
%calculate effect size : Welch test:
% d = Ma - Mb / (sqrt(varA + varB)/2)
Ma= nanmean(alldataIN(condsare(1),:)); Mb=nanmean(alldataIN(condsare(2),:));
varA = nanvar(alldataIN(condsare(1),:)); varB = nanvar(alldataIN(condsare(2),:));
d = (Ma - Mb) / sqrt((varA + varB/2));
 
disp(['ttest (indep): t(' num2str(stat.df) ')=' num2str(stat.tstat) ', p=' num2str(p) ', d= ' num2str(d)])
elseif ttesttypetmp==2
    [h,p,CI,stat] = ttest(alldataIN(condsare(1),:), alldataIN(condsare(2),:));
%calculate effect size: 
% d= mean(diff) / std(diff)
diffs = alldataIN(condsare(1),:) - alldataIN(condsare(2),:);
d= nanmean(diffs)/nanstd(diffs);
disp(['ttest (paired): t(' num2str(stat.df) ')=' num2str(stat.tstat) ', p=' num2str(p) ', d= ' num2str(d)])
end
    

if p<.001
    usem = '***';
elseif p<.01
    usem = '**';
elseif p<.05
    usem = '*';
else
    usem = 'ns';
end

if iBOXtype==3; % dprime
    sigbar(icomp) = 3.9; % place ns at right height.
end
tc=text(mean(condsare), sigbar(icomp), usem,'fontsize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

plot([condsare(1) condsare(2)], [sigbar(icomp) sigbar(icomp)], '-k', 'linew', 1); shg

end
%%

% tc=text(1.5, .97, usem,'fontsize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% plot([1 2], [.95 .95], '-k', 'linew', 2); shg


end
%%
%% plot bar for attn/conf as well.
% 
% cmapB= cbrewer('qual', 'Paired', 12);
% col1 = cmapB(8,:);
% col2 = cmapB(10,:);
% 
% %create data!
% clf
% for it = 1:2
%     switch it
%         case 1
% useD = storeConf;
%         case 2
% useD = storeAttn;
%     end
% types = {'Confidence', 'Attention'};
% figure(2); subplot(1,2,it);
% % subplot(1,2,iBOXtype)
%     set(gcf, 'units', 'normalized', 'position', figsize, 'color', 'w')
%     
%     barD = squeeze(useD(2,:,:));
%     %reshape to group like types together;
% %     barDn(:,1)= barD(:,1); % H
% %     barDn(:,2)= barD(:,3); % CR
% %     barDn(:,3)= barD(:,2); % M
% %     barDn(:,4)= barD(:,4); % FA
%     
%     barM = abs(squeeze(nanmean(barD,1)));
%     stE = CousineauSEM(barD);
%     
%     bh=bar(barM);
%     hold on;
%     errorbar(1:4, barM, stE, 'k','Linestyle', 'none');
% %     set(gca, 'XtickLabel', {'Hit', 'FA', 'CR', 'Miss'}, 'fontsize', 15);
%     set(gca, 'XtickLabel', {'Hit', 'M', 'CR', 'FA'}, 'fontsize', 15);
%     bh.FaceColor = 'flat';
%     bh.CData(1,:) = col1;
%     bh.CData(2,:) = col2;
%     bh.CData(3,:) = col2;
%     bh.CData(4,:) = col1;
% 
%     ylabel(types{it});
% end
%     
