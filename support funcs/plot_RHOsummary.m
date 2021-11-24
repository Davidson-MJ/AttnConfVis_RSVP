

% script to plot pretty Rho summary based on input
% MS fig 2

figure(8); clf;
set(gcf, 'units', 'normalized','position',[0 0 .5 1]);
rhoPLOTS=nan(nppants, 2);
pvals=nan(1,3);
%also plot rho across ppants.
%for present, absent, combined.
pvals=[];
for iTargType=1:2
    for ippant =1:size(alltoScatter,1)
        switch iTargType
            case 1% all target objectively present .
                tmp1 = find(allSDTindex(ippant,1,:));% Hit
                tmp2= find(allSDTindex(ippant,3,:));% FA
                userows = [tmp1; tmp2];
                xcase = 'Perceived present';
            case 2
                tmp1 = find(allSDTindex(ippant,2,:));% M
                tmp2= find(allSDTindex(ippant,4,:));% CR
                userows = [tmp1; tmp2];
                xcase = 'Perceived absent';
            case 3
                userows = 1:size(allSDTindex,3);
                xcase = 'All Targets';
        end
        Xv = squeeze(alltoScatter(ippant, userows,1));
        Yv = squeeze(alltoScatter(ippant, userows,2));
        
        
        storer = corr(Xv',Yv', 'type', 'Spearman');
        rhoPLOTS(ippant, iTargType)=storer;
        
    end
    
    xlabsare{iTargType} = xcase;
    try [~,pvals(iTargType)] = ttest(rhoPLOTS(:,iTargType));
    catch
    end
    %note for one sample ttests, d = (m-u)/s
end
%

rhoPLOTS(:,2)=rhoPLOTS(:,2).*-1;
barM = squeeze(mean(rhoPLOTS,1));

% RHO plots: (use different colour).
cmapB= cbrewer('seq', 'Blues', 4);
%plot parameters (H,M,FA,CR):

col1 = cmapB(4,:);
col2 = cmapB(2,:);
useCols= [col1;col2];

% hold on
% params.cols =[col1; col2]; % present, absent.
% params.plotScatter=1;
% h= box_and_scatter(rhoPLOTS, params);
% 
% rd=cell(2,1);
% rd{1,:}=rhoPLOTS(:,1)';
% rd{2,:}=rhoPLOTS(:,2)';
% rm_raincloud(rd, params.cols(1,:));
% % 

bh=bar(barM); hold on
bh.FaceColor= col1;
b2=bar([nan, barM(2), nan]);
b2.FaceColor = col2;
stD= std(rhoPLOTS,1);
eh=errorbar(1:size(barM,2), barM, stD,'k', 'LineStyle', 'none', 'LineWidth', 2);
% overlay ind data points:
  hold on;
    
    for id=1:size(rhoPLOTS,2)
        useD = rhoPLOTS(:,id);
        jitter = zscore(rand(length(useD), 1))/40;
        sc=scatter(id+jitter, useD);
        sc.LineWidth = 2;
        sc.SizeData = 100;
        sc.MarkerFaceColor = useCols(id,:);%[.9 .9 .9];
        sc.MarkerEdgeColor = 'k';%useCols(id,:);
%         sc.MarkerFaceAlpha = .4;
        %
    end
set(gca, 'xticklabels', xlabsare, 'fontsize', fontsize, 'xtick', [1,2])

xlim([.5 size(rhoPLOTS,2)+.5])

q= fdr(pvals);
%plot adjusted pvals:
for ip=1:length(q)
    
    %             loc = eh.YData(ip) + sign(eh.YData(ip))*eh.YNegativeDelta(ip) + (.1)*sign(eh.YData(ip));
    loc = eh.YData(ip) + 2.2*sign(eh.YData(ip))*eh.UData(ip);
    if q(ip)<.001
        usem = '***';
    elseif q(ip)<.01
        usem = '**';
    elseif q(ip)<.05
        usem = '*';
    else
        usem = '';
    end
    
    tc=text(ip, loc, usem,'fontsize', 40, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% if exp 2, add n/a to plot
if iExp==2
    text(2, 0, 'na','fontsize', 40, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
%% is the diff bw sig?
[h,p,CI,stat]=ttest(rhoPLOTS(:,1), rhoPLOTS(:,2));
if p<.001
    usem = '***';
elseif p<.01
    usem = '**';
elseif p<.05
    usem = '*';
else
    usem = '';
end
%cohens d = mean(diff) / SD(diff). (for paired samples ttest)
%d = m /sd (one sided).
diffs = rhoPLOTS(:,1) - rhoPLOTS(:,2);
d= mean(diffs)/std(diffs);

if p<=.05
    tc=text(1.5, .9, usem,'fontsize', 40, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    plot([1 2], [.8 .8], '-k', 'linew', 2); shg
end
%%


ylabel({['Spearman''s rho']});
ylim([-.15 1])
title(['Experiment ' num2str(iExp) ])
set(gca, 'fontsize', fontsize);
set(gcf, 'color', 'w')
box on