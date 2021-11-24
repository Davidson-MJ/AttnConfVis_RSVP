% showP1N1_SSVEPintrialERP

%called from calcER_featsperrtial

avERPplot = squeeze(nanmean(avERP.avg,1));
%
%our relevant time windows can be plotted:
figure(1); clf
timeax=avERP.time;

plot(timeax, avERPplot), xlabel('Time (s)')
%% dynamically select P1 and N1 window, and plot result.
P1win = onset + [0.09, 0.140]; % initial search windows:
N1win = onset + [0.140, 0.190];

P1t = dsearchn(timeax', P1win');
N1t = dsearchn(timeax', N1win');

%create search vectors
P1tvec = P1t(1):P1t(2);
N1tvec = N1t(1):N1t(2);

%find peaks in these data windows.
P1pk = find(diff(avERPplot(P1tvec))<0, 1,'first');
N1pk = find(diff(avERPplot(N1tvec))>0, 1,'first');
%%
% adapt per ppant +- 30ms
P1pktime= timeax(P1tvec(P1pk));
N1pktime= timeax(N1tvec(N1pk));

%what are the new windows, per ppant:
P1win = P1pktime + [-.03, +.03]; % i
N1win = N1pktime + [-.03, +.03];

%now plot the new windows
P1t = dsearchn(timeax', P1win');
N1t = dsearchn(timeax', N1win');

SSVEPwin = P1pktime + [.2, 1.2];
SSVEPt = dsearchn(timeax', SSVEPwin');
xlim([.5 2]);
hold on;
%% plot the average as sanity check.
hold on; plot(xlim, [0 0 ], ['k-'])
plot(timeax(SSVEPt(1):SSVEPt(2)), avERPplot(SSVEPt(1):SSVEPt(2)), 'k', 'linew', 3)

plot([P1pktime P1pktime], [0 avERPplot(P1tvec(P1pk))], ['b-'], 'linew', 3) 
plot([N1pktime N1pktime], [0 avERPplot(N1tvec(N1pk))], ['r-'], 'linew', 3) 


%%
shg
pause(0.5)