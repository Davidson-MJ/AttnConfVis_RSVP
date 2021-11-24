% measureBasetoPeak;
% outsourcing this script to measure base to peak, in defined windows.

%% Search Negcross and Poscross
pos_index=zeros(length(tmpERP),1); % Create an empty list for positive peaks
pos_index(find(tmpERP>0))=1; % Index of all positive points for EEG
difference=diff(pos_index); % Locate the first positive and negative point (in time) for each series of consecutive points
poscross=find(difference==1); % Return the position of all first positive points
negcross=find(difference==-1); % Return the position of all first negative points
EEGder=smooth(diff(tmpERP),5); % smooth replaced 'Meanfilt', meanfilt is a function that uses a 5 sample moving window to smooth derivative
pos_index=zeros(length(EEGder),1); % Repeat above procedure on smoothed signal (?)
pos_index(find(EEGder>0))=1; % Index of all positive points above minimum threshold %%!!!!!!!!!!!!!!!!changed to 0!!!!!!!!!!!!!
difference=diff(pos_index); % Locate first positive and negative points
peaks=find(difference==-1)+1; % Find pos ZX and neg ZX of the derivative (peaks)
troughs=find(difference==1)+1; % Find pos ZX and neg ZX of the derivative (troughs)
peaks(tmpERP(peaks)<0)=[]; % Rejects peaks below zero
troughs(tmpERP(troughs)>0)=[]; % Rejects troughs above zero

%now we have our pos and neg (peaks and troughs).
%store the P1 and N1, according to window.
peaks_sec = timeax(peaks);
troughs_sec= timeax(troughs);


%%
%% dynamically select P1 and N1 window, and store pk amplitudes.

%find positive peaks in different windows:
winds(1,:)=[0.08, 0.160]; % initial search windows (90 - 160 ms)
% winds(2,:)=[0.161, 0.260]; % initial search windows (90 - 160 ms)


% if idatatype==2 && ippant==9  % seems to be a strange shift for subset of ppants.
% %     winds(1,:) = [.161, .24];
%     winds(2,:) = [.161, .280];
% end
%%
for iwind = 1%:2
    
    P1win = onset + winds(iwind,:);
    
    %% find coincidence of pks/troughs in rel windows.
    availPK = (peaks_sec>P1win(1)) & (peaks_sec<P1win(2));
    
    %% find amplitude of these avail.
    pkamps = tmpERP(peaks(availPK));
    pkloc= peaks(find(availPK));
    
    %% correct in case there are multiple:
    mID = find(pkamps==max(pkamps));
    P1 = pkamps(mID);
    P1loc = pkloc(mID);
    
    
    
    
    
    %
    %% store these peaks (if present, some channels may not).
    if ~isempty(P1)
        if iwind==1
            P1_quintiles(ichan,iquin) = P1;
            P1loct = P1loc;
        else
            P1_quintiles(ichan,iquin) = P1;
            P2loct = P1loc;
        end
        
        
        
    else %no positive peak found, show ERP for confirmation.
        
         if iwind==1
            P1_quintiles(ichan,iquin) = nan;
            P1loct = 301; % marks at onset if missed peak.
        else
            P1_quintiles(ichan,iquin) = nan;
            P2loct = 301; 
        end
             
%         figure(10); clf
%         plot(timeax, tmpERP); hold on;
%         %addsearch window
%         Xs =[P1win(1), P1win(1), P1win(2), P1win(2)];
%         Ys = [-4 4 4 -4];
%         ph=patch(Xs, Ys, 'r');
%         ph.FaceAlpha = .1;
%         plot([onset onset], ylim, ['k:'])
%         title(['ippant ' num2str(ippant) ', no peak found in window ' num2str(iwind) '(red)'])
%         xlim([.6 1.5])
%         pause(0.5);
    end
end % diff windows.
%%
if showPeakssaved
    %%
    P1t=P1_byalpha(ichan,iquin);    
    timeax2=timeax+.3;
    figure(10); clf
    plot(timeax2, tmpERP, 'k', 'linew', 2); hold on;
  
    r=plot([1 1], [ylim ], 'k:', 'linew', 2);
  
    % add p1 peak stored.
    p2=plot([timeax2(P1loct) timeax2(P1loct)], [0 tmpERP(P1loct)], 'b-', 'linew', 3);   
    xlim([timeax2(P1loct)-.2 timeax2(P1loct)+.2])
    % show av
    
    p2av=plot([timeax2(P1loct-5) timeax2(P1loct+5)], [0 0], 'r-', 'linew', 3);    
    xlabel('Time since trial start [s]')
    ylabel('\bf\mu \rmV');
    set(gca, 'fontsize', 20)

  
    legend([r, p2], {'RSVP onset', 'P1 amplitude'});
%     ylim([-8 30])
xlim([0.7 1.5])
end
