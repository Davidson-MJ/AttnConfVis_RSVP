function RESPColours = brewCOLOURS

% select the colour pallete for all plots.

%using the cbrewer package. 

%Confidence
confspec= cbrewer('div', 'Spectral', 5);

%Visiblity
visspec= cbrewer('seq', 'YlGn', 5);

%Attention
attspec= cbrewer('seq', 'Reds', 5);

%Alpha 

alphaspec= cbrewer('seq', 'Purples', 5);

RESPColours=[];
RESPColours(1,:,:) = flipud(confspec);
RESPColours(2,:,:) = flipud(visspec);
RESPColours(3,:,:) = flipud(attspec);
RESPColours(4,:,:) = flipud(alphaspec);

end