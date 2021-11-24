% ANOVAERPfeats
figure(ifeat+10); 
fontsize= 20;
set(gcf, 'units', 'normalized', 'position', [.5 1 .5 1], 'color', 'w'); shg
            
 %mean topoplot of the effect first.
            usechans=CPchans;
            
            subplot(221);
            tpPxCh = squeeze(nanmean(GFX_alphavsERPfeats(:, ifeat, :, :),4));
            %topoplot with subchans.
            topoplot(squeeze(nanmean(tpPxCh,1)), elocs32, 'emarker2', {usechans '*' 'w', 10, 3});
            c=colorbar;
            if normon==1
            ylabel(c, ' normalized amplitude');
            else
                ylabel(c, ' amplitude');
            end                
            hold on
            title(featsare{ifeat});
            set(gca, 'fontsize', fontsize);
            %% now perform rm_anova on data.
%             rmanova_return = rmanova(data,factorarray,subjects,varnames, btw_ss_col)

dataALL = squeeze(GFX_alphavsERPfeats(:, ifeat, :, :));
[nsubs, nchans, nIV]=size(dataALL);

%specify design:
factorarray =[ones(nsubs,1);repmat(2,[nsubs,1]);...
    repmat(3,[nsubs,1]);repmat(4,[nsubs,1]);repmat(5,[nsubs,1])];

subjects = repmat([1:nsubs]', [nIV,1]);


%output: 
pvals = nan(nchans,1);
Fvals = nan(nchans,1);
%%
for ichan = 1:nchans
    
    datatmp = squeeze(dataALL(:,ichan,:)); % needs to be a columnvector:
    
    datacol = reshape(datatmp, [nsubs*nIV, 1]);
    try
    %perform rmanova
    anovtab = rmanova(datacol, factorarray, subjects);
    
    pvals(ichan) = anovtab.table{2,7};
    Fvals(ichan) = anovtab.table{2,6};
    catch
    end
end
%%
subplot(2,2,2);
%mask nonsig: 
allel = ones(1,length(Fvals));
%nonsig
% ns = find(pvals<.05);
% allel(ns)=0;
topoplot(Fvals, elocs32, 'pmask', allel); title('Main effect: Alpha'); c=colorbar;
ylabel(c, 'F stat');
caxis([0 5]);
set(gca, 'fontsize', fontsize);

%%
subplot(2,2,4);
pvals(isnan(pvals))=1.1;
%mask nonsig: 
allel = ones(1,length(Fvals));
%nonsig
ns = find(pvals<.05);
allel(ns)=0;
topoplot(log(pvals), elocs32, 'pmask', allel); title('log(pval)'); c=colorbar;
caxis([log(.001) log(.05)]);

ylabel(c, 'log(p)');
set(gcf, 'color', 'w');
set(gca, 'fontsize', fontsize);

%%
%%
colormap('viridis')
        cd(figuredir)
print('-dpng', ['TOPOrmANOVA for ' featsare{ifeat} ' ' xlabis]);
