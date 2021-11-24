function [RespCats, lengthbins] = splitTrialsintoBins(inputVector, nbins)
%Split vector of trials into n bins. Output is a cell array, of trial
%indices per bin.
            

disp(['Sorting ' num2str(length(inputVector)) ' values into ' num2str(nbins) ' bins ']);
disp(['between values (min)' num2str(min(inputVector)) ' and (max) ' num2str(max(inputVector) )]);
    Qp = quantile(inputVector, nbins);
            %%
disp(['Bin boundaries are ' num2str(Qp) ]);
            
            RespCats={};
            lengthbins = nan(1, nbins);
            % collect index of the trials in each quintile.
            cat1 = find(inputVector<Qp(1));
            
            RespCats{1} = cat1;
            lengthbins(1) = length(cat1);
            for iqp = 2:length(Qp)
            %Cat2
            ind1= find(inputVector>=Qp(iqp-1));
            ind2= find(inputVector<Qp(iqp));
            %find members of both:
            cat2=intersect(ind1,ind2);
            
            RespCats{iqp} = cat2;
             lengthbins(iqp) = length(cat2);
            end
            
            %% add final 
            
            RespCats{iqp+1}= find(inputVector>=Qp(iqp));
            lengthbins(iqp+1) = length(RespCats{iqp+1});
            %%
disp(['End result, ' num2str(length(RespCats)) ' bins,  std for bin counts = ' num2str(std(lengthbins))]); 


%%
end