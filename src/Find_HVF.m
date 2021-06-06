function [HVG,data_selected,ind_proposed,vec] = Find_HVG(data,geneNAME,test,param)
% data: row--samples, column--features
% 'anderson'  : Anderson-Darling test
% 'param.alpha' : Significance level of hypothesis test (default = 1e-4)
% 
% 'ske-kur'   : Skewness and Kurtosis test
% 'param.dim' : The number of dimensions to use selecting genes
%               (default = 3)
% 'param.skewness' : Selecting the principal component in the absolute
%                    range of that parameter (default = 0.5)
% 'param.kurtosis_min' :  
% 'param.s' : The number of selecting genes (default = 200)

if exist('param.dim') == 1
    param.dim = 3;
end
if exist('param.s') == 1
    param.s = 200;
end
[~,vec] = pca(data','Centered',true,'Algorithm','svd');
hp = zeros(2,size(vec,2));
switch test
    case 'anderson'
        if exist('param.alpha') == 1
            param.alpha = 1e-4;
        end
        for i = 1:size(vec,2)
            pd = fitdist(vec(:,i),'norm');
            [hp(1,i),hp(2,i)] = chi2gof(vec(:,i),'Alpha',param.alpha,'CDF',pd);
        end
        ind = find(hp(1,:) == 0);
        A = vec(:,ind(1):ind(1)+param.dim);
        ind(1)
        
    case 'ske-kur'
        if exist('param.skewness') == 1
            param.skewness = 0.5;
        end
        if exist('param.kurtosis_min') == 1
            param.kurtosis_min = 3;
        end
        % skewness
        hp(1,:) = skewness(vec);
        ind_skewness = find(abs(hp(1,:)) <= param.skewness);
        
        % kurtosis
        hp(2,:) = kurtosis(vec);
        ind_kurtosis = find(hp(2,:) >= param.kurtosis_min);
        ind = intersect(ind_skewness, ind_kurtosis);
        A = vec(:,ind(1):ind(1)+param.dim);
        ind(1)
end

B = sum(A.^2,2);
[~,ind_proposed] = sort(B,'descend');
HVG = geneNAME(ind_proposed(1:param.s));
data_selected = data(:,ind_proposed(1:param.s));

end
