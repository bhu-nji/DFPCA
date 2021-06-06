clear; rng(0,'twister')

% Load dataset
load('simulation_dataset.mat');

% Load labels of samples
load sample_labels.mat

% Data normalization
data.scaled = normalize(X,2,'zscore');

% Create the gene names
for i = 1:size(X,2)
    geneNAME.all(i,1) = cellstr(append('G',num2str(i)));
end

% Selecting genes (Distortion-free PCA)
param.skewness = 0.5;
param.kurtosis_min = 0;
param.dim = 20;
param.s = 200;
test = 'ske-kur';
[geneNAME.selected,data.selected,index.sorted_genes] = Find_HVG(data.scaled,...
    geneNAME.all,test,param);
data.selected = data.scaled(:,index.sorted_genes(1:param.s));

% PCA-UMAP
dim = 7;
[~,Z.pca_proposed] = pca(data.selected,'Centered',true);
Z.umap_proposed = run_umap(Z.pca_proposed(:,1:dim),'verbose','none');

% Visualization
Legend = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'};
col = [1 0 0; 0 0 1; 0 1 0; 1 0 1; 0 1 1; 0.8,0.8,0.8];

figure; ax = [1 2 3]; set(gcf,'color','white');
gscatter(Z.umap_proposed(:,ax(1)),Z.umap_proposed(:,ax(2)),label,col,'+o^sd',6)
legend(Legend{:});
title('Distortion-free PCA','FontSize',14);
xlabel(['UMAP ',num2str(ax(1))]); ylabel(['UMAP ',num2str(ax(2))]);
set(gca, 'FontName','Times','FontSize',14); 