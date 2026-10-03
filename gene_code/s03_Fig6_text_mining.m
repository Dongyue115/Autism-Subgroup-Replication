clc;clear;

path = 'PPI_genes/text_mining/';

s1_freq = flipud(load([path 'S1.mat']).freq);
s2_freq = flipud(load([path 'S2.mat']).freq);
s3_freq = flipud(load([path 'S3.mat']).freq);
s4_freq = flipud(load([path 'S4.mat']).freq);

plot_data = [s1_freq';s2_freq';s3_freq';s4_freq'];

figure('Position', [100 100 800 400]);
h = barh(plot_data, 'stacked');
h(1).FaceColor = [0.0784 0.4039 0.3725];
h(2).FaceColor = [0.1922 0.5529 0.5176];
h(3).FaceColor = [0.3804 0.7333 0.698];
h(4).FaceColor = [0.5922 0.8392 0.8275];
h(5).FaceColor = [0.8275 0.3961 0.1569];
h(6).FaceColor = [0.9176 0.5412 0.2588];
h(7).FaceColor = [0.9804 0.6706 0.4353];

for i = 1:length(h)
    h(i).EdgeColor = 'none';
end
set(gca, 'YDir', 'reverse');
set(gca, 'Xtick',[0:0.25:1],'Ytick',[1:4],'Yticklabel',{'S1','S2','S3','S4'},'Fontsize',14)
set(gca, 'TickDir', 'out');
ylim([0.3, size(plot_data,1) + 0.7]);
ylabel('Subgroup');
title('Relative frequency of gene-phenotype association','FontSize',18)
box off;
