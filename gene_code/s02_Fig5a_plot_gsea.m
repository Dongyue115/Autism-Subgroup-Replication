clc;clear;

path = 'GSEA_result_GLM/';
s1_t = readtable([path 'ASD_related/s1_PLS2.xlsx']);
s2_t = readtable([path 'ASD_related/s2.xlsx']);
s3_t = readtable([path 'ASD_related/s3.xlsx']);
s4_t = readtable([path 'ASD_related/s4_PLS2.xlsx']);
s1 = [s1_t(6:8,:);s1_t(4,:);s1_t(2,:);s1_t(1,:);s1_t(5,:)];
s2 = [s2_t(6:8,:);s2_t(4,:);s2_t(2,:);s2_t(1,:);s2_t(5,:)];
s3 = [s3_t(6:8,:);s3_t(4,:);s3_t(2,:);s3_t(1,:);s3_t(5,:)];
s4 = [s4_t(6:8,:);s4_t(4,:);s4_t(2,:);s4_t(1,:);s4_t(5,:)];

% s1_t = readtable([path 'Nonpsychiatric/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Nonpsychiatric/s2.xlsx']);
% s3_t = readtable([path 'Nonpsychiatric/s3.xlsx']);
% s4_t = readtable([path 'Nonpsychiatric/s4_PLS2.xlsx']);
% s1 = s1_t([6,4,3,2,5,1],:);
% s2 = s2_t([6,4,3,2,5,1],:);
% s3 = s3_t([6,4,3,2,5,1],:);
% s4 = s4_t([6,4,3,2,5,1],:);

% s1_t = readtable([path 'Other neuropsychiatric/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Other neuropsychiatric/s2.xlsx']);
% s3_t = readtable([path 'Other neuropsychiatric/s3.xlsx']);
% s4_t = readtable([path 'Other neuropsychiatric/s4_PLS2.xlsx']);
% s1 = s1_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s2 = s2_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s3 = s3_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s4 = s4_t([11,9,1,7,10,5,2,4,6,3,8],:);

% s1_t = readtable([path 'Synaptic signaling/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Synaptic signaling/s2.xlsx']);
% s3_t = readtable([path 'Synaptic signaling/s3.xlsx']);
% s4_t = readtable([path 'Synaptic signaling/s4_PLS2.xlsx']);
% s1 = s1_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s2 = s2_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s3 = s3_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s4 = s4_t([5,4,3,2,7,10,9,1,11,6,8],:);

% s1_t = readtable([path 'Immune signaling/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Immune signaling/s2.xlsx']);
% s3_t = readtable([path 'Immune signaling/s3.xlsx']);
% s4_t = readtable([path 'Immune signaling/s4_PLS2.xlsx']);
% s1 = s1_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s2 = s2_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s3 = s3_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s4 = s4_t([5,1,2,4,9,8,11,3,12,10,7,6],:);

% s1_t = readtable([path 'Protein translation/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Protein translation/s2.xlsx']);
% s3_t = readtable([path 'Protein translation/s3.xlsx']);
% s4_t = readtable([path 'Protein translation/s4_PLS2.xlsx']);
% s1 = s1_t([3,4,5,2,1],:);
% s2 = s2_t([3,4,5,2,1],:);
% s3 = s3_t([3,4,5,2,1],:);
% s4 = s4_t([3,4,5,2,1],:);

val_s1 = -log10(s1.padj).*sign(s1.NES);
val_s2 = -log10(s2.padj).*sign(s2.NES);
val_s3 = -log10(s3.padj).*sign(s3.NES);
val_s4 = -log10(s4.padj).*sign(s4.NES);

plot_val = [val_s1 val_s2 val_s3 val_s4];

%% 
figure
addpath 'user/github_repo/cbrewer2/'
% cmap = cbrewer2('div', 'RdYlBu', 200, 'linear');
cmap = cbrewer2('div', 'RdYlBu');
cmap = flipud(cmap(30:240, :));   
set(gcf, 'unit', 'centimeters', 'position', [2 2 18 10]); 

imagesc(plot_val);
colormap(cmap);
caxis([-4 4])
colorbar('Ticks',[-4:2:4])
set(gca, 'XTick', [1:4], 'XTickLabel',{'S1','S2','S3','S4'},'YTick', [1:size(s1,1)],'YTickLabel',s1.pathway,...
    'Position',[0.4,0.1,0.5,0.5],'FontSize',12);


for i = 1:size(plot_val, 1)     
    for j = 1:size(plot_val, 2) 
        current_val = plot_val(i, j);
        if current_val == 0
            text_str = '0';
        else
            text_str = sprintf('%.2f', current_val);
        end
        text(j, i, text_str, ...
            'HorizontalAlignment', 'center', ...  
            'VerticalAlignment', 'middle', ...   
            'FontSize', 12, 'FontName','Arial','Color', 'k');                     
    end
end
title('Autism gene sets','FontSize', 15);
% title('Nonpsychiatric gene sets','FontSize', 15);
% title('Other neuropsychiatric gene sets','FontSize', 15);
% title('Synaptic signaling gene sets','FontSize', 15)
% title('Immune signaling gene sets','FontSize', 15)
% title('Protein translation gene sets','FontSize',12)
% set(h, 'VerticalAlignment', 'bottom', 'Margin', 5); 


