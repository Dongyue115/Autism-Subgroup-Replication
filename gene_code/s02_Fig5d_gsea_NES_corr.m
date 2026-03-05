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
G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
Buch_G1 = zscore([-3.668 -2.875 -3.606 -2.849; -2.066 -1.680 -2.071 -1.911; -1.899 -1.513 -1.994 -1.662;...
    -2.318 -1.428 -2.119 -1.840; -1.261 -1.024 -1.065 -1.448; 0.984 -1.140 0.778 -1.110;...
    2.975 3.923 3.586 3.862]);


% s1_t = readtable([path 'Nonpsychiatric/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Nonpsychiatric/s2.xlsx']);
% s3_t = readtable([path 'Nonpsychiatric/s3.xlsx']);
% s4_t = readtable([path 'Nonpsychiatric/s4_PLS2.xlsx']);
% s1 = s1_t([6,4,3,2,5,1],:);
% s2 = s2_t([6,4,3,2,5,1],:);
% s3 = s3_t([6,4,3,2,5,1],:);
% s4 = s4_t([6,4,3,2,5,1],:);
% G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
% Buch_G1 = zscore([0.823 -1.068 0.890 -1.016; 0.912 -0.667 0.904 -0.788; 1.096 -0.737 1.061 0.952; ...
%     1.091 -0.836 1.302 -1.142; 1.672 1.363 1.786 1.349; 1.714 1.599 1.978 1.338]);


% s1_t = readtable([path 'Other neuropsychiatric/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Other neuropsychiatric/s2.xlsx']);
% s3_t = readtable([path 'Other neuropsychiatric/s3.xlsx']);
% s4_t = readtable([path 'Other neuropsychiatric/s4_PLS2.xlsx']);
% s1 = s1_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s2 = s2_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s3 = s3_t([11,9,1,7,10,5,2,4,6,3,8],:);
% s4 = s4_t([11,9,1,7,10,5,2,4,6,3,8],:);
% G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
% Buch_G1 = zscore([-2.395 -1.486 -2.051 -2.133; -1.468 -1.425 -1.251 -1.453; -1.567 -1.167 -1.210 -1.572;...
%     -1.549 -1.556 -1.211 -1.871; -1.520 -1.376 -0.984 -1.701; -1.286 -1.229 -0.923 -1.5; -1.323 -1.199 -1.139 -1.452;...
%     -1.086 -0.879 1.114 -1.342; 0.878 -1.006 0.835 -0.924; 0.97 -0.943 0.772 -0.985; 0.956 -0.812 1.235 -0.815]);


% s1_t = readtable([path 'Synaptic signaling/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Synaptic signaling/s2.xlsx']);
% s3_t = readtable([path 'Synaptic signaling/s3.xlsx']);
% s4_t = readtable([path 'Synaptic signaling/s4_PLS2.xlsx']);
% s1 = s1_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s2 = s2_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s3 = s3_t([5,4,3,2,7,10,9,1,11,6,8],:);
% s4 = s4_t([5,4,3,2,7,10,9,1,11,6,8],:);
% G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
% Buch_G1 = zscore([-1.943 -1.624 -1.953 -2.018; -2.237 -1.751 -2.355 -2.053; -2.162 -1.62 -1.993 -2.053;...
%     -2.1 -1.346 -1.668 -1.914; -2.408 -1.746 -2.249 -2.167; -2.439 -1.537 -2.192 -2.142; -2.334 -1.939 -2.355 -2.032;...
%     -1.703 -1.453 -1.306 -1.921;-2.457 -1.741 -2.259 -2.004; -1.953 -1.467 -1.246 -2.025;-2.214 -1.357 -1.773 -1.795]);


% s1_t = readtable([path 'Immune signaling/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Immune signaling/s2.xlsx']);
% s3_t = readtable([path 'Immune signaling/s3.xlsx']);
% s4_t = readtable([path 'Immune signaling/s4_PLS2.xlsx']);
% s1 = s1_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s2 = s2_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s3 = s3_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% s4 = s4_t([5,1,2,4,9,8,11,3,12,10,7,6],:);
% G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
% Buch_G1 = zscore([1.511 1.571 1.908 1.430; 1.590 1.937 1.869 1.609; 1.566 1.827 1.885 1.532; 1.608 2.096 1.985 1.767;...
%     1.641 2.263 1.976 1.883; 1.656 2.093 1.973 1.770; 1.871 1.775 1.971 1.269; 1.958 1.754 2.169 1.846;...
%     1.885 2.426 1.993 2.391; 1.985 1.974 2.044 1.997; 2.047 2.221 2.151 2.194; 1.974 1.812 2.166 1.858]);


% s1_t = readtable([path 'Protein translation/s1_PLS2.xlsx']);
% s2_t = readtable([path 'Protein translation/s2.xlsx']);
% s3_t = readtable([path 'Protein translation/s3.xlsx']);
% s4_t = readtable([path 'Protein translation/s4_PLS2.xlsx']);
% s1 = s1_t([3,4,5,2,1],:);
% s2 = s2_t([3,4,5,2,1],:);
% s3 = s3_t([3,4,5,2,1],:);
% s4 = s4_t([3,4,5,2,1],:);
% G1_NES = zscore([s1.NES s2.NES s3.NES s4.NES]);
% Buch_G1 = zscore([1.520 -0.504 0.937 0.749; 1.728 -0.475 1.163 0.75;1.766 -0.483 1.194 0.757;...
%     1.875 0.536 1.439 1.230; 1.836 0.713 1.275 0.998]);


[r_g1,p_g1] = corr(reshape(G1_NES,[],1), reshape(Buch_G1,[],1));
cos_sim = 1 - pdist2(G1_NES(:)',Buch_G1(:)', 'cosine');

figure
set(gcf,'position',[100 100 400 350]);
[P_plot,~] = polyfit(reshape(G1_NES,[],1),reshape(Buch_G1,[],1),1);
colors = [0.753 0.259 0.357;0.282 0.718 0.388;0.420 0.310 0.616;0.369 0.678 0.878];
groupnames = {'S1','S2','S3','S4'};


for i = 1:size(G1_NES,2)
    h1(i) = scatter(G1_NES(:,i),Buch_G1(:,i),32,'filled','MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:));
    hold on
    for j = 1:size(G1_NES,1)
        text(G1_NES(j,i),Buch_G1(j,i),num2str(j),"FontSize",10,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    end
end
xlim([-2.2 3])
ylim([-2.2 3])
xlabel('LEAP NES')
ylabel('ABIDE NES')
set(gca,'Fontsize',12)

plot(reshape(G1_NES,[],1),polyval(P_plot,reshape(G1_NES,[],1)),'Color',[0.388 0.388 0.388],'LineWidth',1.5);

if p_g1<0.001
    p_str = 'p<0.001';
else
    p_str = strcat('p= ',num2str(p_g1,'%.3f'));
end
text('string',p_str,'Units','normalized','position',[0.05,0.83],'Color',[0.388 0.388 0.388],'Fontsize',12,'FontWeight','bold');
text('string',strcat('r =',num2str(r_g1,'%.3f')),'Units','normalized','position',[0.05,0.91],'Color',[0.388 0.388 0.388],'Fontsize',12,'FontWeight','bold');

legend(h1,groupnames,'Location','southeast');
legend('boxoff')
box(gca,'on')   

title('Autism gene sets','FontSize',15);
% title('Nonpsychiatric gene sets','FontSize',15);
% title('Other neuropsychiatric gene sets','FontSize',15);
% title('Synaptic signaling gene sets','FontSize',15);
% title('Immune signaling gene sets','FontSize',15);
% title('Protein translation gene sets','FontSize',15);
