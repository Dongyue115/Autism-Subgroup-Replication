clc;clear;

addpath 'user/MultipleTestingToolbox'
path ='/project/';

LEAP_ROI = readtable([path 'atlas_group.xlsx'],'Sheet','LEAP_243ROI');
reorderT = sortrows(LEAP_ROI,{'network_order','Group_Order_Buch','Number'},{'ascend','ascend','ascend'});

asdpath = '/user/FC_asd_data/';
conpath = '/user/FC_con_data/';

load('/user/T_sub.mat');
T.t1_sex(T.t1_sex==-1) = 0;%sex:0female;1male

asd_idx = find(T.t1_diagnosis==2);
con_idx = find(T.t1_diagnosis==1);
T_asd = T(asd_idx,:);
T_con = T(con_idx,:);

mask1 = flipud(tril(ones(243),-1));
for a1 = 1:size(T_asd,1)
    load([asdpath strcat(num2str(T_asd.subjects(a1)),'_corr.mat')]);
    asd_fc_conn{a1} = fc;
    vector_fc = fc(mask1==1);
    asd_fc(a1,:) = vector_fc';
end
for a2 = 1:size(T_con,1)
    load([conpath strcat(num2str(T_con.subjects(a2)),'_corr.mat')]);
    con_fc_conn{a2} = fc;
    vector_fc = fc(mask1==1);
    con_fc(a2,:) = vector_fc';
end
newdata = zscore([asd_fc;con_fc]);
z_asd_fc = newdata([1:size(T_asd,1)],:);
z_con_fc = newdata([1+size(T_asd,1):end],:);

respath = '/user/RCCA_results_243ROI/';
load([respath 'idx_all_4cluster.mat']);
c1 = idx_all{3};
c2 = idx_all{2};
c3 = idx_all{4};
c4 = idx_all{1};

%% GLM
addpath '/user/GRETNA/Stat/';
roi_num = 243;

T_asd.t1_diagnosis = 1 * ones(length(T_asd.t1_diagnosis),1);%asd1,con0
T_con.t1_diagnosis = 0 * ones(length(T_con.t1_diagnosis),1);

T1_group = [T_asd.t1_diagnosis(c1);T_con.t1_diagnosis];%S1-CON
T1_age = [T_asd.t1_ageyrs(c1);T_con.t1_ageyrs];
T1_sex = [T_asd.t1_sex(c1);T_con.t1_sex];
T1_site = dummyvar([T_asd.t1_site(c1);T_con.t1_site]);
T1_fiq = [T_asd.t1_fsiq(c1);T_con.t1_fsiq];
T1_FD = [T_asd.meanFD(c1);T_con.meanFD];
regressor1 = [ones(length(T1_group),1) T1_group T1_age T1_sex T1_fiq T1_FD T1_site(:,1:4)];

Contrast_t = zeros(1, size(regressor1, 2));
Contrast_t(2) = 1;
[~,~, t1,~] = gretna_GroupAnalysis([z_asd_fc(c1,:);z_con_fc], regressor1, Contrast_t, 'T');
TTest2_P1 = 2 *(1-tcdf(abs(t1), size([z_asd_fc(c1,:);z_con_fc],1)-size(regressor1, 2))); 

T2_group = [T_asd.t1_diagnosis(c2);T_con.t1_diagnosis];%S2-CON
T2_age = [T_asd.t1_ageyrs(c2);T_con.t1_ageyrs];
T2_sex = [T_asd.t1_sex(c2);T_con.t1_sex];
T2_site = dummyvar([T_asd.t1_site(c2);T_con.t1_site]);
T2_fiq = [T_asd.t1_fsiq(c2);T_con.t1_fsiq];
T2_FD = [T_asd.meanFD(c2);T_con.meanFD];
regressor2 = [ones(length(T2_group),1) T2_group T2_age T2_sex T2_fiq T2_FD T2_site(:,1:4)];
[~,~, t2,~] = gretna_GroupAnalysis([z_asd_fc(c2,:);z_con_fc], regressor2, Contrast_t, 'T');
TTest2_P2 = 2 *(1-tcdf(abs(t2), size([z_asd_fc(c2,:);z_con_fc],1)-size(regressor2, 2))); 

T3_group = [ones(length(c1),1);zeros(length(c2),1)];%S1-S2
T3_age = [T_asd.t1_ageyrs(c1);T_asd.t1_ageyrs(c2)];
T3_sex = [T_asd.t1_sex(c1);T_asd.t1_sex(c2)];
T3_site = dummyvar([T_asd.t1_site(c1);T_asd.t1_site(c2)]);
T3_fiq = [T_asd.t1_fsiq(c1);T_asd.t1_fsiq(c2)];
T3_FD = [T_asd.meanFD(c1);T_asd.meanFD(c2)];
regressor3 = [ones(length(T3_group),1) T3_group T3_age T3_sex T3_fiq T3_FD T3_site(:,1:4)];
[~,~, t3,~] = gretna_GroupAnalysis([z_asd_fc(c1,:);z_asd_fc(c2,:)], regressor3, Contrast_t, 'T');
TTest2_P3 = 2 *(1-tcdf(abs(t3), size([z_asd_fc(c1,:);z_asd_fc(c2,:)],1)-size(regressor3, 2)));

T4_group = [T_asd.t1_diagnosis(c3);T_con.t1_diagnosis];%S3-CON
T4_age = [T_asd.t1_ageyrs(c3);T_con.t1_ageyrs];
T4_sex = [T_asd.t1_sex(c3);T_con.t1_sex];
T4_site = dummyvar([T_asd.t1_site(c3);T_con.t1_site]);
T4_fiq = [T_asd.t1_fsiq(c3);T_con.t1_fsiq];
T4_FD = [T_asd.meanFD(c3);T_con.meanFD];
regressor4 = [ones(length(T4_group),1) T4_group T4_age T4_sex T4_fiq T4_FD T4_site(:,1:4)];
[~,~, t4,~] = gretna_GroupAnalysis([z_asd_fc(c3,:);z_con_fc], regressor4, Contrast_t, 'T');
TTest2_P4 = 2 *(1-tcdf(abs(t4), size([z_asd_fc(c3,:);z_con_fc],1)-size(regressor4, 2))); 

T5_group = [T_asd.t1_diagnosis(c4);T_con.t1_diagnosis];%S4-CON
T5_age = [T_asd.t1_ageyrs(c4);T_con.t1_ageyrs];
T5_sex = [T_asd.t1_sex(c4);T_con.t1_sex];
T5_site = dummyvar([T_asd.t1_site(c4);T_con.t1_site]);
T5_fiq = [T_asd.t1_fsiq(c4);T_con.t1_fsiq];
T5_FD = [T_asd.meanFD(c4);T_con.meanFD];
regressor5 = [ones(length(T5_group),1) T5_group T5_age T5_sex T5_fiq T5_FD T5_site(:,1:4)];
[~,~, t5,~] = gretna_GroupAnalysis([z_asd_fc(c4,:);z_con_fc], regressor5, Contrast_t, 'T');
TTest2_P5 = 2 *(1-tcdf(abs(t5), size([z_asd_fc(c4,:);z_con_fc],1)-size(regressor5, 2)));

T6_group = [ones(length(c1),1);zeros(length(c3),1)];%S1-S3
T6_age = [T_asd.t1_ageyrs(c1);T_asd.t1_ageyrs(c3)];
T6_sex = [T_asd.t1_sex(c1);T_asd.t1_sex(c3)];
T6_site = dummyvar([T_asd.t1_site(c1);T_asd.t1_site(c3)]);
T6_fiq = [T_asd.t1_fsiq(c1);T_asd.t1_fsiq(c3)];
T6_FD = [T_asd.meanFD(c1);T_asd.meanFD(c3)];
regressor6 = [ones(length(T6_group),1) T6_group T6_age T6_sex T6_fiq T6_FD T6_site(:,1:4)];
[~,~, t6,~] = gretna_GroupAnalysis([z_asd_fc(c1,:);z_asd_fc(c3,:)], regressor6, Contrast_t, 'T');
TTest2_P6 = 2 *(1-tcdf(abs(t6), size([z_asd_fc(c1,:);z_asd_fc(c3,:)],1)-size(regressor6, 2)));

T7_group = [ones(length(c1),1);zeros(length(c4),1)];%S1-S4
T7_age = [T_asd.t1_ageyrs(c1);T_asd.t1_ageyrs(c4)];
T7_sex = [T_asd.t1_sex(c1);T_asd.t1_sex(c4)];
T7_site = dummyvar([T_asd.t1_site(c1);T_asd.t1_site(c4)]);
T7_fiq = [T_asd.t1_fsiq(c1);T_asd.t1_fsiq(c4)];
T7_FD = [T_asd.meanFD(c1);T_asd.meanFD(c4)];
regressor7 = [ones(length(T7_group),1) T7_group T7_age T7_sex T7_fiq T7_FD T7_site(:,1:4)];
[~,~, t7,~] = gretna_GroupAnalysis([z_asd_fc(c1,:);z_asd_fc(c4,:)], regressor7, Contrast_t, 'T');
TTest2_P7 = 2 *(1-tcdf(abs(t7), size([z_asd_fc(c1,:);z_asd_fc(c4,:)],1)-size(regressor7, 2)));

T8_group = [ones(length(c2),1);zeros(length(c3),1)];%S2-S3
T8_age = [T_asd.t1_ageyrs(c2);T_asd.t1_ageyrs(c3)];
T8_sex = [T_asd.t1_sex(c2);T_asd.t1_sex(c3)];
T8_site = dummyvar([T_asd.t1_site(c2);T_asd.t1_site(c3)]);
T8_fiq = [T_asd.t1_fsiq(c2);T_asd.t1_fsiq(c3)];
T8_FD = [T_asd.meanFD(c2);T_asd.meanFD(c3)];
regressor8 = [ones(length(T8_group),1) T8_group T8_age T8_sex T8_fiq T8_FD T8_site(:,1:4)];
[~,~, t8,~] = gretna_GroupAnalysis([z_asd_fc(c2,:);z_asd_fc(c3,:)], regressor8, Contrast_t, 'T');
TTest2_P8 = 2 *(1-tcdf(abs(t8), size([z_asd_fc(c2,:);z_asd_fc(c3,:)],1)-size(regressor8, 2)));

T9_group = [ones(length(c2),1);zeros(length(c4),1)];%S2-S4
T9_age = [T_asd.t1_ageyrs(c2);T_asd.t1_ageyrs(c4)];
T9_sex = [T_asd.t1_sex(c2);T_asd.t1_sex(c4)];
T9_site = dummyvar([T_asd.t1_site(c2);T_asd.t1_site(c4)]);
T9_fiq = [T_asd.t1_fsiq(c2);T_asd.t1_fsiq(c4)];
T9_FD = [T_asd.meanFD(c2);T_asd.meanFD(c4)];
regressor9 = [ones(length(T9_group),1) T9_group T9_age T9_sex T9_fiq T9_FD T9_site(:,1:4)];
[~,~, t9,~] = gretna_GroupAnalysis([z_asd_fc(c2,:);z_asd_fc(c4,:)], regressor9, Contrast_t, 'T');
TTest2_P9 = 2 *(1-tcdf(abs(t9), size([z_asd_fc(c2,:);z_asd_fc(c4,:)],1)-size(regressor9, 2)));

T10_group = [ones(length(c3),1);zeros(length(c4),1)];%S3-S4
T10_age = [T_asd.t1_ageyrs(c3);T_asd.t1_ageyrs(c4)];
T10_sex = [T_asd.t1_sex(c3);T_asd.t1_sex(c4)];
T10_site = dummyvar([T_asd.t1_site(c3);T_asd.t1_site(c4)]);
T10_fiq = [T_asd.t1_fsiq(c3);T_asd.t1_fsiq(c4)];
T10_FD = [T_asd.meanFD(c3);T_asd.meanFD(c4)];
regressor10 = [ones(length(T10_group),1) T10_group T10_age T10_sex T10_fiq T10_FD T10_site(:,1:4)];
[~,~, t10,~] = gretna_GroupAnalysis([z_asd_fc(c3,:);z_asd_fc(c4,:)], regressor10, Contrast_t, 'T');
TTest2_P10 = 2 *(1-tcdf(abs(t10), size([z_asd_fc(c3,:);z_asd_fc(c4,:)],1)-size(regressor10, 2)));

pnum = length(TTest2_P1);
c_pval_t = mafdr([TTest2_P1;TTest2_P2;TTest2_P3;TTest2_P4;TTest2_P5;TTest2_P6;TTest2_P7;TTest2_P8;TTest2_P9;...
    TTest2_P10],'BHFDR',true);

c_pval1 = c_pval_t(1:pnum);
c_pval2 = c_pval_t(1+pnum:pnum*2);
c_pval3 = c_pval_t(1+pnum*2:pnum*3);
c_pval4 = c_pval_t(1+pnum*3:pnum*4);
c_pval5 = c_pval_t(1+pnum*4:pnum*5);
c_pval6 = c_pval_t(1+pnum*5:pnum*6);
c_pval7 = c_pval_t(1+pnum*6:pnum*7);
c_pval8 = c_pval_t(1+pnum*7:pnum*8);
c_pval9 = c_pval_t(1+pnum*8:pnum*9);
c_pval10 = c_pval_t(1+pnum*9:pnum*10);

c_pval1 = c_pval1;%% change
t1 = t1;%% change

% thr = 0.001;
thr = 0.05;
cpval_new1 = ones(roi_num ,roi_num);
cpval_new1(mask1==1) = c_pval1;
t1_new = zeros(roi_num ,roi_num);
t1_new(mask1==1) = t1;
cpval_new1_plot = zeros(roi_num ,roi_num);
cpval_new1_plot(cpval_new1<thr)=1;
t1_plot = cpval_new1_plot .* t1_new + fliplr(flipud((cpval_new1_plot .* t1_new)'));

%% merge to 39 regions
to_merge = t1_plot;
reverse_fc = flipud(to_merge);%positive_order

mapping = reorderT.Group_Order_Buch;
region_num = length(unique(reorderT.Group_Order_Buch));
newfc = zeros(region_num,region_num);

for i = 1:region_num
    for j = i:region_num
        rows = find(mapping == i);
        cols = find(mapping == j);
        smallfc = reverse_fc(rows, cols);
        newfc(i, j) = mean(smallfc(smallfc~=0));
        newfc(j, i) = newfc(i, j); 
    end
end

mask3 = ones(39);
mask3(eye(39) == 1) = 0;
fc_dis_r = newfc .* mask3;
fc_dis_r(isnan(fc_dis_r)) = 0;
plot_matrix = flipud(fc_dis_r);

%% PLOT matrix
load('colorbar(bwr).mat');
figure
set(gcf, 'unit', 'centimeters', 'position', [2 2 19 19]); 
ax3 = axes('Position', [0.16, 0.23, 0.7, 0.7]); 
imagesc(plot_matrix);
clim([-5 5]);
colormap(ax3,map);

% h = colorbar;
% h.FontSize = 14;
title('S1 VS CON','FontSize',18);
set(gca, 'XTick', [], 'YTick', []);

hold on
%% network assignment
NumNet = 39;
nLim = 9;
nDMN = 8;
nFPTC = 3;
nSN = 4;
nCOTC = 3;
nSMN = 2;
nAud = 1;
nVis = 5;
nCBM = 1;
nBS = 3;

xline([nLim+0.5 nLim+nDMN+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+0.5],'-','LineWidth',1.2)
yline([nSMN+nAud+nVis+nCBM+nBS+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+nDMN+0.5],'-','LineWidth',1.2)

load('colorbar_brain.mat')
ax1 = axes('Position', [0.12, 0.23, 0.03, 0.7]);
imagesc([1:39]'); 
colors = [repmat(brain_map(10,:),nLim,1);repmat(brain_map(9,:),nDMN,1);repmat(brain_map(8,:),nFPTC,1);...
    repmat(brain_map(7,:),nSN,1);repmat(brain_map(6,:),nCOTC,1);repmat(brain_map(5,:),nSMN,1);...
    repmat(brain_map(4,:),nAud,1);repmat(brain_map(3,:),nVis,1);repmat(brain_map(2,:),nCBM,1);repmat(brain_map(1,:),nBS,1)];
colormap(ax1, colors);

region_name = unique(reorderT.Group_Buch,'stable');
net_name = unique(reorderT.network_Buch,'stable');
set(gca, 'XTick',[], 'YTick',[1:39],'YTicklabel',region_name,'TickLength',[0,0],'FontSize',10);
axis tight;
set(gca, 'YDir', 'normal'); 

ax2 = axes('Position', [0.16, 0.19, 0.7, 0.03]);
imagesc([1:39]); 
colormap(ax2, colors);
set(gca, 'XTick', [nLim/2+0.5 nLim+nDMN/2+0.5 nLim+nDMN+nFPTC/2+0.5 nLim+nDMN+nFPTC+nSN/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis+nCBM/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis+nCBM+nBS/2+0.5],'XTicklabel',net_name, 'XTickLabelRotation', 45,...
    'YTick', [],'TickLength',[0,0],'FontSize',12);
axis tight;

