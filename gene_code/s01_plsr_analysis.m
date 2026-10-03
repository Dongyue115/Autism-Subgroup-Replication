clc;clear;

path ='project/';
genepath = 'replicate/';

gene_names = readcell([genepath 'expression_243ROI.csv'], 'Range', '1:1'); 
gene_data = readmatrix([genepath 'expression_243ROI.csv'], 'NumHeaderLines', 1);

LEAP_ROI = readtable([path 'atlas_group.xlsx'],'Sheet','LEAP_243ROI');

nan_row = any(isnan(gene_data),2);
del_roi = find(nan_row==1);
reorder_gene = gene_data(nan_row==0,:);
saveT = LEAP_ROI;
saveT(del_roi,:) = [];

asdpath = '/user/FC_asd_data/';
conpath = '/user/FC_con_data/';

load('/user/T_sub.mat');
T.t1_sex(T.t1_sex==-1) = 0;%sex:0female;1male

asd_idx = find(T.t1_diagnosis==2);
con_idx = find(T.t1_diagnosis==1);
T_asd = T(asd_idx,:);
T_con = T(con_idx,:);

mask = tril(ones(size(reorder_gene,1)),-1);
for s = 1:size(T_asd,1)
    load([asdpath strcat(num2str(T_asd.subjects(s)),'_corr.mat')]);
    zFC(del_roi,:) = [];
    zFC(:,del_roi) = [];
    vector_fc = zFC(mask==1);
    allfc_asd(s,:) = vector_fc';
end
for s = 1:size(T_con,1)
    load([conpath strcat(num2str(T_con.subjects(s)),'_corr.mat')]);
    zFC(del_roi,:) = [];
    zFC(:,del_roi) = [];
    vector_fc = zFC(mask==1);
    allfc_con(s,:) = vector_fc';
end

newdata = zscore([allfc_asd;allfc_con]);
z_asd_fc = newdata([1:size(T_asd,1)],:);
z_con_fc = newdata([1+size(T_asd,1):end],:);

respath = '/user/RCCA_results_243ROI/';
load([respath 'idx_all_4cluster.mat']);
c1 = idx_all{3};
c2 = idx_all{2};
c3 = idx_all{4};
c4 = idx_all{1};

%% GLM
T_asd.t1_diagnosis = 1 * ones(length(T_asd.t1_diagnosis),1);%asd1,con0
T_con.t1_diagnosis = 0 * ones(length(T_con.t1_diagnosis),1);

T1_group = [T_asd.t1_diagnosis(c1);T_con.t1_diagnosis];
T1_age = [T_asd.t1_ageyrs(c1);T_con.t1_ageyrs];
T1_sex = [T_asd.t1_sex(c1);T_con.t1_sex];
T1_site = dummyvar([T_asd.t1_site(c1);T_con.t1_site]);
T1_fiq = [T_asd.t1_fsiq(c1);T_con.t1_fsiq];
T1_FD = [T_asd.meanFD(c1);T_con.meanFD];
regressor1 = [ones(length(T1_group),1) T1_group T1_age T1_sex T1_fiq T1_FD T1_site(:,1:4)];

T2_group = [T_asd.t1_diagnosis(c2);T_con.t1_diagnosis];
T2_age = [T_asd.t1_ageyrs(c2);T_con.t1_ageyrs];
T2_sex = [T_asd.t1_sex(c2);T_con.t1_sex];
T2_site = dummyvar([T_asd.t1_site(c2);T_con.t1_site]);
T2_fiq = [T_asd.t1_fsiq(c2);T_con.t1_fsiq];
T2_FD = [T_asd.meanFD(c2);T_con.meanFD];
regressor2 = [ones(length(T2_group),1) T2_group T2_age T2_sex T2_fiq T2_FD T2_site(:,1:4)];

T3_group = [T_asd.t1_diagnosis(c3);T_con.t1_diagnosis];
T3_age = [T_asd.t1_ageyrs(c3);T_con.t1_ageyrs];
T3_sex = [T_asd.t1_sex(c3);T_con.t1_sex];
T3_site = dummyvar([T_asd.t1_site(c3);T_con.t1_site]);
T3_fiq = [T_asd.t1_fsiq(c3);T_con.t1_fsiq];
T3_FD = [T_asd.meanFD(c3);T_con.meanFD];
regressor3 = [ones(length(T3_group),1) T3_group T3_age T3_sex T3_fiq T3_FD T3_site(:,1:4)];

T4_group = [T_asd.t1_diagnosis(c4);T_con.t1_diagnosis];
T4_age = [T_asd.t1_ageyrs(c4);T_con.t1_ageyrs];
T4_sex = [T_asd.t1_sex(c4);T_con.t1_sex];
T4_site = dummyvar([T_asd.t1_site(c4);T_con.t1_site]);
T4_fiq = [T_asd.t1_fsiq(c4);T_con.t1_fsiq];
T4_FD = [T_asd.meanFD(c4);T_con.meanFD];
regressor4 = [ones(length(T4_group),1) T4_group T4_age T4_sex T4_fiq T4_FD T4_site(:,1:4)];

Contrast = zeros(1, size(regressor1, 2));
Contrast(2) = 1;

addpath '/user/GRETNA/Stat/';
[~,~, t1,~] = gretna_GroupAnalysis([z_asd_fc(c1,:);z_con_fc], regressor1, Contrast, 'T');
[~,~, t2,~] = gretna_GroupAnalysis([z_asd_fc(c2,:);z_con_fc], regressor2, Contrast, 'T');
[~,~, t3,~] = gretna_GroupAnalysis([z_asd_fc(c3,:);z_con_fc], regressor3, Contrast, 'T');
[~,~, t4,~] = gretna_GroupAnalysis([z_asd_fc(c4,:);z_con_fc], regressor4, Contrast, 'T');

roi_num = size(reorder_gene,1);
t1_new = zeros(roi_num ,roi_num);
t1_new(mask==1) = t1;
t1_fc = t1_new + t1_new';
t1_region_fc = sum(t1_fc,1);
t2_new = zeros(roi_num ,roi_num);
t2_new(mask==1) = t2;
t2_fc = t2_new + t2_new';
t2_region_fc = sum(t2_fc,1);
t3_new = zeros(roi_num ,roi_num);
t3_new(mask==1) = t3;
t3_fc = t3_new + t3_new';
t3_region_fc = sum(t3_fc,1);
t4_new = zeros(roi_num ,roi_num);
t4_new(mask==1) = t4;
t4_fc = t4_new + t4_new';
t4_region_fc = sum(t4_fc,1);

%% plot explained variance 
dim = 10;
[pXL,pYL,pXS,pYS,pBETA,pPCTVAR,pMSE,pstats] = plsregress(zscore(reorder_gene),zscore(t1_region_fc'),dim);

figure
plot(1:dim,100*pPCTVAR(2,1:dim),'-o','LineWidth',2,'Color',[140/255,0,0]);
set(gca,'linewidth',2);
set( gca, 'Position',[ 0.3, 0.3, 0.5, 0.5]);
box off
xlabel('Number of PLS Components','fontsize',12);
ylabel('Variance Explained in Y','fontsize',12);
set(gca,'Fontsize',12);
ylim([0 30])
title('S1')

%%  assess the significance of PLS component
perm_num = 10000;
X = zscore(reorder_gene);
Y = zscore(t1_region_fc');%change subgroup
for dim = 1:5
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,dim);
    temp = cumsum(100*PCTVAR(2,1:dim));
    Rsquared = temp(dim);
    for j = 1:perm_num
        order = randperm(size(Y,1));
        Yp = Y(order,:);
        [pXL,pYL,pXS,pYS,pBETA,pPCTVAR,pMSE,pstats] = plsregress(X,Yp,dim);
        ptemp = cumsum(100*pPCTVAR(2,1:dim));
        Rsq(j) = ptemp(dim);
    end
dim
R_pls(dim) = Rsquared;
p_pls(dim) = length(find(Rsq>=Rsquared))/perm_num;
end

%% correaltion between pls component and FC
% plot in s01_spintest.m

%% output genelist
output_dir = 'PLSR_result/output_gene_order_fromGLM/';
output_filename = {'S1_PLS2_geneWeights.csv'};
fc_stat = t1_region_fc';
[geneindex1 PLS1 Z1] = d_PLS_bootstrap_forlist(fc_stat,reorder_gene,roi_num,gene_names);

%% print out results
fid1 = fopen(fullfile(output_dir,output_filename{1}),'w');
for i=1:length(gene_names)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1);
