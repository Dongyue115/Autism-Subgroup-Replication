clc;clear;

addpath '/user/munkres';

load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
behv = [T_asd.t1_viq T_asd.t1_sa_css_all T_asd.t1_rrb_css_all];

respath = '/user/RCCA_results_243ROI/';
load([respath 'para_com_out.mat']);
para_com = para_com_out;

%% CV alignment
load('/user/sub1000_split.mat');
ref_behv = para_com{1,1}.train_coeff_B;
ref_brain = para_com{1,1}.train_coeff_A;
ref_cv_behv = para_com{1,1}.out_train_cv_behv;%ref
ref_cv_brain = para_com{1,1}.out_train_cv_brain;%ref

for c = 1:length(para_com)
    cv_behv = para_com{1,c}.out_train_cv_behv;
    cv_brain = para_com{1,c}.out_train_cv_brain;

    trsub_behv = behv(train_idx1(c,:),:);
    m_trsub_behv = mean(trsub_behv);
    std_trsub_behv = std(trsub_behv);
    z_trsub_behv = (trsub_behv - m_trsub_behv)./std_trsub_behv;

    [r_old,p_old] = corr(cv_behv,z_trsub_behv,'type','Spearman');
    r_z = 0.5*log((1+r_old)./(1-r_old));
    r_recovered = (exp(2 * r_z) - 1) ./ (exp(2 * r_z) + 1);

    S_behv = para_com{1,c}.train_coeff_B;
    S_brain = para_com{1,c}.train_coeff_A;

    out = d_Hungarian(S_brain, S_behv, ref_brain, ref_behv,cv_brain,cv_behv,ref_cv_brain,ref_cv_behv);
   
    Sbehv_perm = out.v;
    Sbrain_perm = out.u;
    Coeff_behv_perm = out.v_c;
    Coeff_brain_perm = out.u_c;

    Stest_brain_perm = para_com{1,c}.out_test_cv_brain * out.matrix;
    Stest_behv_perm = para_com{1,c}.out_test_cv_behv * out.matrix;

    train_ccr_permuted{c} = para_com{1,c}.out_train_ccr(out.pos);
    test_ccr_permuted{c} = diag(corr(Stest_brain_perm, Stest_behv_perm))';

    Sbehv_permuted{c} = Sbehv_perm;
    Sbrain_permuted{c} = Sbrain_perm;
    Coeff_behv_permuted{c} = Coeff_behv_perm;
    Coeff_brain_permuted{c} = Coeff_brain_perm;
    Stest_behv_permuted{c} = Stest_behv_perm;
    Stest_brain_permuted{c} = Stest_brain_perm;
    pos{c} = out.pos;
    sign_change{c} = out.all;

    [r_new,p_new] = corr(Sbehv_perm,z_trsub_behv,'type','Spearman');
    rnew_z = 0.5*log((1+r_new)./(1-r_new));
    r_recovered_new{c} = (exp(2 * rnew_z) - 1) ./ (exp(2 * rnew_z) + 1);

    permuted_para_com{c}.out_train_ccr = train_ccr_permuted{c};
    permuted_para_com{c}.out_test_ccr = test_ccr_permuted{c};
    permuted_para_com{c}.out_train_cv_brain = Sbrain_permuted{c};
    permuted_para_com{c}.out_train_cv_behv = Sbehv_permuted{c};
    permuted_para_com{c}.out_test_cv_brain = Stest_brain_permuted{c};
    permuted_para_com{c}.out_test_cv_behv = Stest_behv_permuted{c};
    permuted_para_com{c}.train_coeff_B = Coeff_behv_permuted{c};
    permuted_para_com{c}.train_coeff_A = Coeff_brain_permuted{c};
end
all_r = cat(3,r_recovered_new{:});
mean_loading = mean(all_r,3);

%% plot fig1 matrix
figure
load('M:/Desktop/colorbar(bbr).mat');
ax1 = axes('Position', [0.16, 0.23, 0.1, 0.6]); 
imagesc(mean_loading(1,:)')
% imagesc(mean_loading(2,:)')
% imagesc(mean_loading(3,:)')
clim([-0.75 0.75]);
colormap(ax1,map);
set(gca, 'XTick', [], 'YTick',[1:3],'YTicklabel',{'VIQ','SA','RRB'},'TickLength',[0,0],'FontSize',14);
%% mean train_cv and test_cv by person
result_train_behv = zeros(length(asd_idx),size(behv,2));
result_train_brain = zeros(length(asd_idx),size(behv,2));
result_test_behv = zeros(length(asd_idx),size(behv,2));
result_test_brain = zeros(length(asd_idx),size(behv,2));
count_per_train = zeros(length(asd_idx),1);
count_per_test = zeros(length(asd_idx),1);

for m = 1:1000
    current_order_train = train_idx1(m,:);
    current_order_test = test_idx1(m,:);
    current_behv = Sbehv_permuted{1,m};
    current_brain = Sbrain_permuted{1,m};
    current_behv_test = Stest_behv_permuted{1,m};
    current_brain_test = Stest_brain_permuted{1,m};
    for n = 1:size(train_idx1,2)
        per_idx = current_order_train(n);
        result_train_behv(per_idx,:) = result_train_behv(per_idx,:) + current_behv(n,:);
        result_train_brain(per_idx,:) = result_train_brain(per_idx,:) + current_brain(n,:);
        count_per_train(per_idx) = count_per_train(per_idx) + 1;
    end
    for n = 1:size(test_idx1,2)
        per_idx = current_order_test(n);
        result_test_behv(per_idx,:) = result_test_behv(per_idx,:) + current_behv_test(n,:);
        result_test_brain(per_idx,:) = result_test_brain(per_idx,:) + current_brain_test(n,:);
        count_per_test(per_idx) = count_per_test(per_idx) + 1;
    end
end
for l = 1:length(asd_idx)
    if count_per_train(l)>0
        result_train_behv(l,:) = result_train_behv(l,:) / count_per_train(l);
        result_train_brain(l,:) = result_train_brain(l,:) / count_per_train(l);
    end
     if count_per_test(l)>0
        result_test_behv(l,:) = result_test_behv(l,:) / count_per_test(l);
        result_test_brain(l,:) = result_test_brain(l,:) / count_per_test(l);
    end
end


%% plot Fig1
mean_train_cv_behv = zscore(result_train_behv);
mean_train_cv_brain = zscore(result_train_brain);
mean_test_cv_behv = zscore(result_test_behv);
mean_test_cv_brain = zscore(result_test_brain);

cv_num = 3;
[r3,p3] = corr(mean_train_cv_brain(:,cv_num),mean_train_cv_behv(:,cv_num));
[P3,S3] = polyfit(mean_train_cv_brain(:,cv_num),mean_train_cv_behv(:,cv_num),1);

set(gcf,'position',[100 100 400 400]);
% s2 = scatter(mean_test_cv_brain(:,cv_num),mean_test_cv_behv(:,cv_num),15,'filled','Marker','o',...
%       'MarkerEdgeColor',[0.953 0.494 0.573],'MarkerFaceColor',[0.953 0.494 0.573]);
% s2 = scatter(mean_test_cv_brain(:,cv_num),mean_test_cv_behv(:,cv_num),15,'filled','Marker','o',...
%       'MarkerEdgeColor',[0.428 0.714 0.831],'MarkerFaceColor',[0.428 0.714 0.831]);
s2 = scatter(mean_test_cv_brain(:,cv_num),mean_test_cv_behv(:,cv_num),15,'filled','Marker','o',...
      'MarkerEdgeColor',[0.404 0.749 0.428],'MarkerFaceColor',[0.404 0.749 0.428]);


hold on
% s1 = scatter(mean_train_cv_brain(:,cv_num),mean_train_cv_behv(:,cv_num),15,'filled','Marker','o',...
%       'MarkerEdgeColor',[0.631 0.294 0.361],'MarkerFaceColor',[0.631 0.294 0.361]);
% s1 = scatter(mean_train_cv_brain(:,cv_num),mean_train_cv_behv(:,cv_num),15,'filled','Marker','o',...
%       'MarkerEdgeColor',[0.192 0.345 0.478],'MarkerFaceColor',[0.192 0.345 0.478]);
s1 = scatter(mean_train_cv_brain(:,cv_num),mean_train_cv_behv(:,cv_num),15,'filled','Marker','o',...
      'MarkerEdgeColor',[0.137 0.357 0.243],'MarkerFaceColor',[0.137 0.357 0.243]);

hold on
xlabel('Connectivity score','FontSize',14)
ylabel('Clinical score','FontSize',14)
title('CV3 - Social Affect related','FontSize',15)
% title('CV2 - Verbal IQ related','FontSize',15)
% title('CV1 - RRB related','FontSize',15)
ylim([-3.2 3.2])
xlim([-3.2 3.2])

