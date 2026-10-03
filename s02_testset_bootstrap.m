clc;clear;

datapath = '/user/FC_asd_data/';

load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
behv = [T_asd.t1_viq T_asd.t1_sa_css_all T_asd.t1_rrb_css_all];

%% vector fc
ROI_num = 243;
mask = flipud(tril(ones(ROI_num,ROI_num),-1));
for s = 1:length(asd_idx)
    load([datapath strcat(num2str(T_asd.subjects(s)),'_corr.mat')]);
    lowertril_fc = fc(mask==1); 
    allfc(s,:) = lowertril_fc';
end

Time1 = 1000;
load('/user/sub1000_split.mat');

respath = '/user/RCCA_results_243ROI/';
load([respath 'para_com_in.mat']);

lambda1 = [1,2,3,4,5];
lambda2 = [0.001,0.01];
ft_num = [100:10:400];
test_ccr_orig_all = zeros(Time1, size(behv,2));

%% build 1000 robust ensemble estimates testset CCR with replacement
for i = 1:Time1
    trsub1_fc1 = allfc(train_idx1(i,:),:);
    tesub1_fc1 = allfc(test_idx1(i,:),:);
    trsub1_behv = behv(train_idx1(i,:),:);
    tesub1_behv = behv(test_idx1(i,:),:);

    %% standard fc and behv
    m_trsub1_fc = mean(trsub1_fc1);
    std_trsub1_fc1 = std(trsub1_fc1);
    z_trsub1_fc = (trsub1_fc1 - m_trsub1_fc)./std_trsub1_fc1;
    z_tesub1_fc = (tesub1_fc1 - m_trsub1_fc)./std_trsub1_fc1;

    m_trsub1_behv = mean(trsub1_behv);
    std_trsub1_behv = std(trsub1_behv);
    z_trsub1_behv = (trsub1_behv - m_trsub1_behv)./std_trsub1_behv;
    z_tesub1_behv = (tesub1_behv - m_trsub1_behv)./std_trsub1_behv;

    bootstrap_num = 1000;
    test_ccr_new = nan(bootstrap_num, size(behv,2));
    for b = 1:bootstrap_num
        myresample = randsample(size(z_trsub1_fc,1),size(z_trsub1_fc,1),1);
        r_fc = z_trsub1_fc(myresample,:);
        r_behv = z_trsub1_behv(myresample,:);

        fc_sel = para_com{1,i}.fc_out_idx(1:ft_num(para_com_in{1,i}.max_z));
        [results_out_new] = rcc_matlab(r_fc(:,fc_sel), r_behv, ...
            lambda1(para_com_in{1,i}.max_x), lambda2(para_com_in{1,i}.max_y));

        test_cv_brain_new = z_tesub1_fc(:,fc_sel)*results_out_new.coeff_A;
        test_cv_behv_new = z_tesub1_behv*results_out_new.coeff_B;
        test_ccr_new(b,:) = diag(corr(test_cv_brain_new, test_cv_behv_new))';
    end
     test_ccr_ensemble_cell{i} = test_ccr_new;
     test_ccr_ensemble(i,:) = mean(test_ccr_new,1);
end

%% permuted for null distribution
numPermutations = 1000;
testsub_num = 14;

for p = 1:numPermutations
    shuffled_labels = randperm(length(asd_idx));
    behv_perm = behv(perm_labels, :);

    trsub1_fc1 = allfc(train_idx1(p,:),:);
    tesub1_fc1 = allfc(test_idx1(p,:),:);
    trsub1_behv = behv_perm(train_idx1(p,:),:);
    tesub1_behv = behv_perm(test_idx1(p,:),:);

    %% standard fc and behv
    m_trsub1_fc = mean(trsub1_fc1);
    std_trsub1_fc1 = std(trsub1_fc1);
    z_trsub1_fc = (trsub1_fc1 - m_trsub1_fc)./std_trsub1_fc1;
    z_tesub1_fc = (tesub1_fc1 - m_trsub1_fc)./std_trsub1_fc1;

    m_trsub1_behv = mean(trsub1_behv);
    std_trsub1_behv = std(trsub1_behv);
    Y_permuted_tr = (trsub1_behv - m_trsub1_behv)./std_trsub1_behv;
    Y_permuted_te = (tesub1_behv - m_trsub1_behv)./std_trsub1_behv;

    bootstrap_num = 1000;
    for b = 1:bootstrap_num
        myresample = randsample(size(z_trsub1_fc,1),size(z_trsub1_fc,1),1);
        r_fc = z_trsub1_fc(myresample,:);
        r_behv = Y_permuted_tr(myresample,:);

        fc_sel = para_com{1,p}.fc_out_idx(1:ft_num(para_com_in{1,p}.max_z));
        [results_out_new] = rcc_matlab(r_fc(:,fc_sel), r_behv, ...
            lambda1(para_com_in{1,i}.max_x), lambda2(para_com_in{1,i}.max_y));

        test_cv_brain_new = z_tesub1_fc(:,fc_sel)*results_out_new.coeff_A;
        test_cv_behv_new = Y_permuted_te*results_out_new.coeff_B;
        test_ccr_new(b,:) = diag(corr(test_cv_brain_new, test_cv_behv_new))';
    end
    test_ccr_perm_cell{p} = test_ccr_new;
    test_ccr_perm(p,:) = mean(test_ccr_new,1);
end

save([respath 'mean_test_ccr_perm.mat'],'test_ccr_perm');
save([respath 'mean_test_ccr_ensemble.mat'],'test_ccr_ensemble');


%% corrected variance paired ttest
delta_ccr = test_ccr_ensemble - test_ccr_perm;
p_empirical = nan(1, size(behv,2));
cohens_d = nan(1, size(behv,2));
t_corrected = nan(1, size(behv,2));
p_corrected = nan(1, size(behv,2));

train_size = size(train_idx1, 2);
test_size = size(test_idx1, 2);
corrected_factor = (1 / Time1) + (test_size / train_size);

for cv = 1:size(behv,2)
    p_empirical(cv) = (1 + sum(delta_ccr(:, cv) <= 0)) / (numPermutations + 1);
    cohens_d(cv) = mean(delta_ccr(:, cv)) / std(delta_ccr(:, cv), 0, 1);
    t_corrected(cv) = mean(delta_ccr(:, cv)) / sqrt(corrected_factor * var(delta_ccr(:, cv), 0, 1));
    p_corrected(cv) = 1 - tcdf(t_corrected(cv), Time1 - 1);
end

