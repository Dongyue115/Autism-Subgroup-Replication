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
    for b = 1:bootstrap_num
        myresample = randsample(size(z_trsub1_fc,1),size(z_trsub1_fc,1),1);
        r_fc = z_trsub1_fc(myresample,:);
        r_behv = z_trsub1_behv(myresample,:);

        fc_sel = para_com_in{1,i}.fc_idx(1:ft_num(para_com_in{1,i}.max_z));
        [results_out_new] = rcc_matlab(r_fc(:,fc_sel), r_behv, ...
            lambda1(para_com_in{1,i}.max_x), lambda2(para_com_in{1,i}.max_y));

        test_cv_brain_new = z_tesub1_fc(:,fc_sel)*results_out_new.coeff_A;
        test_cv_behv_new = z_tesub1_behv*results_out_new.coeff_B;
        test_ccr_new(b,:) = diag(corr(test_cv_brain_new, test_cv_behv_new))';
    end
     test_ccr_ensemble_cell{i} = test_ccr_new;
end

%% permuted for null distribution
numPermutations = 1000;
testsub_num = 14;

for p = 1:numPermutations
    shuffled_labels = randperm(length(asd_idx));
    test_labels = shuffled_labels(1:testsub_num);
    train_labels = shuffled_labels(testsub_num+1:end);

    trsub1_behv = behv(train_labels,:);
    tesub1_behv = behv(test_labels,:);
    trsub1_fc1 = allfc(train_idx1(p,:),:);
    tesub1_fc1 = allfc(test_idx1(p,:),:);

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

        fc_sel = para_com_in{1,i}.fc_idx(1:ft_num(para_com_in{1,i}.max_z));
        [results_out_new] = rcc_matlab(r_fc(:,fc_sel), r_behv, ...
            lambda1(para_com_in{1,i}.max_x), lambda2(para_com_in{1,i}.max_y));

        test_cv_brain_new = z_tesub1_fc(:,fc_sel)*results_out_new.coeff_A;
        test_cv_behv_new = Y_permuted_te*results_out_new.coeff_B;
        test_ccr_new(b,:) = diag(corr(test_cv_brain_new, test_cv_behv_new))';
    end
    test_ccr_perm_cell{p} = test_ccr_new;
end

save([respath 'test_ccr_perm_cell.mat'],'test_ccr_perm_cell');
save([respath 'test_ccr_ensemble_cell.mat'],'test_ccr_ensemble_cell');
