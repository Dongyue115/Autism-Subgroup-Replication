clc;clear;

datapath = '/user/FC_asd_data/';
savepath = '/user/RCCA_results_243ROI/';

load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
behv = [T_asd.t1_viq T_asd.t1_sa_css_all T_asd.t1_rrb_css_all];

%% vector fc
for s = 1:length(asd_idx)
    load([datapath strcat(num2str(T_asd.subjects(s)),'_corr.mat')]);
    mask = flipud(tril(ones(size(fc)),-1));
    lowertril_fc = fc(mask==1); 
    allfc(s,:) = lowertril_fc';
end

Time1 = 1000;
% [train_idx1,test_idx1] = d_split_samples(Time1,size(T_asd,1),round(size(T_asd,1)*0.05));
load('/user/sub1000_split.mat'); %train_idx1,test_idx1:1000times

for i = 1:Time1
    trsub1_fc1 = allfc(train_idx1(i,:),:);
    tesub1_fc1 = allfc(test_idx1(i,:),:);
    trsub1_behv = behv(train_idx1(i,:),:);
    tesub1_behv = behv(test_idx1(i,:),:);
    %% standard behaviour and fc
    m_trsub1_behv = mean(trsub1_behv);
    std_trsub1_behv = std(trsub1_behv);
    z_trsub1_behv = (trsub1_behv - m_trsub1_behv)./std_trsub1_behv;
    z_tesub1_behv = (tesub1_behv - m_trsub1_behv)./std_trsub1_behv;

    m_trsub1_fc = mean(trsub1_fc1);
    std_trsub1_fc = std(trsub1_fc1);
    z_trsub1_fc = (trsub1_fc1 - m_trsub1_fc)./std_trsub1_fc;
    z_tesub1_fc = (tesub1_fc1 - m_trsub1_fc)./std_trsub1_fc;

    %% choose parameters in RCCA
    Time2 = 30;
    [train_idx2,test_idx2] = d_split_samples(Time2,size(train_idx1,2),round(size(train_idx1,2)*0.05));
    for j = 1:Time2
        trsub2_fc1  = trsub1_fc1(train_idx2(j,:),:);
        tesub2_fc1  = trsub1_fc1(test_idx2(j,:),:);
        trsub2_behv = trsub1_behv(train_idx2(j,:),:);
        tesub2_behv = trsub1_behv(test_idx2(j,:),:);
        %% standard behaviour and fc
        m_trsub2_behv = mean(trsub2_behv);
        std_trsub2_behv = std(trsub2_behv);
        z_trsub2_behv = (trsub2_behv - m_trsub2_behv)./std_trsub2_behv;
        z_tesub2_behv = (tesub2_behv - m_trsub2_behv)./std_trsub2_behv;

        m_trsub2_fc = mean(trsub2_fc1);
        std_trsub2_fc = std(trsub2_fc1);
        z_trsub2_fc = (trsub2_fc1 - m_trsub2_fc)./std_trsub2_fc;
        z_tesub2_fc = (tesub2_fc1 - m_trsub2_fc)./std_trsub2_fc;

        %% feature selection
        Time3 = 100;
        [train_idx3,test_idx3] = d_split_samples(Time3,size(train_idx2,2),round(size(train_idx2,2)*0.05));
        
        fc_r = zeros(size(trsub2_fc1, 2),Time3);
        for k = 1:Time3
            trsub3_fc1 = z_trsub2_fc(train_idx3(k,:),:);
            trsub3_behv = z_trsub2_behv(train_idx3(k,:),:);
            [r_all,p_all] = corr(trsub3_fc1,trsub3_behv,'type','Spearman');
            fc_r(:,k) = mean(abs(r_all),2);
        end
        fc_rmean = mean(fc_r,2);
        [fc_list,fc_idx] = sort(fc_rmean,'descend'); %ranked FC list

        %% parameter combination
        lambda1 = [1,2,3,4,5];
        lambda2 = [0.001,0.01];
        ft_num = [100:10:400];
        
        for l1 = 1:length(lambda1)
            for l2 = 1:length(lambda2)
                for l3 = 1:length(ft_num)
                    %rcc_matlab.m from Buch et al.(2023)
                    [results] = rcc_matlab(z_trsub2_fc(:,fc_idx(1:ft_num(l3))), z_trsub2_behv, lambda1(l1), lambda2(l2));
                    test_ccr(l3,:) = diag(corr(z_tesub2_fc(:,fc_idx(1:ft_num(l3)))*results.coeff_A, z_tesub2_behv*results.coeff_B))';
                end
                test_ccr_all{l1,l2} = test_ccr;
            end
        end
        test_ccr_crosscub{j} = test_ccr_all;
    end
    
    
    for l1 = 1:length(lambda1)
        for l2 = 1:length(lambda2)
            for m = 1:Time2
                ccr1(:,m) = test_ccr_crosscub{m}{l1,l2}(:,1);
            end
            m_ccr1 = median(abs(ccr1),2); %first canonical variate correlation across subsamples
            para_ccr1(l1,l2,:) = m_ccr1;
        end
    end
    [max_absr,max_index] = max(para_ccr1(:));
    [max_x,max_y,max_z] = ind2sub(size(para_ccr1),max_index); %optimal parameter combination

    para_rcca_in{i} = struct('fc_idx',fc_idx, 'para_ccr1',para_ccr1, 'max_absr',max_absr, 'max_x',max_x , 'max_y',max_y , 'max_z',max_z);
  
    %% calculate ccr in test set
    %% feature selection
        Time3_out = 100;
        [train_idx,test_idx] = d_split_samples(Time3_out,size(train_idx1,2),round(size(train_idx1,2)*0.05));
       
         for k = 1:Time3_out
            out_trsub_fc1 = z_trsub1_fc(train_idx(k,:),:);
            out_trsub_behv = z_trsub1_behv(train_idx(k,:),:);
            [r_out,p_out] = corr(out_trsub_fc1,out_trsub_behv,'type','Spearman');
            fc_rout(:,k) = mean(abs(r_out),2);
         end
        fc_routmean = mean(fc_rout,2);
        [fc_out_list,fc_out_idx] = sort(fc_routmean,'descend'); %ranked FC list
        [results_out] = rcc_matlab(z_trsub1_fc(:,fc_out_idx(1:ft_num(max_z))), z_trsub1_behv, lambda1(max_x), lambda2(max_y));

        out_test_cv_brain = z_tesub1_fc(:,fc_out_idx(1:ft_num(max_z)))*results_out.coeff_A;
        out_test_cv_behv = z_tesub1_behv*results_out.coeff_B;
        out_test_ccr = diag(corr(out_test_cv_brain, out_test_cv_behv))';

        out_train_cv_brain = results_out.variate_U;
        out_train_cv_behv = results_out.variate_V;
        out_train_ccr = results_out.r;

        train_coeff_A = results_out.coeff_A;
        train_coeff_B = results_out.coeff_B;

        para_rcca_out{i} = struct('out_test_cv_brain',out_test_cv_brain, 'out_test_cv_behv',out_test_cv_behv, 'out_test_ccr',out_test_ccr,...
            'out_train_cv_brain',out_train_cv_brain, 'out_train_cv_behv',out_train_cv_behv, 'out_train_ccr',out_train_ccr,...
            'fc_out_idx',fc_out_idx,'train_coeff_A',train_coeff_A,'train_coeff_B',train_coeff_B);
        
        save([savepath 'in/' strcat('para_rcca_in_',num2str(i,'%04d'),'.mat')],'para_rcca_in');
        save([savepath 'out/' strcat('para_rcca_out_',num2str(i,'%04d'),'.mat')],'para_rcca_out');
end


savepath_in = [savepath 'in/'];
savepath_out = [savepath 'out/'];

filename_in = {dir(savepath_in).name};
filename_in = filename_in(3:end);
filename_out = {dir(savepath_out).name};
filename_out = filename_out(3:end);

for i = 1:length(filename_in)
    para_mat_in = importdata([savepath_in filename_in{i}]);
    para_com_in{i} = para_mat_in{i};
    para_mat_out = importdata([savepath_out filename_out{i}]);
    para_com_out{i} = para_mat_out{i};
end

save([savepath 'para_com_in.mat'],'para_com_in');
save([savepath 'para_com_out.mat'],'para_com_out');