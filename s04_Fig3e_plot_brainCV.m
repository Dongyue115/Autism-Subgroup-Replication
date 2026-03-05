clc;clear;

addpath 'user/munkres';

respath = '/user/RCCA_results_243ROI/';
load([respath 'CV_brain_permuted_1000replicates.mat']);
load('/user/sub1000_split.mat');

avg_CV = cell(1,270);

for i = 1:1000
    cv_brain = Sbrain_permuted{i};
    for j = 1:size(cv_brain,1)
        k = train_idx1(i,j);
        if isempty(avg_CV{k})
            avg_CV{k} = cv_brain(j,:);
        else
            avg_CV{k} = [avg_CV{k};cv_brain(j,:)];
        end
    end
end

meanCols = cellfun(@(x) mean(x,1), avg_CV, 'UniformOutput', false);
avg270_CV = cell2mat(meanCols');

%% plot
mean_cv = avg270_CV;

respath = '/user/RCCA_results_243ROI/';
load([respath 'idx_all_4cluster.mat']);
clusterLabels = zeros(270,1);
clusterLabels(idx_all{3}) = 1;
clusterLabels(idx_all{2}) = 2;
clusterLabels(idx_all{4}) = 3;
clusterLabels(idx_all{1}) = 4;

%% plot 
cv_pick = 1;
GroupColors = [0.753 0.259 0.357;0.282 0.718 0.388;0.420 0.310 0.616;0.369 0.678 0.878]; 
offset = [7.5,5.1,2.7,0];%CV1

for g = 1:4
    idx = find(clusterLabels == g);
    all_data = [];

    for i = idx'
        plot_data = avg_CV{i}(:, cv_pick);%CV1
        [pdf_values, x_values] = ksdensity(plot_data);
        all_data = [all_data;plot_data];

        plot(x_values, pdf_values + offset(g), 'Color', [GroupColors(g,:) 0.5], 'LineWidth', 0.5);
        hold on;
    end
    [pdf_values_all,x_values_all] = ksdensity(all_data); 
    plot(x_values_all, pdf_values_all + offset(g), 'Color', GroupColors(g,:), 'LineWidth', 2.5);
    hold on

    [~, maxIdx] = max(pdf_values_all);
    center_value = x_values_all(maxIdx);

    line([center_value center_value], [offset(g) offset(g) + max(pdf_values_all)], ...
        'Color', GroupColors(g,:), 'LineStyle', '--', 'LineWidth', 1.5);
end

ylim([0 9.7])%CV1
xlim([-4.4 4.4])%CV1

set(gca,'YTick',[0,2.48,2.72,4.88,5.12,7.28,7.52,9.7],'YTickLabels',{'0','2.5','0','2.2','0','2.2','0','2.2'},'fontsize',12);%CV1
set(gcf,'position',[100 100 650 700]); 
xlabel('CV Score','FontSize',16);
ylabel('Density in 4 subgroups','FontSize',16);
title('CV1','FontSize',18)
box off;