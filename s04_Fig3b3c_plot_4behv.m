clc;clear;

addpath '/user/munkres';

respath = '/user/RCCA_results_243ROI/';
load([respath 'para_com_out.mat']);
load('/user/sub1000_split.mat');

num_clusters = 4;

for i = 1:length(para_com_out)
    cv_brain = para_com_out{1,i}.out_train_cv_brain;
    cosine_distance = pdist(zscore(cv_brain), 'cosine');
    Z = linkage(cosine_distance, 'average');

    clusterLabels = cluster(Z, 'maxclust', num_clusters);
    all_labels(:,i) = clusterLabels;
end

%% top up to all 270subs label
for i = 1:1000
    current_labels = all_labels(:, i);
    current_labels_allsub = zeros(270,1);
    for j = 1:length(clusterLabels)
        current_labels_allsub(train_idx1(i,j))  = current_labels(j);
    end
    allsubsamples_label(i,:) = current_labels_allsub';
end

%% match labels
allsub_aligned_label = zeros(1000,270);

base_labels = allsubsamples_label(1,:); %ref_label
allsub_aligned_label(1,:) = base_labels;

for i = 2:size(allsubsamples_label,1)
    label_for_align = allsubsamples_label(i,:);
    
    label1 = base_labels(base_labels~=0 & label_for_align~=0);
    label2 = label_for_align(base_labels~=0 & label_for_align~=0);
    
    confusion_matrix = zeros(num_clusters, num_clusters);
    for j = 1:num_clusters
        for k = 1:num_clusters
            confusion_matrix(j, k) = sum((label1 == j) & (label2 == k));
        end
    end
    [assignment,~] = munkres(-confusion_matrix);
    for a = 1:num_clusters
        assign_order(a) = find(assignment(a, :));
    end

    new_labels2 = zeros(size(label2));
    new_labels1 = zeros(size(label_for_align));
    for p = 1:num_clusters
        new_labels2(label2 == p) = assign_order(p);
        new_labels1(label_for_align == p) = assign_order(p);
    end
    allsub_aligned_label(i,:) = new_labels1;
end

%% evaluate the cluster assignment
sub_cluster_count = zeros(num_clusters,size(allsub_aligned_label,2));
for i = 1:size(allsub_aligned_label,2)
    for j = 1:size(allsub_aligned_label,1)
        if allsub_aligned_label(j,i)==1
            sub_cluster_count(1,i) = sub_cluster_count(1,i)+1;
        elseif allsub_aligned_label(j,i)==2
            sub_cluster_count(2,i) = sub_cluster_count(2,i)+1;
        elseif allsub_aligned_label(j,i)==3
            sub_cluster_count(3,i) = sub_cluster_count(3,i)+1;
        elseif allsub_aligned_label(j,i)==4
            sub_cluster_count(4,i) = sub_cluster_count(4,i)+1; 
        end
    end
end

for i = 1:size(sub_cluster_count,2)
    assign_pro(:,i) = sub_cluster_count(:,i)/sum(sub_cluster_count(:,i));
end

idx_all = cell(1,num_clusters);
count = [];
for i = 1:size(assign_pro,2)
    idx = find(assign_pro(:,i)==max(assign_pro(:,i)));
    if length(idx)~=1
        idx = idx(1); %random assignment
        count = [count;i];
    end
        if idx==1
            idx_all{1} = [idx_all{1};i];
        elseif idx==2
            idx_all{2} = [idx_all{2};i];
        elseif idx==3
            idx_all{3} = [idx_all{3};i];
        elseif idx==4
            idx_all{4} = [idx_all{4};i]; 
        end
end
save([respath 'idx_all_4cluster.mat'],'idx_all');

%% plot
addpath(genpath('github_repo'),'-end')
RGB = cbrewer2('seq', 'OrRd', 200, 'linear');

resort_idx = [idx_all{3};idx_all{2};idx_all{4};idx_all{1}];%alter order based Buch's results
new_assign_pro = [assign_pro(3,:);assign_pro(2,:);assign_pro(4,:);assign_pro(1,:)];
resort_pro = new_assign_pro(:,resort_idx);
imagesc(resort_pro)
colormap(RGB); 
caxis([0 0.8]);
colorbar('Ticks',0:0.2:0.8);

hold on
annotation('textbox', [0.16,0.75, 0.1, 0.1], ...
    'String',num2str(length(idx_all{3})), 'EdgeColor', 'none','HorizontalAlignment', 'center', 'Color', 'black','FontSize',14);
annotation('textbox', [0.36,0.55, 0.1, 0.1], ...
    'String',num2str(length(idx_all{2})), 'EdgeColor', 'none','HorizontalAlignment', 'center', 'Color', 'black','FontSize',14);
annotation('textbox', [0.56,0.35, 0.1, 0.1], ...
    'String',num2str(length(idx_all{4})), 'EdgeColor', 'none','HorizontalAlignment', 'center', 'Color', 'black','FontSize',14);
annotation('textbox', [0.73,0.15, 0.1, 0.1], ...
    'String',num2str(length(idx_all{1})), 'EdgeColor', 'none','HorizontalAlignment', 'center', 'Color', 'black','FontSize',14);
x_lines = [length(idx_all{3})+0.5,length(idx_all{3})+length(idx_all{2})+0.5,length(idx_all{3})+length(idx_all{2})+length(idx_all{4})+0.5];    

y_range = [0,5];
for i = 1:length(x_lines)
    line([x_lines(i) x_lines(i)], y_range, 'Color', 'k', 'LineWidth', 1.2);
end

title('Label proportion in 1000 training replicates')
xlabel('Subjects','FontSize',14)
set(gca, 'xtick', [],'ytick',[1,2,3,4],'yticklabel',{'S1','S2','S3','S4'},'Fontsize',14);
set(gcf,'position',[100 100 800 300]);  

%% subtype clinical
load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
behv = [T_asd.t1_sa_css_all T_asd.rrb_css T_asd.t1_viq T_asd.t1_css_total_all];
z_behv = zscore(behv);

clusterLabels = zeros(size(assign_pro,2),1);
clusterLabels(idx_all{3}) = 1;
clusterLabels(idx_all{2}) = 2;
clusterLabels(idx_all{4}) = 3;
clusterLabels(idx_all{1}) = 4;

%% figure
symptom_labels = {'Social Affect','RRB','Verbal IQ','Total Severity'};
GroupColors = [0.753 0.259 0.357;0.282 0.718 0.388;0.420 0.310 0.616;0.369 0.678 0.878]; 
subtype = unique(clusterLabels);

t = tiledlayout(1, length(symptom_labels), 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
for i = 1:size(z_behv,2)
    nexttile;
    symptomData = z_behv(:,i);

    for j = 1:length(subtype)
        yPos(j) = median(symptomData(clusterLabels == j));
        bar(j,yPos(j),'FaceColor',GroupColors(j, :));
        hold on
    end

    box = boxplot(symptomData, clusterLabels,'Colors', 'k', 'Symbol', 'o', 'Widths', 0.25);
    h = findobj(box, 'Tag', 'Box'); 
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), 'k', 'FaceAlpha', 1); 
    end

    for q = 1:size(box,2)
        medianValue = median(symptomData(clusterLabels == q));
        set(box(3,q), 'Visible', 'off'); 
        set(box(4,q), 'Visible', 'off');
        verticalLines = findobj(box(:,q), 'Type', 'line', 'LineStyle', '--');
        set(verticalLines, 'LineStyle', '-');
        set(box(6,q), 'Visible', 'off');

        hold on;
        plot(q, medianValue, 'wo', 'MarkerSize', 2, 'LineWidth', 1.2);
        outliers = findobj(box(:,q), 'Tag', 'Outliers'); 
        set(outliers, 'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'k'); 
    end

    line([-0.2,5],[0,0],'linestyle','--','Color','black');
    % ylim([-2.8 2.8])
    ylim([-2.5 2.5])
    ylabel('Z')
    xlabel('Subgroups')
    title(symptom_labels{i}) 
    set(gca,'Fontsize',14,'box','off')
end
set(gcf,'unit','centimeters','position',[3 2 30 10])
