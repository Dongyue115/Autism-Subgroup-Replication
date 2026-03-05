function [all_train_splits,all_test_splits] = d_split_samples(split_nums,allsub_size,test_size)

    train_size = allsub_size - test_size;

    all_test_splits = zeros(split_nums,test_size);
    all_train_splits = zeros(split_nums,train_size);

    test_splits_set = containers.Map;

    for i = 1:split_nums
        %% no repeated series
        is_unique_split = false;

        while ~is_unique_split
            shuffled_labels = randperm(allsub_size);
            test_labels = shuffled_labels(1:test_size);
            train_labels = shuffled_labels(test_size+1:end);
            
            test_labels_str = mat2str(sort(test_labels));
            
            if ~isKey(test_splits_set,test_labels_str)
                test_splits_set(test_labels_str) = true;

                is_unique_split = true;

                all_test_splits(i,:) = test_labels;
                all_train_splits(i,:) = train_labels;
            end

        end
    end

end