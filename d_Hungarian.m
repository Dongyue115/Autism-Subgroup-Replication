function out = d_Hungarian(res_u, res_v, org_u, org_v, real_u, real_v, ref_u, ref_v)

%%%%%%%%%%%%%%%%
% Hungarian method for matched label in 1000 times cross-validation
% Reference Benjamin B. Risk et al.(2014), 'An evaluation of independent component analyses with an application to resting‐state fMRI'
%%%%%%%%%%%%%%%%

    [n_sub,n_ft] = size(res_v);
    
    l2_mat1 = zeros(n_ft,n_ft);
    l2_mat2 = zeros(n_ft,n_ft);

    for j = 1:n_ft
        for i = 1:n_ft
            l2_mat1(i,j) = sum((org_v(:,i)-res_v(:,j)).^2)/n_sub;%row:ref;col:S;
            l2_mat2(i,j) = sum((org_v(:,i)+res_v(:,j)).^2)/n_sub;
        end
    end
    l2_mat1 = sqrt(l2_mat1);
    l2_mat2 = sqrt(l2_mat2);

%% choose the minimum distance for each component
l2_mat = l2_mat1 .* (l2_mat1<=l2_mat2) + l2_mat2 .* (l2_mat2<l2_mat1);
[assign,cost] = munkres(l2_mat);
for q = 1:n_ft
    matching_vector(q) = find(assign(q, :) == 1);
end

sign_change = -1 * (l2_mat2 < l2_mat1) + 1 *(l2_mat1 <= l2_mat2);
sign = sign_change .* assign;

S_behv_perm = real_v * sign';
S_brain_perm = real_u * sign';

coef_behv_perm = res_v * sign';
coef_brain_perm = res_u * sign';

for a = 1:size(sign,2)
    col_sign = sign(:,a);
    sign_simple(:,a) = col_sign(col_sign~=0);
end
allchange = sign_simple .* matching_vector;

res_reorder.u = S_brain_perm;
res_reorder.v = S_behv_perm;
res_reorder.u_c = coef_brain_perm;
res_reorder.v_c = coef_behv_perm;
res_reorder.pos = matching_vector;
res_reorder.all = allchange;
res_reorder.matrix = sign';

%% Output the reordered structure
out = res_reorder;
end