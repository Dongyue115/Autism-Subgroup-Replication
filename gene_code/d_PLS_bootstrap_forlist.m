% function [geneindex1 PLS1 Z1] = d_PLS_bootstrap(response_var_file, predictor_var_file,ROInum,genes_name)
function [geneindex2 PLS2 Z2] = d_PLS_bootstrap(response_var_file, predictor_var_file,ROInum,genes_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS bootstrap function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% response_var_file ------ full path to the PLS_MRI_response_vars.csv file
%%%                           that is created by the NSPN_CorticalMyelination 
%%%                           wrapper script
%%% predictor_var_file ----- full path to the PLS_gene_predictor_vars.csv file
%%%                           that is provided as raw data
%%% output_dir ------------- where to save the PLS_geneWeights and PLS_ROIscores 
%%%                           files (for PLS1 and PLS2 separately)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Running PLS')

%unwrap and tidy MRI response variable names and store the response variables in matrix Y
MRIdata = response_var_file;

%import predictor variables
GENEdata = predictor_var_file;
geneindex = 1:size(GENEdata,2);

%number of bootstrap iterations
bootnum = 1000;

%DO PLS in 2 dimensions (with 2 components) 
X = zscore(GENEdata);
Y = zscore(MRIdata);
dim = 10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,dim);

%store regions IDs and weights in descending order of weight for both
%components
[R1,p1] = corr([XS(:,1)],MRIdata);
[R2,p2] = corr([XS(:,2)],MRIdata);

%align PLS components with desired direction%
if R1(1,1)<0
    stats.W(:,1) = -1*stats.W(:,1);
    XS(:,1) = -1*XS(:,1);
end
if R2(1,1)<0
    stats.W(:,2) = -1*stats.W(:,2);
    XS(:,2) = -1*XS(:,2);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids = genes_name(x1);
geneindex1 = geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids = genes_name(x2);
geneindex2 = geneindex(x2);

region_ind = 1:ROInum;
PLS1_ROIscores = zeros(ROInum,1);
PLS2_ROIscores = zeros(ROInum,1);

PLS1_score = XS(:,1);
PLS2_score = XS(:,2);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights = [];
PLS2weights = [];

%start bootstrap
disp('  Bootstrapping - could take a while')
for i = 1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:) = myresample; %store resampling out of interest
    Xr = X(myresample,:); % define X for resampled regions
    Yr = Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Xr,Yr,dim); %perform PLS for resampled data
    
    temp = stats.W(:,1);%extract PLS1 weights
    newW = temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW = -1*newW;
    end
    PLS1weights = [PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp = stats.W(:,2);%extract PLS2 weights
    newW = temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW = -1*newW;
    end
    PLS2weights = [PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw = std(PLS1weights');
PLS2sw = std(PLS2weights');

%get bootstrap weights
temp1 = PLS1w./PLS1sw';
temp2 = PLS2w./PLS2sw';

nan_idx = find(isnan(temp1));%when brainspan, this situation will come,maybe cuz the same row inX
candidates = (-0.05:0.001:0.05)';
random_indices = randperm(length(candidates), length(nan_idx));
temp1(nan_idx) = candidates(random_indices);

%order bootstrap weights (Z) and names of regions
[Z1 ind1] = sort(temp1,'descend');
PLS1 = PLS1ids(ind1);
geneindex1 = geneindex1(ind1);
[Z2 ind2] = sort(temp2,'descend');
PLS2 = PLS2ids(ind2);
geneindex2 = geneindex2(ind2);
