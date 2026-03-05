function [results]=rcc_matlab(X,Y,lambda1,lambda2)
%%%%%%%%%%%%%%%%
% RCC for Matlab with L2-norm translated from R MixOmics
% From Buch et al. (2023)
%%%%%%%%%%%%%%%%

[nrowX,ncolX]=size(X);
[nrowY,ncolY]=size(Y);

% Covariance and variance (with L2-norm penalties lambda1 and lambda2)
Cxx=cov(X)+ diag(lambda1*ones(ncolX,1));
Cyy=cov(Y)+ diag(lambda2*ones(ncolY,1));

% Mean center
Cxy=( (X-mean(X))'*(Y-mean(Y)) )/(nrowX-1);

Cxx_dim=size(Cxx);
Cyy_dim=size(Cyy);
p=Cxx_dim(1);
q=Cyy_dim(1);

Cxx_new=(Cxx+Cxx')/2;
Cyy_new=(Cyy+Cyy')/2;

% Cholesky factorization
Cxxfac=chol(Cxx_new);
Cyyfac=chol(Cyy_new);

% Matrix inverse
Cxxfacinv=inv(Cxxfac);
Cyyfacinv=inv(Cyyfac);
covmat=Cxxfacinv'*Cxy*Cyyfacinv;

% Singular value decomposition,
[a,d,b] = svd(covmat,'econ');

% canonical coefficients
A=Cxxfacinv * a;
B=Cyyfacinv * b;
r=diag(d)';

% standardize weights
X_aux = X-mean(X);
Y_aux = Y-mean(Y);
X_aux(isnan(X_aux)) = 0;
Y_aux(isnan(Y_aux)) = 0;

% Canonical Variates
U = X_aux*A;
V = Y_aux*B;

results=struct('r',r,'variate_U',U, 'variate_V',V,'coeff_A',A, 'coeff_B',B);

