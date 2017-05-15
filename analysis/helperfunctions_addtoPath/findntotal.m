%% function to derive N_total,as described in the MAPseq methods section
function [N_t]=findntotal(matrix,cutoff)

%binarize input
matrix_binary=matrix>=cutoff;
matrix_binary=matrix_binary(sum(matrix_binary,2)~=0,:);

%find counts
N_a=sum(matrix_binary,1);
N_obs=size(matrix_binary,1);

%calculate terms

term5=N_obs-sum(N_a);
term4=sum(prod(nchoosek(N_a,2),2));
term3=sum(prod(nchoosek(N_a,3),2));
term2=sum(prod(nchoosek(N_a,4),2));
term1=sum(prod(nchoosek(N_a,5),2));
term0=prod(N_a);

% solve the polynomial
p = [term5 term4 -1*term3 term2 -1*term1 term0];
r = roots(p);

eq= @(x) term5*x.^5 + term4*x.^4 - term3*x.^3 +term2*x.^2 - term1*x + term0;

%return the largest real root.
N_t=max(r);
% N_t=round(max(r));