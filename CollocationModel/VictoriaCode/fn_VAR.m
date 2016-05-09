%This function estimates the VAR coefficients (A_hat), constant matrix (A0_hat), and
%error variance matrix (SIGMA_hat) for the system Z(t) = A0 + A*Z(t-1) + E(t),
%where E(t) ~ N(0,SIGMA).  Z is the matrix with observations for the
%variables of interest (MxT, with M variables and T observations).  The
%estimates are computed for observations 2 to T (because of the one-period
%lag).

function [A0_hat A_hat SIGMA_hat] = fn_VAR(Z)

M = size(Z,1);  %M is number of variables
T = size(Z,2);  %T is number of observations

X = ones(T-1,M+1);
X(:,2:M+1) = Z(:,1:T-1)';

for i = 1:M
    [B(1:M+1,i),I,R(1:T-1,i)] = regress(Z(i,2:T)',X);
end

A0_hat = B(1,1:M)';   %estimated VAR constant matrix

A_hat = B(2:M+1,:)';    %estimated VAR coefficient matrix

SIGMA_hat = zeros(M,M); %setup for estimated error variance matrix

for i = 1:M
    for j = 1:M
        SIGMA_hat(i,j) = (1/(T-1))*(R(:,i)'*R(:,j));
    end
end %yields estimated error variance matrix
