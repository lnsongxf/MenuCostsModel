function [Pr_mat,Pr_mat_key,zbar] = fn_var_to_markov(A0,A1,A2,SIGMA,N,random_draws,method)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GENERALIZED MARKOV APPROXIMATIONS TO VAR PROCESSES
%
% This function converts a VAR to a discretized Markov process,
% generalizing the approach in Tauchen (1986) by allowing for more general
% var./cov. structure in the error term.  The required multivariate normal
% probabilities are calculated using a Monte Carlo-type technique
% implemented in the function qscmvnv.m, developed by and available on the
% website of Alan Genz: http://www.math.wsu.edu/faculty/genz/homepage.
%
% Original VAR: A0*Z(t) = A1 + A2*Z(t-1) + e(t), e(t) ~ N(0,SIGMA)
%
% INPUTS:
% 1 - A0, A1, A2 are the VAR coefficients, as indicated above, with
%     A0 assumed non-singular
% 2 - N is n x 1, where n = # of vars, N(i) = # grid points for ith var.
% 3 - SIGMA is the arbitrary positive semi-definite error var./cov. matrix
% 4 - random_draws is the number of random draws used in the required
%     Monte Carlo-type integration of the multivariate normal
% 5 - method switch determines the grid selection method
%       - method = 1 uses a uniformly spaced grid covering a fixed number
%         of std. dev. of the relevant component variables.  This is the
%         grid spacing strategy proposed in Tauchen (1986).
%       - method = 2 selects grid points based on approximately equal
%         weighting from the UNIVARIATE normal cdf.  This method is adapted
%         from code written by Jonathan Willis.  (Note that method = 2
%         requires the use of the MATLAB statistics toolbox.)
%
% OUTPUTS:
% 1 - Pr_mat is the Prod(N) x Prod(N) computed transition probability matrix
%     for the discretized Markov chain
% 2 - Pr_mat_key is n x Prod(N) matrix s.t. if Z* is the discretized Markov
%     approximation to the VAR Z, then Z*(state i) = Pr_mat_key(:,i)
% 3 - zbar is the n x max(N) matrix s.t. zbar(i,1:N(i)) is the univariate
%     grid for the ith component of Z*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(N,1);  %number of variables in VAR

%compute reduced form parameters & steady-state mean
A1bar = inv(A0)*A1;
A2bar = inv(A0)*A2;
SIGMAbar = inv(A0)*SIGMA*(inv(A0)');

sstate_mean = inv(eye(n)-A2bar)*A1bar;

m = 2;  %number std deviations of the VAR covered by grid

%iterate to obtain var./cov. structure of the PROCESS (not error term)
SIGMAprocess = SIGMAbar;
SIGMAprocess_last = SIGMAprocess;
dif = 1;
while dif>0.00000001;
    SIGMAprocess = A2bar*SIGMAprocess_last*(A2bar') + SIGMAbar;
    dif = max(max(SIGMAprocess-SIGMAprocess_last));
    SIGMAprocess_last = SIGMAprocess;
end;

%This block equally spaces grid points bounded by m*(std.deviation of
%process) on either side of the unconditional mean of the process.  Any
%more sophisticated spacing of the grid points could be implemented by
%changing the definition of zbar below.
zbar = zeros(n,max(N)); 
grid_stdev = diag(SIGMAprocess).^0.5;
if method==1;
    grid_increment = zeros(n,1);
    for i = 1:n;
         grid_increment(i) = 2*m*grid_stdev(i)/(N(i)-1);
         zbar(i,1) = -m*grid_stdev(i) + sstate_mean(i);
        for j = 1:N(i)-1;
            zbar(i,j+1) = zbar(i,j) + grid_increment(i);
        end;
    end;
elseif method==2;
    d = zeros(n,max(N));
    b = -4:.005:4;
    c = normcdf(b,0,1);
    for i = 1:n;
        a = (1/(2*N(i))):(1/N(i)):1;
        for j = 1:N(i);
            [d1,d(i,j)] = min((a(j)-c).^2);
        end;
        zbar(i,1:N(i)) = grid_stdev(i)*b(d(i,:))+sstate_mean(i);
    end;
end;

%compute key matrix & pos matrix
Pr_mat_key = zeros(length(N),prod(N));
Pr_mat_key_pos = zeros(length(N),prod(N));
Pr_mat_key(length(N),:) = repmat(zbar(length(N),1:N(length(N))),[1 prod(N)/N(length(N))]);
Pr_mat_key_pos(length(N),:) = repmat(1:N(length(N)),[1 prod(N)/N(length(N))]);
for i=length(N)-1:-1:1;
    Pr_mat_key(i,:) = repmat(kron(zbar(i,1:N(i)),ones(1,prod(N(i+1:length(N))))),[1 prod(N)/prod(N(i:length(N)))]);
    Pr_mat_key_pos(i,:) = repmat(kron(1:N(i),ones(1,prod(N(i+1:length(N))))),[1 prod(N)/prod(N(i:length(N)))]);
end;

nstate = prod(N);
Pr_mat_intervals = zeros(n,nstate,2);   %this will store the unadjusted limits of integration for each variable in each state, for input into the Genz code
if method==1;
    for i = 1:nstate;  %number of states
        for j = 1:n;    %number of variables
            if Pr_mat_key_pos(j,i)==1;
                Pr_mat_intervals(j,i,1) = -inf;
                Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + (grid_increment(j)/2);
            elseif Pr_mat_key_pos(j,i)==N(j);
                Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - (grid_increment(j)/2);
                Pr_mat_intervals(j,i,2) = inf;
            else
                Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - (grid_increment(j)/2);
                Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + (grid_increment(j)/2);
            end;
        end;
    end;
elseif method==2;
    for i = 1:nstate;  %number of states
        for j = 1:n;    %number of variables
            if Pr_mat_key_pos(j,i)==1;
                Pr_mat_intervals(j,i,1) = -inf;
                Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + (zbar(j,Pr_mat_key_pos(j,i)+1)-zbar(j,Pr_mat_key_pos(j,i)))/2;
            elseif Pr_mat_key_pos(j,i)==N(j);
                Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - (zbar(j,Pr_mat_key_pos(j,i))-zbar(j,Pr_mat_key_pos(j,i)-1))/2;
                Pr_mat_intervals(j,i,2) = inf;
            else
                Pr_mat_intervals(j,i,1) = zbar(j,Pr_mat_key_pos(j,i)) - (zbar(j,Pr_mat_key_pos(j,i))-zbar(j,Pr_mat_key_pos(j,i)-1))/2;
                Pr_mat_intervals(j,i,2) = zbar(j,Pr_mat_key_pos(j,i)) + (zbar(j,Pr_mat_key_pos(j,i)+1)-zbar(j,Pr_mat_key_pos(j,i)))/2;
            end;
        end;
    end;
end;

error_est = zeros(nstate,nstate);
Pr_mat_intervals_adjusted = zeros(n,nstate,2);
Pr_mat = zeros(nstate,nstate);
for i = 1:nstate; %rows of Pr_mat
    Pr_mat_intervals_adjusted(:,:,1) = Pr_mat_intervals(:,:,1) - repmat((A1bar + A2bar*Pr_mat_key(:,i)),1,nstate);
    Pr_mat_intervals_adjusted(:,:,2) = Pr_mat_intervals(:,:,2) - repmat((A1bar + A2bar*Pr_mat_key(:,i)),1,nstate);
    for j = 1:nstate;   %columns of Pr_mat
        %Pr_mat(i,j) = P(state j|state i)
        [Pr_mat(i,j), error_est(i,j)] = qscmvnv(random_draws,SIGMAbar,Pr_mat_intervals_adjusted(:,j,1),eye(n),Pr_mat_intervals_adjusted(:,j,2));
    end;
end;

%rounding error adjustment
round_sum = sum(Pr_mat,2);
for i = 1:size(Pr_mat,2);
    Pr_mat(i,:) = Pr_mat(i,:)/round_sum(i);
end