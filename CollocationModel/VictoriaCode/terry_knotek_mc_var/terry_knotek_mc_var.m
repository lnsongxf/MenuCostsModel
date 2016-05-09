%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates results for: "Markov-Chain Approximations of Vector
% Autoregressions: Application of General Multivariate-Normal Integration
% Techniques"
%
% Citation:
% Terry, Stephen J. and Edward S. Knotek II (2011) "Markov-Chain
% Approximations of Vector Autoregressions: Application of General
% Multivariate-Normal Integration Techniques" Economics Letters 110(1): 4-6.
%
% Working paper version: 
% Terry, Stephen J. and Edward S. Knotek II (2008) "Markov-Chain
% Approximations of Vector Autoregressions: Application of General
% Multivariate-Normal Integration Techniques" Federal Reserve Bank of
% Kansas City Research Working Paper 08-02.
% Available online at:
% http://www.kansascityfed.org/Publicat/Reswkpap/rwpmain.htm
%
% Calls files:
% * fn_var_to_markov.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;      % closes open windows/figures
clear all;      % clears old variables from memory

%seed normal random number generator
seed_no = 16349618410;
randn('seed',seed_no);

%seed uniform random number generator
rand('seed',1111111);

sims = 1000;  %number of simulations
T = 100; %length of simulated series
disp(['sims = ', num2str(sims)]);
disp(['T = ', num2str(T)]);

random_draws = 1000;
disp(['random_draws=',num2str(random_draws)]);

%In each example coeffs will be:
A1bar = [-.5; 0.9; 0.6];
A2bar = [.25 .1 0.5; -.5 .09 -.75; 0.6 0 .15];

disp('In each example, the actual coefficients are:');
disp('A1bar=');
disp(num2str(A1bar));
disp('A2bar=');
disp(num2str(A2bar));
disp('Eigenvalue of A2bar');
disp(num2str(eig(A2bar)));

%grid points will be
N = [5; 5; 5];

%general setup/storage
M = size(A1bar,1);  %number of variables in simulated VARs
a = 1;  %2= a+1 =Markov approx. vs. direct sim
Pr_mat = zeros(prod(N),prod(N));    %used to store the probability matrices from each of the alternative functions
Pr_mat_key = zeros(M,prod(N));    %used to store the keys from each of the alternative functions
zbar = zeros(M,max(N));   %univariate grid key

%Non-diagonal and non-singular error covariance
tic
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Non-diagonal and non-singular error covariance')
SIGMAbar = [.4 .18 .3; .18 .2 .1; .3 .1 .7];
disp('SIGMAbar=');
disp(num2str(SIGMAbar));
disp('Eigenvalue of SIGMAbar');
disp(num2str(eig(SIGMAbar)));

[Pr_mat(:,:),Pr_mat_key(:,:),zbar(:,:)] = fn_var_to_markov(eye(M),A1bar,A2bar,SIGMAbar,N,random_draws,1);
toc

%obtain uniform random draws for use in Markov simulation
U = zeros(T,sims);
U(2:T,:) = unifrnd(0,1,T-1,sims);

%matrix for storage of simulated series, a+1 dimension because it will also
%store directly simulated series
Z = zeros(M,T,sims,a+1);

%initial state is middle state
initial_state = ((prod(N)+1)/2);    %=steady-state mean

%matrix for use in simulation of the Markov chains
Pr_sum = zeros(prod(N),prod(N));

%actually simulate the Markov chains
for j = 1:a
    Pr_sum(:,1) = Pr_mat(:,1);
    for i = 2:prod(N);
        Pr_sum(:,i) = Pr_sum(:,i-1) + Pr_mat(:,i);
    end

    old_state = initial_state;
    
    for i = 1:sims
        Z(:,1,i,j) = Pr_mat_key(:,initial_state);
    end
    
    for i = 1:sims
        for t = 2:T
            new_state=find(U(t,i)<=Pr_sum(old_state,:),1);
            Z(:,t,i,j) = Pr_mat_key(:,new_state);
            old_state = new_state;
        end
    end
end

%directly simulate the VAR process

%draw normal errors (easy case)
E = zeros(M,T,sims);
for sim_count = 1:sims
    E(:,:,sim_count) = mvnrnd([0; 0; 0],SIGMAbar,T)';
end;

%initialize at steady state
for j = 1:sims
    Z(:,1,j,a+1) = (inv(eye(M)-A2bar)*A1bar);
end

%obtain the rest of the VAR process
for t=2:T;
    for i=1:sims;
        Z(:,t,i,a+1) = A1bar(:,:) + A2bar*Z(:,t-1,i,a+1) + E(:,t,i);
    end;
end;

%now we will estimate A, A0, and SIGMAZ from the simulations in the large Z
%matrix, store the estimates, mean estimates, and st dev of estimates in
%the matrices below

A1bar_hat = zeros(M,sims,a+1);
A1bar_mean = zeros(M,1,a+1);

A2bar_hat = zeros(M,M,sims,a+1);
A2bar_mean = zeros(M,M,a+1);

SIGMAbar_hat = zeros(M,M,sims,a+1);
SIGMAbar_mean = zeros(M,M,a+1);

%obtain estimates from Stephen's function fn_VAR
for j = 1:a+1
    for i = 1:sims
        [A1bar_hat(:,i,j), A2bar_hat(:,:,i,j), SIGMAbar_hat(:,:,i,j)] = fn_VAR(Z(:,:,i,j));
    end
end


for j = 1:a+1
    for k = 1:M
        A1bar_mean(k,1,j) = mean(A1bar_hat(k,:,j));
        for l = 1:M
            A2bar_mean(k,l,j) = mean(A2bar_hat(k,l,:,j));
            SIGMAbar_mean(k,l,j) = mean(SIGMAbar_hat(k,l,:,j));
        end
    end
end

disp('   ');
disp('Vector Process Comparison:');
disp('   ');
disp('Actual constant matrix');
disp(A1bar);
disp('Mean estimates for constant matrix');
disp('    Markov    Direct Simulation');
disp([A1bar_mean(:,1,1),A1bar_mean(:,1,2)]);
disp('Actual coefficient matrix');
disp(A2bar);
disp('Mean estimates for coefficient matrix');
disp('    Markov                        Direct Simulation');
disp([A2bar_mean(:,:,1),A2bar_mean(:,:,2)]);
disp('Actual error variance matrix');
disp(SIGMAbar);
disp('Mean estimates for error variance matrix');
disp('    Markov                        Direct Simulation');
disp([SIGMAbar_mean(:,:,1),SIGMAbar_mean(:,:,2)]);

















































%seed normal random number generator
seed_no = 16349618410;
randn('seed',seed_no);

%seed uniform random number generator
rand('seed',1111111);

%general setup/storage
M = size(A1bar,1);  %number of variables in simulated VARs
a = 1;  %2= a+1 =Markov approx. vs. direct sim
Pr_mat = zeros(prod(N),prod(N));    %used to store the probability matrices from each of the alternative functions
Pr_mat_key = zeros(M,prod(N));    %used to store the keys from each of the alternative functions
zbar = zeros(M,max(N));   %univariate grid key

%Non-diagonal and singular covariance
tic
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Non-diagonal and singular covariance')
A0 = [1 0 0; 1 -1 -1; 0 0 1];
SIGMAbar = inv(A0)*[.1^2 0 0; 0 0 0; 0 0 .3^2]*(inv(A0)');
disp('A0=');
disp(num2str(A0));
disp('SIGMA=');
disp(num2str([.1^2 0 0; 0 0 0; 0 0 .3^2]));
disp('SIGMAbar=');
disp(num2str(SIGMAbar));
disp('Eigenvalue of SIGMAbar (check for stability)');
disp(num2str(eig(SIGMAbar)));

[Pr_mat(:,:),Pr_mat_key(:,:),zbar(:,:)] = fn_var_to_markov(eye(M),A1bar,A2bar,SIGMAbar,N,random_draws,1);
toc

%obtain uniform random draws for use in Markov simulation
U = zeros(T,sims);
U(2:T,:) = unifrnd(0,1,T-1,sims);

%matrix for storage of simulated series, a+1 dimension because it will also
%store directly simulated series
Z = zeros(M,T,sims,a+1);

%initial state is middle state
initial_state = ((prod(N)+1)/2);    %=steady-state mean

%matrix for use in simulation of the Markov chains
Pr_sum = zeros(prod(N),prod(N));

%actually simulate the Markov chains
for j = 1:a
    Pr_sum(:,1) = Pr_mat(:,1);
    for i = 2:prod(N);
        Pr_sum(:,i) = Pr_sum(:,i-1) + Pr_mat(:,i);
    end

    old_state = initial_state;
    
    for i = 1:sims
        Z(:,1,i,j) = Pr_mat_key(:,initial_state);
    end
    
    for i = 1:sims
        for t = 2:T
            new_state=find(U(t,i)<=Pr_sum(old_state,:),1);
            Z(:,t,i,j) = Pr_mat_key(:,new_state);
            old_state = new_state;
        end
    end
end

%directly simulate the VAR process

%draw normal errors (easy case)
E = zeros(M,T,sims);
for sim_count = 1:sims
    E(:,:,sim_count) = mvnrnd([0; 0; 0],SIGMAbar,T)';
end;
        
%initialize at steady state
for j = 1:sims
    Z(:,1,j,a+1) = (inv(eye(M)-A2bar)*A1bar);
end

%obtain the rest of the VAR process
for t=2:T;
    for i=1:sims;
        Z(:,t,i,a+1) = A1bar(:,:) + A2bar*Z(:,t-1,i,a+1) + E(:,t,i);
    end;
end;

%now we will estimate A, A0, and SIGMAZ from the simulations in the large Z
%matrix, store the estimates, mean estimates, and st dev of estimates in
%the matrices below

A1bar_hat = zeros(M,sims,a+1);
A1bar_mean = zeros(M,1,a+1);

A2bar_hat = zeros(M,M,sims,a+1);
A2bar_mean = zeros(M,M,a+1);

SIGMAbar_hat = zeros(M,M,sims,a+1);
SIGMAbar_mean = zeros(M,M,a+1);

%obtain estimates from Stephen's function fn_VAR
for j = 1:a+1
    for i = 1:sims
        [A1bar_hat(:,i,j), A2bar_hat(:,:,i,j), SIGMAbar_hat(:,:,i,j)] = fn_VAR(Z(:,:,i,j));
    end
end


for j = 1:a+1
    for k = 1:M
        A1bar_mean(k,1,j) = mean(A1bar_hat(k,:,j));
        for l = 1:M
            A2bar_mean(k,l,j) = mean(A2bar_hat(k,l,:,j));
            SIGMAbar_mean(k,l,j) = mean(SIGMAbar_hat(k,l,:,j));
        end
    end
end

disp('   ');
disp('Vector Process Comparison:');
disp('   ');
disp('Actual constant matrix');
disp(A1bar);
disp('Mean estimates for constant matrix');
disp('    Markov    Direct Simulation');
disp([A1bar_mean(:,1,1),A1bar_mean(:,1,2)]);
disp('Actual coefficient matrix');
disp(A2bar);
disp('Mean estimates for coefficient matrix');
disp('    Markov                        Direct Simulation');
disp([A2bar_mean(:,:,1),A2bar_mean(:,:,2)]);
disp('Actual error variance matrix');
disp(SIGMAbar);
disp('Mean estimates for error variance matrix');
disp('    Markov                        Direct Simulation');
disp([SIGMAbar_mean(:,:,1),SIGMAbar_mean(:,:,2)]);


























































%seed normal random number generator
seed_no = 16349618410;
randn('seed',seed_no);

%seed uniform random number generator
rand('seed',1111111);

%general setup/storage
M = size(A1bar,1);  %number of variables in simulated VARs
a = 1;  %2= a+1 =Markov approx. vs. direct sim
Pr_mat = zeros(prod(N),prod(N));    %used to store the probability matrices from each of the alternative functions
Pr_mat_key = zeros(M,prod(N));    %used to store the keys from each of the alternative functions
zbar = zeros(M,max(N));   %univariate grid key

%Diagonal and singular covariance
tic
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Diagonal and singular covariance')
SIGMA1 = [0.8 0 0; 0 0.6 0; 0 0 0.7];
SIGMA2 = [0.55 0 0; 0 0 0; 0 0 0.3];
disp('SIGMA1=');
disp(num2str(SIGMA1));
disp('SIGMA2=');
disp(num2str(SIGMA2));
SIGMAbar = SIGMA2*SIGMA1*(SIGMA2');
disp('SIGMAbar=');
disp(num2str(SIGMAbar));
disp('Eigenvalue of SIGMAbar (check for stability)');
disp(num2str(eig(SIGMAbar)));

[Pr_mat(:,:),Pr_mat_key(:,:),zbar(:,:)] = fn_var_to_markov(eye(M),A1bar,A2bar,SIGMAbar,N,random_draws,1);
toc

%obtain uniform random draws for use in Markov simulation
U = zeros(T,sims);
U(2:T,:) = unifrnd(0,1,T-1,sims);

%matrix for storage of simulated series, a+1 dimension because it will also
%store directly simulated series
Z = zeros(M,T,sims,a+1);

%initial state is middle state
initial_state = ((prod(N)+1)/2);    %=steady-state mean

%matrix for use in simulation of the Markov chains
Pr_sum = zeros(prod(N),prod(N));

%actually simulate the Markov chains
for j = 1:a
    Pr_sum(:,1) = Pr_mat(:,1);
    for i = 2:prod(N);
        Pr_sum(:,i) = Pr_sum(:,i-1) + Pr_mat(:,i);
    end

    old_state = initial_state;
    
    for i = 1:sims
        Z(:,1,i,j) = Pr_mat_key(:,initial_state);
    end
    
    for i = 1:sims
        for t = 2:T
            new_state=find(U(t,i)<=Pr_sum(old_state,:),1);
            Z(:,t,i,j) = Pr_mat_key(:,new_state);
            old_state = new_state;
        end
    end
end

%directly simulate the VAR process

%draw normal errors (hard cases)
E = zeros(M,T,sims);
for sim_count = 1:sims
    E(:,:,sim_count) = mvnrnd([0; 0; 0],SIGMA1,T)';
end;

for t = 1:t
    for sim_count = 1:sims
        E(:,t,sim_count) = SIGMA2*E(:,t,sim_count);
    end
end
        
%initialize at steady state
for j = 1:sims
    Z(:,1,j,a+1) = (inv(eye(M)-A2bar)*A1bar);
end

%obtain the rest of the VAR process
for t=2:T;
    for i=1:sims;
        Z(:,t,i,a+1) = A1bar(:,:) + A2bar*Z(:,t-1,i,a+1) + E(:,t,i);
    end;
end;

%now we will estimate A, A0, and SIGMAZ from the simulations in the large Z
%matrix, store the estimates, mean estimates, and st dev of estimates in
%the matrices below

A1bar_hat = zeros(M,sims,a+1);
A1bar_mean = zeros(M,1,a+1);

A2bar_hat = zeros(M,M,sims,a+1);
A2bar_mean = zeros(M,M,a+1);

SIGMAbar_hat = zeros(M,M,sims,a+1);
SIGMAbar_mean = zeros(M,M,a+1);

%obtain estimates from Stephen's function fn_VAR
for j = 1:a+1
    for i = 1:sims
        [A1bar_hat(:,i,j), A2bar_hat(:,:,i,j), SIGMAbar_hat(:,:,i,j)] = fn_VAR(Z(:,:,i,j));
    end
end


for j = 1:a+1
    for k = 1:M
        A1bar_mean(k,1,j) = mean(A1bar_hat(k,:,j));
        for l = 1:M
            A2bar_mean(k,l,j) = mean(A2bar_hat(k,l,:,j));
            SIGMAbar_mean(k,l,j) = mean(SIGMAbar_hat(k,l,:,j));
        end
    end
end

disp('   ');
disp('Vector Process Comparison:');
disp('   ');
disp('Actual constant matrix');
disp(A1bar);
disp('Mean estimates for constant matrix');
disp('    Markov    Direct Simulation');
disp([A1bar_mean(:,1,1),A1bar_mean(:,1,2)]);
disp('Actual coefficient matrix');
disp(A2bar);
disp('Mean estimates for coefficient matrix');
disp('    Markov                        Direct Simulation');
disp([A2bar_mean(:,:,1),A2bar_mean(:,:,2)]);
disp('Actual error variance matrix');
disp(SIGMAbar);
disp('Mean estimates for error variance matrix');
disp('    Markov                        Direct Simulation');
disp([SIGMAbar_mean(:,:,1),SIGMAbar_mean(:,:,2)]);