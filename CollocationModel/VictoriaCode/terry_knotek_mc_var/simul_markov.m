function [path path_ind] = simul_markov(grid,trans,T)

%{
This function simulates a Markov chain

Created by:     Victoria Gregory
Date:           9/14/2015

Inputs:         grid = grid of states
                trans = transition matrix
                T = length of simulation
                
Outputs:        path = simulated path along the states

%}

path_index = zeros(T,1);
path_val = zeros(T,1);

% draw an intial state (all states equal probability)
initial = randsample(length(grid),1);

path_index(1) = initial ;        % initialize position
path_val(1) = grid(path_index(1));

% simulate
for t=2:T
    index = path_index(t-1);
    pi_row = trans(index,:);
    path_val(t) = randsample(grid,1,true,pi_row);
    path_vec = grid == path_val(t);
    path_index(t) = find(path_vec==1);
    
end

path = path_val;
path_ind = path_index;

end