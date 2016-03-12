% Model parameters
param.beta      = 0.99;     % discount factor
param.delta     = 0.3;      % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phi       = 0.5;      % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigmazeta = 0.225;    % stddev of productivity shocks
param.Phi       = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom       = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
param.tauc      = 0.005;    % tolerance for forecasting rule
param.n         = 5000;     % number of firms
param.T         = 96;       % simulation length
param.S         = 25;       % simulations for computing forecasting coeffs
param.s         = 100;      % simulations for moment computations

% Law of motion
b0     = 0.015;
b1     = 0.3;
b2     = 0.25;

% Sizes of grids
Np     = 20;
Na     = 5;
Ny     = 3;
Nm     = 3;

% COMPARE WITH JAMES
% Folding a, y, Dm into a single var to get entire state space in one go
Nvar = [Na; ...
    Ny; ...
    Nm];
muvar = [0; ...
    b0 + b2*param.mu*(1-param.rhom); ...
    param.mu*(1-param.rhom)];
Avar = [param.rhoa, 0, 0; ...
    0, b1, b2*param.rhom; ...
    0, 0, param.rhom];
Svar = [param.sigmazeta; ...
    b2*param.sigmaeps; ...
    param.sigmaeps];

% Create grid and transition matrix for A, Y, dM
[grid,A,~]=tauchenvar(Nvar,muvar,Avar,Svar);
grid = exp(grid');
A = A';

% Create implied grid for inflation
agrid       = unique(grid(:,1));
ygrid       = unique(grid(:,2));
mgrid       = unique(grid(:,3));
yypmp       = gridmake(ygrid, ygrid, mgrid);
pigrid      = yypmp(:,3) - log(yypmp(:,2)) + log(yypmp(:,1));
yypmppp     = [yypmp pigrid];
Npi         = Ny*Ny*Nm;

% Form the matrix that determines whether you can transit from:
% Y --> Y', m', pi'
Pryypmppp   = sparse(Ny,Ny*Nm*Npi);
p = 1;
for m = 1:Nm*Ny
    for y = 1:Ny
        temp = sparse(Ny,Ny*Nm);
        temp(y,m) = 1;
        Pryypmppp(:,p:p+Ny*Nm-1) = temp;
        p=p+Ny*Nm;
    end
end

% First Kroenecker product: repeat vertically across m dimension
Prymypmppp     = kron(ones(Nm,1),Pryypmppp);
% Second Kroenecker product: repeat horizontally and vertically for a dim
Praymypmppp    = kron(ones(Na),Prymypmppp);
% Third Kroenecker product: repeat horizontally and vertically for p dim
Prpaympypmppp  = kron(ones(Np),Praymypmppp);

% Normalize to reflect that each transition is equally probable
Prpaympypmppp = sparse(diag(1./sum(Prpaympypmppp,2))*Prpaympypmppp);

% This is the matrix for transitioning between the N states
N_trans = kron(A,speye(Np));

% I *think* this is the matrix that should be used to form expectations
exp_matrix = sparse(N_trans*Prpaympypmppp);

