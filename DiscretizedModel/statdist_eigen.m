function [ stat_density ] = statdist_eigen( parms, pPfun )
%STATDIST_EIGEN Compute the stationary distribution of prices using the
%eigenvector method
%------------------------------------------------------
%   INPUTS
%   - parms      = structure; model parameters
%       * .grid  = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i   = scalar; size of particular state variable vector e.g. Na
%   - pPfun      = function; price of a firm given states and current price
%
%   OUTPUT
%   - stat_density = Stationary distribution of prices
%-------------------------------------------------


Qpp = zeros(parms.Npp*parms.Na,parms.Npp);
row = 1;
for a = 1:parms.Na
    for pP = 1:parms.Npp
        a_val = parms.grid(1,a);
        pP_val = parms.pPgrid(1,pP);
        
        lower = find(pPfun(pP_val,a_val) >= parms.pPgrid);
        lower = lower(end);
        upper = find(pPfun(pP_val,a_val) < parms.pPgrid);        
        
        if isempty(upper)
            upper = length(parms.pPgrid);
            Qpp(row,upper) = 1;
        else
            upper = upper(1);
            Qpp(row,lower) = (parms.pPgrid(upper) - pPfun(pP_val,a_val))/...
                (parms.pPgrid(upper) - parms.pPgrid(lower));
            Qpp(row,upper) = (pPfun(pP_val,a_val) - parms.pPgrid(lower))/...
                (parms.pPgrid(upper) - parms.pPgrid(lower));
        end        
        row = row + 1;
        
    end
end

Qa = kron(parms.trans_a,ones(parms.Npp,1));
col = 1;
Q = zeros(parms.Npp*parms.Na,parms.Npp*parms.Na);
for a = 1:parms.Na
    for pP = 1:parms.Npp
        Q(:,col) = Qpp(:,pP).*Qa(:,a);        
        col = col + 1;
    end
end
    
% Pertub the Q matrix so eigenvalues are unique
eta = min(nonzeros(Q))/(2*parms.Npp);
index = find(Q == 0);
Q(index) = eta;

for i = 1:size(Q,1)
    Q(i,:) = Q(i,:)/(sum(Q(i,:)));
end

% Find eigenvector ==> stationary distribution
[V,D,W] = eig(Q);
V = real(V);
D = real(D);
W = real(W);
index = find(diag(D) > 0.999999);

% Compute the stationary distribution of prices
stat_density = W(:,index)'/sum(W(:,index));


%% Sample from stationary distribution and plot
% Draw from stationary distribution to find distribution of prices
I = 100000;         % Sample size
stat_dist_pPgrid = repmat(parms.pPgrid,1,parms.Na);      % 

sample_dist = randsample(stat_dist_pPgrid,I,true,stat_density);

figure
subplot(1,2,1)
hist(sample_dist,parms.Npp);
title('Distribution of prices (histogram)','fontsize',12)
set(gca,'fontsize',12)
grid on
subplot(1,2,2)
bandwidth = 0.1;
[f,xi] = ksdensity(sample_dist);  %,'width',bandwidth);
plot(xi,f,'linewidth',2);
title('Distribution of prices (kernel density)','fontsize',12)
set(gca,'fontsize',12)
grid on


end

