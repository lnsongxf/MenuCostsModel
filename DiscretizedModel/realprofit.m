function [ prof ] = realprofit(flag,parms,pP,a,Y)
%REALPROFIT Real profit function for Terry & Knotek (2008)
%-------------------------------------------------
%   Compute the firm's profit as a function of the relative price,
%   productivity, aggregate output/demand, money, inflation, and 
%   decision about whether to pay menu cost
%
%   INPUTS
%   - flag       = 'K' if keep old price, 'C' if change to new price
%   - parms      = structure; model parameters
%       * .grid  = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i   = scalar; size of particular state variable vector e.g. Na
%   - pP         = real price: old if flag = 'K', new if flag = 'C'
%   - a          = current productivty
%   - Y          = current output
%
%-------------------------------------------------

theta = parms.theta;
delta = parms.delta;
sigma = parms.sigma;
alpha = parms.alpha;
phi   = parms.phi;
Phi   = parms.Phi;
Npp   = parms.Npp;


switch flag
    case 'K'  % Keep old price
        prof = (Y')*(pP).^(1-theta) - ...
            delta*(Y').^sigma.*( (Y'./a').^(1/alpha)*(pP).^(-theta/alpha) ).^(1+phi);
        
    case 'C' % Change to new price, pay menu cost Phi
        prof = ((Y')*(pP)).^(1-theta) - ...
            delta*repmat(Y',1,Npp).^sigma.*( (Y'./a').^(1/alpha)*(pP).^(-theta/alpha) + Phi).^(1+phi);

end



end

