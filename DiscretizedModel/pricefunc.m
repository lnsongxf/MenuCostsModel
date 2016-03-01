function [ pPout ] = pricefunc(parms,pP,astate,dmstate,Ystate,pistate,Vk,Vc,V)
%VALFUN Summary of this function goes here
%-------------------------------------------------
%   Compute the firm's value function as function of past value function,
%   parameters, the old relative price, productivity, 
%   aggregate output/demand, money, inflation, and decision about whether 
%   to pay menu cost.
%
%   Do this one state-space element at a time, because different points in
%   the space lead to either 'K' or 'C'
%
%   INPUTS
%   - flag          = 'K' if keep old price, 'C' if change to new price
%   - parms         = structure; model parameters
%       * .Ppgrid   = vector; grid points for real prices (Npp x 1) 
%       * .grid     = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i      = scalar; size of particular state variable vector e.g. Na
%       * .trans    = matrix; transition prob matrix (Na*Ndm*Ny*Npi x Na*Ndm*Ny*Npi)
%   - pP            = current real price
%   - statenum      = index number representing the current state 
%   - V             = matrix; value function (Na*Ndm*Ny*Npi x Np)
%
%   OUTPUTS
%   - pPout = new price. If 'K', same as before. If 'C', new argmax.
%-------------------------------------------------

% Find statenumber in grid
statenum = find(ismember(parms.grid(1,:),astate) & ismember(parms.grid(2,:),dmstate) ...
    & ismember(parms.grid(3,:),Ystate) & ismember(parms.grid(4,:),pistate));

trans = parms.trans(statenum,:);
pPgrid = parms.pPgrid;

% Check whether firm keeps or changes price
tmp = abs( pPgrid - pP );   % how far from grid?
[~, idx] = min(tmp);     % index of closest values   
pP = pPgrid(idx);       % Put price on the grid

if Vk(statenum,idx) > Vc(statenum,idx)  % Keep price
    pPout = pP;

elseif Vk(statenum,idx) <= Vc(statenum,idx)  % Change price
    val2max =  realprofit('C',parms,pP,astate,Ystate)  + parms.beta*trans*V; 
    [~, idx] = max(val2max); % max value and pPgrid index of max value
    pPout = pPgrid(idx);
    

end




end


