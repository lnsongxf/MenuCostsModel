function out = menufun_menucosts_KS(flag,s,x,param,glob,options)
% x=candidate new price

% Parameters
theta   = param.theta;
delta   = param.delta;
sigma   = param.sigma;
alpha   = param.alpha;
Phi     = param.Phi;
phi     = param.phi;

switch flag
    case 'change'
        a   = s(:,2);
        Y   = s(:,4);
        out = Y.*x.^(1-theta) - delta.*Y.^(sigma).*...
            ((Y./a).^(1./alpha).*x.^(-theta./alpha) + Phi).^(1+phi);
    case 'keep'
        p   = s(:,1);
        a   = s(:,2);
        Y   = s(:,4);
        out = Y.*p.^(1-theta) - delta.*Y.^(sigma).*...
            ((Y./a).^(1./alpha).*p.^(-theta./alpha)).^(1+phi);
end

end