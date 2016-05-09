function out = menufun(flag,s,x,cbar,param,glob,options)

    % x=candidate new price

    % Parameters
    epsilon = param.epsilon;
    gamma   = param.gamma;
    alpha   = param.alpha;

    switch flag
        case 'keep'
            x   = s(:,1);
            nu  = s(:,2);
            out = cbar.^(1-epsilon*gamma).*(alpha.*x).^(-epsilon).*(x-1./nu);
        case 'change'
            nu  = s(:,2);
            out = cbar.^(1-epsilon*gamma).*(alpha.*x).^(-epsilon).*(x-1./nu);
    end

end