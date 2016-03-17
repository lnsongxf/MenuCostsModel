function [c,v] = solve_cKS(cold,cKS,param,glob,options)

    s           = glob.s;
    totaltic    = tic;
    
    %% Bellman iteration
    for citer = (1:options.Nbell)
        glob.citer  = citer;
        % 1. Compute values;
        v           = solve_valfuncKS(cold,s,param,glob,options); 
        % 2. Update c
        ck          = glob.Phi\full(v.vk);   
        cc          = glob.Phi\full(v.vc); 
        ce          = glob.Phi\full(v.ve); 
        c           = [ck;cc;ce];
        % 3. Compute distance and update
        dc          = norm(c-cold)/norm(cold); 
        cold        = c;
        if strcmp(options.print,'Y');
            fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
        end
    end
    
    %% Newton iterations
    if strcmp(options.print,'Y');
        fprintf('~~~~~ Newton iterations ~~~~~\n');
    end
    eq.flag.cconv = false;
    for citer = (1:options.Nnewt)
        % 1. Compute values
        [v,jac]     = solve_valfuncKS(cold,s,param,glob,options);  
        % 2. Update c 
        ckold = cold(1:end/3);  
        ccold = cold(end/3+1:2*end/3);
        ceold = cold(2*end/3+1:end);
        c               = cold - jac\([glob.Phi*ckold - full(v.vk) ; 
                                   glob.Phi*ccold - full(v.vc) ;
                                   glob.Phi*ceold - full(v.ve)]);   
        % 3. Compute distances and update
        dc          = norm(c-cold)/norm(cold);
        cold        = c;
        if strcmp(options.print,'Y');
            fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
        end
        % 4. Check convergence
        if (dc<options.tolc)
            eq.flag.cconv = true;
        end
        if eq.flag.cconv,break,end;
    end


end