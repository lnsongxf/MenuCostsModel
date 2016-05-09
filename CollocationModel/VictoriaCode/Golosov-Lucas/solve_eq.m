function [eq,val] = solve_eq(param,glob,options)

cbarlb0    = options.cbarlb;
cbarub0    = options.cbarub;
cbarlb     = cbarlb0;
cbarub     = cbarub0;

% Storage
cbarinvec          = zeros(options.itermaxcbar,1); 
cbaroutvec         = zeros(options.itermaxcbar,1);
options.cresult = [];
tictic = tic;

for tt = (1:options.itermaxcbar)
    % 1. Update cbar 
    cbar                = (1/2)*(cbarlb+cbarub);
    % 2. Solve economy given p
    eq                  = solve_xL(cbar,param,glob,options);  
    options.cresult     = eq.c;         % Save to use as starting guess
    % 3. Record output and print
    cbarinvec(tt)      = cbar;
    cbaroutvec(tt)     = eq.cbar;
    if strcmp(options.eqprint,'Y') 
        fprintf('%2i. cbarin:\t%2.6f\tcbarout:\t%2.6f\tt:%2.1f\n',tt,cbar,eq.cbar,toc(tictic));
    end
    % 4. Set all flags
    d               = cbarinvec-cbaroutvec;
    eq.flag.exist   = ~all(sign(d)==max(sign(d))); 
    eq.flag.equi    = (abs(cbarinvec(tt)-cbaroutvec(tt))<options.tolcbar);
    eq.flag.down    = (cbarinvec(tt)>cbaroutvec(tt));
    eq.flag.up      = (cbarinvec(tt)<cbaroutvec(tt));
    % 5. Shift bounds
    cbarlb             = (eq.flag.up)*cbar    + (eq.flag.down)*cbarlb;
    cbarub             = (eq.flag.up)*cbarub  + (eq.flag.down)*cbar;
    % 6. Break if equilibrium
    if eq.flag.equi,break,end
    % 7. Plot option
    if strcmp(options.eqplot,'Y') 
       figure(888);       
       subplot(2,2,3);
       plot(cbarinvec(cbarinvec~=0),'bo','markersize',6,'markerfacecolor','b');hold on;grid on;
       plot(cbaroutvec(cbaroutvec~=0),'ro','markersize',6,'markerfacecolor','r');
       xlim([0,max(tt,12)]);
       legend('cbarin','cbarout','Location','NorthEast');
       title('Equilibrium');
       drawnow;  
    end
end



end