function x = My_GPSR_Basic(y, A, AT, tau, tolA)

% Set the defaults for the optional parameters




maxiter = 2500;
max_linear_search = 5;
miniter = 5;

verbose = 1;      %report
% sufficient decrease parameter for GP line search
mu = 0.1;
% backtracking parameter for line search;
alpha_backtrack = 0.5;


%init x 是对y做一次模糊后的正变换？
x = AT(y);


% Precompute A'*y since it'll be used a lot
Aty = AT(y);

% initialize u and v
u =  x.*(x >= 0);
v = -x.*(x <  0);

% define the indicator vector or matrix of nonzeros in x
nz_x = (x ~= 0.0);
num_nz_x = sum(nz_x(:));

% Compute and store initial value of the objective function
resid =  y - A(x);
f = 0.5*(resid(:)'*resid(:)) + sum(tau(:).*u(:)) + sum(tau(:).*v(:));


% store given tau, because we're going to change it in the
% continuation procedure
final_tau = tau;

% store given stopping criterion and threshold, because we're going 
% to change them in the continuation procedure
final_stopCriterion = 1;
final_tolA = tolA;

cont_factors = 1;
cont_steps = 1;

iter = 1;

% loop for continuation
for cont_loop = 1:cont_steps

    tau = final_tau * cont_factors(cont_loop);
    
    if verbose
        fprintf(1,'\nSetting tau = %8.4f\n',tau)
    end
    
    if cont_loop == cont_steps
       stopCriterion = final_stopCriterion;
       tolA = final_tolA;
    else 
       stopCriterion = 3;
       tolA = 1e-3;
    end
    
    % Compute and store initial value of the objective function
    resid =  y - A(x);   %x = AT(y)
    f = 0.5*(resid(:)'*resid(:)) + sum(tau(:).*u(:)) + sum(tau(:).*v(:));

    objective(iter) = f;
    
    % Compute the useful quantity resid_base
    resid_base = y - resid;    %resid_base = A(x) , resid = y-A(x),  x = AT(y)
    
    keep_going = 1;

    if verbose
       fprintf(1,'\nInitial obj=%10.6e, nonzeros=%7d\n', f,num_nz_x);
    end
    
    while keep_going

      x_previous = x;

      % compute gradient
      temp = AT(resid_base);
      term  =  temp - Aty;
      %gradu = 目标函数对u求导， term = -AT(y-A(u-v)) = -AT(y-Ax))
      gradu =  term + tau; 
      gradv = -term + tau;

      % set search direction
      %du = -gradu; dv = -gradv; dx = du-dv; 
      dx = gradv-gradu;
      old_u = u; old_v = v;


      % use instead a first guess based on the "conditional" direction
      condgradu = ((old_u>0) | (gradu<0)) .* gradu;
      condgradv = ((old_v>0) | (gradv<0)) .* gradv;
      auv_cond = A(condgradu-condgradv);
      dGd_cond = auv_cond(:)'*auv_cond(:);
      alpha0 = (gradu(:)'*condgradu(:) + gradv(:)'*condgradv(:)) / (dGd_cond + realmin);
                
      % loop to determine steplength, starting wit the initial guess above.
      alpha =alpha0; 
      
      for i = 1:max_linear_search
        % calculate step for this lambda and candidate point
        du = max(u-alpha*gradu,0.0) - u; 
        u_new = u + du;
        dv = max(v-alpha*gradv,0.0) - v; 
        v_new = v + dv;
        dx = du-dv; 
        x_new = x + dx;
        
        % evaluate function at the candidate point
        resid_base = A(x_new);
        resid = y - resid_base;
        f_new = 0.5*(resid(:)'*resid(:)) + sum(tau(:).*u_new(:)) + sum(tau(:).*v_new(:));    
        if f_new <= f + mu * (gradu'*du + gradv'*dv)
          %disp('OK')  
          break
        end
        alpha = alpha * alpha_backtrack;
        fprintf(1,'    reducing lambda to %6.2e\n', alpha)
      end
      
      
      u = u_new; 
      v = v_new; 
      prev_f = f; 
      f = f_new;
      uvmin = min(u,v); 
      u = u - uvmin; 
      v = v - uvmin; 
      x = u-v;
      
      % calculate nonzero pattern and number of nonzeros (do this *always*)
      nz_x_prev = nz_x;
      nz_x = (x~=0.0);
      num_nz_x = sum(nz_x(:));
      
      iter = iter + 1;
      objective(iter) = f;
      alphas(iter) = alpha;
      
      % continue if not yeat reached target value tolA
      criterionObjective = abs((f-prev_f)/(f));
      keep_going = (criterionObjective > tolA);
      
      %keep_going = (f > tolA);
      if verbose
         fprintf(1,'It = %d Objective = %e (target = %e)\n',iter, criterionObjective, tolA) 
      end
      
      % take no less than miniter... 
    if iter<=miniter
       keep_going = 1;
    else %and no more than maxiter iterations  
        if iter > maxiter
            keep_going = 0;
        end
    end
      
    end % end of the main loop of the GP algorithm

end % end of the continuation loop
