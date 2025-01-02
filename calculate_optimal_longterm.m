% finds the optimal long term solution from f_0 = 1 to equilibrium 
function [times, doses, fs] = calculate_optimal_longterm(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, points)

   lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
   mu_c = @(c) mu + k.*c;
   nu_c = @(c) nu - m.*c;

   % calculate limit properties
   [rho_l, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, max_dose);
   f_l = get_equilibrium(c_l);

   options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',10e-6,'ConstraintTolerance',10e-6, 'Display', 'none');

   f_start = get_equilibrium(0);
   f_end = f_l+10^(-10);

   % calculate all optimal doses
   fs = linspace(f_start, f_end, points);
   doses = [];
   times = [];
   tot_time = 0;
   for i = 1:numel(fs)
      f = fs(i);
      if f_start > f_end
         obj = @(c) -calculate_objective(c, f);  
      else
         obj = @(c) calculate_objective(c, f);
      end

      A = m*(f-1)-k*f;
      B = delta_d_0*f*(f-1);
      D = (lambda_1-lambda_0)*f*(f-1)-(mu+nu)*f+nu;

      if D == 0
          lb = -(A+B)/A;
      else
          lb = max(roots([A A+B+D D]))+10^(-6);
      end
      ub = max_dose;

      best_val = Inf; 
      for ell=1:10
         [dose, val] = fmincon(obj, rand * max_dose, [], [], [], [], lb, ub, [], options);
         if val < best_val
            best_dose = dose;
            best_val = val;
         end
      end
      dose = best_dose;

      doses = [doses ; dose];     
      times = [times ; tot_time];

      time = abs(1 / f0_deriv(dose, f) * (f_end - f_start) / points);
      tot_time = tot_time + time;
   end
   
   % the objective function
   function obj = calculate_objective(c, f) 
      rho_c = (lambda_0_c(c) - lambda_1) * f + lambda_1;
      f_prime = f0_deriv(c, f);
      obj = (rho_c - rho_l) ./ f_prime;
   end

   % helper function to get equilibrium ratio
   function f_0 = get_equilibrium(c) 
      A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
      
      sigma = max(eig(A));     
      f_0 = A(2,1)/(sigma-A(1,1)+A(2,1));
   end

   % a function to calcluate f0'
   function dydt = f0_deriv(c, f0)
      a = lambda_1-lambda_0_c(c);
      b = lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c);
      
      dydt = a.*f0.^2-b.*f0+nu_c(c);
   end
end
