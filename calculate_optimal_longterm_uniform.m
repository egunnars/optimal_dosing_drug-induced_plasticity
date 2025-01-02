% finds the optimal long term solution from f_0 = 1 to equilibrium 
function [times, doses, fs, c_l, f_l] = calculate_optimal_longterm_uniform(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu, points)

   lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
   mu_c = @(c) mu + delta_mu.*(c>0);
   nu_c = @(c) nu - delta_nu.*(c>0);

   lambda_c_phi = @(c, phi) lambda_0 - delta_d_0 * c / (c + 1) * phi;
   mu_c_phi = @(c, phi) mu + delta_mu * phi;
   nu_c_phi = @(c, phi) nu - delta_nu * phi;

   % calculate limit properties
   [rho_l, c_l] = best_constant_dose_new_phi(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu);
   f_l = get_equlibrium(max_dose,c_l);

   options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',10e-6,'ConstraintTolerance',10e-6, 'Display', 'none');

   f_start = 1;
   f_end = f_l+10^(-10);

   f_end;

   % calculate all optimal doses
   fs = linspace(f_start, f_end, points);
   doses = [];
   times = [];
   tot_time = 0;
   for i = 1:numel(fs)
      f = fs(i);
      dose = max_dose;

      doses = [doses ; dose];     
      times = [times ; tot_time];
      f0_deriv(dose, f)

      time = abs(1 / f0_deriv(dose, f) * (f_end - f_start) / points);
      tot_time = tot_time + time;
   end

   % helper function to get equilibrium ratio
   function f_0 = get_equlibrium(c,phi) 
      %A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
      A = [lambda_c_phi(c,phi)-mu_c_phi(c,phi),mu_c_phi(c,phi);nu_c_phi(c,phi),lambda_1-nu_c_phi(c,phi)];
      
      sigma = max(eig(A));
      %[V, ~] = eigs(A');
      
      f_0 = A(2,1)/(sigma-A(1,1)+A(2,1));
      %f_0 = V(1, 2) / (V(1, 2) + V(2, 2)); % TODO/WARNING: change this, this is hacky
   end

   % a function to calcluate f0'
   function dydt = f0_deriv(c, f0)
      a = lambda_1-lambda_0_c(c);
      b = lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c);
      
      dydt = a.*f0.^2-b.*f0+nu_c(c);
   end
end
