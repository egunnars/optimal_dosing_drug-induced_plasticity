% finds the optimal constant dose and the corresponding constant rate
function [best_rate, best_dose] = best_constant_dose_new_phi(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu)

lambda_c_phi = @(c, phi) lambda_0 - delta_d_0 * c / (c + 1) * phi;
mu_c_phi = @(c, phi) mu + delta_mu * phi;
nu_c_phi = @(c, phi) nu - delta_nu * phi;

rate_c_phi = @(c, phi) 1/2 * (lambda_c_phi(c, phi) - mu_c_phi(c, phi) + lambda_1 - nu_c_phi(c, phi) + sqrt((lambda_c_phi(c, phi) - mu_c_phi(c, phi) - lambda_1 + nu_c_phi(c, phi))^2 + 4 * mu_c_phi(c, phi) * nu_c_phi(c, phi)));

lb = 0;
ub = 1;
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',10e-10,'ConstraintTolerance',10e-10, 'Display', 'none');
best_rate = Inf;
for ell=1:1000
    [dose, rate] = fmincon(@(phi) rate_c_phi(max_dose,phi),rand,[],[],[],[],lb,ub,[],options);
    if rate < best_rate
        best_rate = rate;
        best_dose = dose;
    end
end

end
