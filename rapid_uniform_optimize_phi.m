function rate = rapid_uniform_optimize_phi(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu)

global lambda_1 lambda_c_phi mu_c_phi nu_c_phi;

c_val = linspace(1, 100, 1000); % values for c 
phi_val = linspace(0.55, 1, 1000); % values for phi 

lambda_c_phi = @(c, phi) lambda_0 - delta_d_0 * c / (c + 1) * phi;
mu_c_phi = @(c, phi) mu + delta_mu * phi;
nu_c_phi = @(c, phi) nu - delta_nu * phi;

rate_c_phi = @(c, phi) 1/2 * (lambda_c_phi(c, phi) - mu_c_phi(c, phi) + lambda_1 - nu_c_phi(c, phi) + sqrt((lambda_c_phi(c, phi) - mu_c_phi(c, phi) - lambda_1 + nu_c_phi(c, phi))^2 + 4 * mu_c_phi(c, phi) * nu_c_phi(c, phi)));

opt_val = Inf;
phi_val = [0:0.001:1];
for k=1:size(phi_val,2)
    temp = rate_c_phi(max_dose,phi_val(k));
    if temp < opt_val
        opt_val = temp;
        opt_phi = phi_val(k);
    end
end

rate = opt_phi;
end