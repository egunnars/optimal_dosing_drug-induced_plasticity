% finds the optimal constant dose and the corresponding constant rate
function [best_rate, best_dose] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, c_ub)

lb = 0;
ub = min(c_ub,max_dose);
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',10e-10,'ConstraintTolerance',10e-10, 'Display', 'none');
best_rate = Inf;
for ell=1:100
    [dose, rate] = fmincon(@(c) growthrate_new(c,lambda_0,delta_d_0,lambda_1,mu,k,nu,m),rand*max_dose,[],[],[],[],lb,ub,[],options);
    if rate < best_rate
        best_rate = rate;
        best_dose = dose;
    end
end

function x = growthrate_new(c,lambda_0,delta_d_0,lambda_1,mu,k,nu,m)
    lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
    mu_c = @(c) mu + k.*c;
    nu_c = @(c) nu - m.*c;

    A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
    x = max(eig(A));
end

end
