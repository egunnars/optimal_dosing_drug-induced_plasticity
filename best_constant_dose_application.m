% finds the optimal constant dose and the corresponding constant rate
function [best_rate, best_dose] = best_constant_dose_application(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m)

lb = 0;
ub = max_dose;
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',10e-10,'ConstraintTolerance',10e-10, 'Display', 'none');
best_rate = Inf;
for ell=1:1000
    [dose, rate] = fmincon(@(c) growthrate(c,lambda_0,delta_d_0,lambda_1,mu,k,nu,m),rand*max_dose,[],[],[],[],lb,ub,[],options);
    if rate < best_rate
        best_rate = rate;
        best_dose = dose;
    end
end

function x = growthrate(c,lambda_0,delta_d_0,lambda_1,mu,k,nu,m)
    lambda_0_c = @(c) lambda_0 - delta_d_0.*(1-exp(-c.*log(2)));
    mu_c = @(c) mu + k.*c;
    nu_c = @(c) nu - m.*c;

    A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
    x = max(eig(A));
end

end
