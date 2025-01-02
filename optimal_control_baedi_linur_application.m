% updated version of optimal controll where mu and nu are functions with slope 
function [c_t,c_t_iter,c_t_orig_iter,f_0_t,f_0_t_iter,f_0_t_orig_iter,gamma_t,final_size,final_size_vec] = optimal_control_baedi_linur_application(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,ic_in)

iter = 1000;
omega = 0.95;

x_t = 0:0.01:T;
c_t = zeros(1, size(x_t, 2));
c_t_iter = zeros(iter,size(x_t,2));
c_t_orig_iter = zeros(iter,size(x_t,2));
f_0_t_iter = zeros(iter,size(x_t,2));

final_size_vec = size(iter,1);

lambda_0_c = @(c) lambda_0 - delta_d_0.*(1-exp(-c.*log(2)));
mu_c = @(c) mu + k.*c;
nu_c = @(c) nu - m.*c;

for ell=1:iter
    options = odeset('RelTol',1e-6,'AbsTol',1e-7);
    [~,y] = ode23s(@(t,y) f0_ode(t,y,x_t,c_t), x_t, ic_in, options);
    f_0_t = y;
    f_0_t_iter(ell,:) = y;

    final_size_vec(ell) = integral(@(t) fin_size(t,x_t,c_t,f_0_t),1,T);
    
    ic = 0;
    [~,y] = ode23s(@(t,y) gamma_ode(t,y,x_t,c_t,f_0_t), flip(x_t), ic, options);
    gamma_t = flip(y);

    for emm=1:size(x_t,2)
        c_t_old = c_t(emm);
        
        a = f_0_t(emm)*(gamma_t(emm)*f_0_t(emm)-1-gamma_t(emm))*delta_d_0*(log(2));
        b = -gamma_t(emm)*k*f_0_t(emm) + gamma_t(emm)*m*(f_0_t(emm)-1);

        r = log(-b/a)/(-log(2));
        count = 0;

        for err = 1:size(r,1)
            root = r(err);
            if abs(imag(root)) < 10^(-10) && real(root)>0 && real(root)<max_dose
                    c_t(emm) = omega*c_t_old+(1-omega)*root;
                    count = count+1;
                    c_t_orig_iter(ell,emm) = root;
            end
        end

        if count == 0
            if H(f_0_t(emm),0,gamma_t(emm)) < H(f_0_t(emm),max_dose,gamma_t(emm))
                c_t(emm) = omega*c_t_old+(1-omega)*0;
                c_t_orig_iter(ell,emm) = 0;
            else
                c_t(emm) = omega*c_t_old+(1-omega)*max_dose;
                c_t_orig_iter(ell,emm) = max_dose;
            end
        end
    end
end

f_0_t_orig_iter = f_0_t_iter;

[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t), x_t, ic_in, options);
f_0_t = y;
final_size = integral(@(t) fin_size(t,x_t,c_t,f_0_t),min(x_t),T);

function H = H(f_0,c,gamma)
    u = (lambda_0_c(c)-lambda_1)*f_0+lambda_1;
    g = (lambda_1-lambda_0_c(c))*f_0^2-(lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c))*f_0+nu_c(c);
    H = u+gamma*g;
end

function dydt = f0_ode(t,y,x_t,c_t)
    f = @(c) lambda_1-lambda_0_c(c);
    g = @(c) lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c);
    
    c = interp1(x_t,c_t,t);
    
    dydt = f(c).*y.^2-g(c).*y+nu_c(c);
end

function dydt = gamma_ode(t,y,x_t,c_t,f_0_t)
    f = @(c,f_0) 2*(lambda_1-lambda_0_c(c))*f_0-(lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c));
    g = @(c,f_0) lambda_0_c(c)-lambda_1;
    
    c = interp1(x_t,c_t,t);
    f_0 = interp1(x_t,f_0_t,t);
    
    dydt = -y.*f(c,f_0)-g(c,f_0);
end

function x = fin_size(t,x_t,c_t,f_0_t)
    c = interp1(x_t,c_t,t);
    f_0 = interp1(x_t,f_0_t,t);

    x = (lambda_0_c(c)-lambda_1).*f_0+(lambda_1);
end

end
