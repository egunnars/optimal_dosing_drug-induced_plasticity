% updated version of optimal controll where mu and nu are functions with slope 
function [c_t,c_t_iter,c_t_orig_iter,f_0_t,f_0_t_iter,f_0_t_orig_iter,gamma_t,final_size,final_size_vec] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,ic_in)

iter = 500;
omega = 0.95;

x_t = 0:1:T;
c_t = zeros(1, size(x_t, 2));
c_t_iter = zeros(iter,size(x_t,2));
c_t_orig_iter = zeros(iter,size(x_t,2));
f_0_t_iter = zeros(iter,size(x_t,2));

final_size_vec = size(iter,1);
error_vec = zeros(iter,1);

lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
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

    for err=1:size(x_t,2)
        c_t_old = c_t(err);
        
        a = f_0_t(err)*(delta_d_0)*(gamma_t(err)*f_0_t(err)-1-gamma_t(err));
        b = -gamma_t(err)*k*f_0_t(err) + gamma_t(err)*m*(f_0_t(err)-1);

        p_2 = b;
        p_1 = 2*b;
        p_0 = a+b;
        p = [p_2, p_1, p_0];

        r = roots(p);
        count = 0;

        for em = 1:size(r,1)
            root = r(em);
            if abs(imag(root)) < 10^(-10) && real(root)>0 && real(root)<max_dose
                second_der = 2*root*p_2+p_1;
                if real(second_der)>0
                    c_t(err) = omega*c_t_old+(1-omega)*root;
                    count = count+1;
                    c_t_orig_iter(ell,err) = root;
                end
            end
        end

        if count == 0
            if H(f_0_t(err),0,gamma_t(err)) < H(f_0_t(err),max_dose,gamma_t(err))
                c_t(err) = omega*c_t_old+(1-omega)*0;
                c_t_orig_iter(ell,err) = 0;
            else
                c_t(err) = omega*c_t_old+(1-omega)*max_dose;
                c_t_orig_iter(ell,err) = max_dose;
            end
        end
    end
   
    c_t_iter(ell,:) = c_t;
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
