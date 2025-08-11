max_dose = 10;

results_mat = zeros(90,2);
doses = zeros(90,1);

lambda_0 = 0.048;
delta_d_0 = 1.095;

mu = 0.001; k = 0.03256;
lambda_1 = 0;
nu = 0.073; m = 0;

figure(1);
tiledlayout(2,2);

T = 100;
[a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur_application(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,1);

results = oc_maeling;

results.c_t = a;
results.c_t_iter = b;
results.c_t_orig_iter = c;
results.f_0_t = d;
results.f_0_t_iter = i;
results.f_0_t_orig_iter = j;
results.gamma_t = e;
results.final_size = f;
results.final_size_vec = g;

nexttile(1);
plot([0:0.01:T],results.c_t,'LineWidth',3)
ylim([0 10^1])

[rho_l, c_l] = best_constant_dose_application(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m);
hold on
yline(c_l,'--black','LineWidth',3);

set(gca,'fontsize', 14)
xlabel('Time $t$','Interpreter','latex','FontSize',19);
ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','latex','FontSize',19);
legend('Optimal dosing schedule ($T=100$)','Optimal equilibrium dose','Interpreter','latex','FontSize',19);

T = 10;
[a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur_application(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,1);

results = oc_maeling;

results.c_t = a;
results.c_t_iter = b;
results.c_t_orig_iter = c;
results.f_0_t = d;
results.f_0_t_iter = i;
results.f_0_t_orig_iter = j;
results.gamma_t = e;
results.final_size = f;
results.final_size_vec = g;

c_max = 10;
c_av = c_max/2;

x_t = 0:0.1:T;
c_t = [0:0.01:1]*c_max;
c_t_new = c_av*ones(1,size(x_t,2));

for ell=1:960
    c_max = (5+ell-1)/1000*10;
    c_av = c_max/2;
    doses(ell) = c_av;

    c_t = [0:0.01:1]*c_max;
    c_t_new = c_av*ones(1,size(x_t,2));

    ic = 1;
    [~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
    f_0_t = y;
    final_size = integral(@(t) fin_size(t,x_t,c_t,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),T);
    results_mat(ell,1) = exp(final_size)*(1-y(end));
    ic = 1;
    [~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t_new,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
    f_0_t = y;
    final_size = integral(@(t) fin_size(t,x_t,c_t_new,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),T);
    results_mat(ell,2) = exp(final_size)*(1-y(end));
end

size_traj = zeros(1001,1);
x_t = [0:0.01:T];
nexttile(3);

c_t_new = doses(856)*ones(1,size(x_t,2)); 
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t_new,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t_new,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,10^10*exp(size_traj),'LineWidth',3);
set(gca,'YScale','log');
ylim([3*10^8 10^10]);

hold on

c_t = [0:0.001:1]*doses(856)*2;
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,10^10*exp(size_traj),'LineWidth',3);
set(gca,'YScale','log');

c_t_const = doses(166)*ones(1,size(x_t,2));
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t_const,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t_const,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,10^10*exp(size_traj),'LineWidth',3);
set(gca,'YScale','log');

dose_opt = results.c_t;
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,dose_opt,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
size_traj = zeros(1001,1);
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,dose_opt,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,10^10*exp(size_traj),'LineWidth',3);
set(gca,'YScale','log');

set(gca,'fontsize', 14)
legend('Constant $c=4.2$','Linear $\bar{c}=4.2$','Constant $c=0.85$','Optimal ($T=10$)','Interpreter','latex','FontSize',19)
xlabel('Time $t$','Interpreter','latex','FontSize',19);
ylabel('Tumor size (log scale)','Interpreter','latex','FontSize',19)

nexttile(2);
plot(x_t,c_t_new,'LineWidth',3);
hold on
plot(x_t,c_t,'LineWidth',3);
plot(x_t,c_t_const,'LineWidth',3);
plot(x_t,dose_opt,'LineWidth',3);
set(gca,'fontsize', 14)
legend('Constant $c=4.2$','Linear $\bar{c}=4.2$','Constant $c=0.85$','Optimal ($T=10$)','Interpreter','latex','FontSize',19)
xlabel('Time $t$','Interpreter','latex','FontSize',19);
ylabel('Dose ($c$) (relative to $EC_{50}^d$)','Interpreter','latex','FontSize',19)

lambda_1 = 0;
nu = 0; m = 0;

c_t_new = doses(856)*ones(1,size(x_t,2));
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t_new,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t_new,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
nexttile(4);
plot(x_t,exp(size_traj).*(1-f_0_t),'LineWidth',3);
ylim([0 0.15]);
hold on

c_t = [0:0.001:1]*doses(856)*2;
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,exp(size_traj).*(1-f_0_t),'LineWidth',3);

c_t_const = doses(166)*ones(1,size(x_t,2));
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,c_t_const,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,c_t_const,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,exp(size_traj).*(1-f_0_t),'LineWidth',3);

dose_opt = results.c_t;
ic = 1;
[~,y] = ode45(@(t,y) f0_ode(t,y,x_t,dose_opt,lambda_0,delta_d_0,mu,k,lambda_1,nu,m), x_t, ic);
f_0_t = y;
size_traj = zeros(1001,1);
for ell=1:size(size_traj,1)
    size_traj(ell) = integral(@(t) fin_size(t,x_t,dose_opt,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m),min(x_t),x_t(ell));
end
plot(x_t,exp(size_traj).*(1-f_0_t),'LineWidth',3);

set(gca,'fontsize', 14)
legend('Constant $c=4.2$','Linear $\bar{c}=4.2$','Constant $c=0.85$','Optimal ($T=10$)','Interpreter','latex','Fontsize',19)
xlabel('Time $t$','Interpreter','latex','Fontsize',19);
ylabel('Number of persister cells','Interpreter','latex','Fontsize',19)

function dydt = f0_ode(t,y,x_t,c_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m)
    lambda_0_c = @(c) lambda_0 - delta_d_0*(1-exp(-c*log(2)));
    mu_c = @(c) mu + k.*c;
    nu_c = @(c) nu - m.*c;

    f = @(c) lambda_1-lambda_0_c(c);
    g = @(c) lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c);
    
    c = interp1(x_t,c_t,t);
    
    dydt = f(c).*y.^2-g(c).*y+nu_c(c);
end

function x = fin_size(t,x_t,c_t,f_0_t,lambda_0,delta_d_0,mu,k,lambda_1,nu,m)
    c = interp1(x_t,c_t,t);
    f_0 = interp1(x_t,f_0_t,t);

    lambda_0_c = @(c) lambda_0 - delta_d_0*(1-exp(-c*log(2)));
    mu_c = @(c) mu + k.*c;
    nu_c = @(c) nu - m.*c;

    x = (lambda_0_c(c)-lambda_1).*f_0+(lambda_1);
end