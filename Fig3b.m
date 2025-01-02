global lambda_0_c lambda_1 mu_c nu_c rho_l;

k_0 = 0.0004; m_0 = 0.0004;

max_dose = 10;
lambda_0 = 0.04; delta_d_0 = 0.08;
mu = 0.0004; k = k_0;
lambda_1 = 0.001;
nu = 0.004; m = 0;

lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
mu_c = @(c) mu + k.*c;
nu_c = @(c) nu - m.*c;

T = 1200;
x_t = 0:1:T;

% calculate limit properties
[rho_l, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, max_dose);
f_l = get_equilibrium(c_l);

points = 1000;
[times, doses, ~] = calculate_optimal_longterm(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, points);

figure;
tiledlayout(1,3);

nexttile(1);
hold on,
ic_in = get_equilibrium(0);
[a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,ic_in);
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

plot(times,doses,'Color',[78 107 166]/255,'LineWidth',5);
plot(x_t,results.c_t,'--','LineWidth',3,'Color',[0 0 0]);
color = 0.3;
rectangle('Position',[900,1,300,99],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)
legend('Transient phase via (3)','Forward-backward sweep','Interpreter','Latex','FontSize',19)
xlim([0 T]);
ylim([1 max_dose]);
set(gca,'YScale','log');
set(gca,'fontsize', 14)
xlabel('Time $t$','Interpreter','Latex','FontSize',19);
ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','Latex','FontSize',19);

nexttile(2);

k = 0;
m = m_0;

lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
mu_c = @(c) mu + k.*c;
nu_c = @(c) nu - m.*c;

% calculate limit properties
[rho_l, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, max_dose);
f_l = get_equilibrium(c_l);

[times, doses, ~] = calculate_optimal_longterm(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, points);

hold on
ic_in = get_equilibrium(0);
[a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,ic_in);

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
          
plot(times,doses,'Color',[78 107 166]/255,'LineWidth',5);
plot(x_t,results.c_t,'--','LineWidth',3,'Color',[0 0 0]);
color = 0.3;
rectangle('Position',[900,1,300,99],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)
legend('Transient phase via (3)','Forward-backward sweep','Interpreter','Latex','FontSize',19)
xlim([0 T]);
ylim([1 max_dose]);
set(gca,'fontsize', 14)
set(gca,'YScale','log');
xlabel('Time $t$','Interpreter','Latex','FontSize',19);
ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','latex','Interpreter','Latex','FontSize',19);

nexttile(3);

k = k_0;
m = m_0;

lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
mu_c = @(c) mu + k.*c;
nu_c = @(c) nu - m.*c;

% calculate limit properties
[rho_l, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, max_dose);
f_l = get_equilibrium(c_l);

[times, doses, ~] = calculate_optimal_longterm(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, points);

hold on
ic_in = get_equilibrium(0);
[a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k,nu,m,T,max_dose,ic_in);

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

plot(times,doses,'Color',[78 107 166]/255,'LineWidth',5);
plot(x_t,results.c_t,'--','LineWidth',3,'Color',[0 0 0]);
color = 0.3;
rectangle('Position',[900,1,300,99],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)
legend('Transient phase via (3)','Forward-backward sweep','Interpreter','Latex','FontSize',19)
xlim([0 T]);
ylim([1 max_dose]);
set(gca,'fontsize', 14)
set(gca,'YScale','log');
xlabel('Time $t$','Interpreter','Latex','FontSize',19);
ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','Latex','FontSize',19);

% the objective function
function obj = calculate_objective(c, f) 
   global lambda_0_c lambda_1 rho_l;
   rho_c = (lambda_0_c(c) - lambda_1) * f + lambda_1;
   f_prime = f0_deriv(c, f);
   obj = (rho_c - rho_l) ./ f_prime;
end

% helper function to get equilibrium ratio
function f_0 = get_equilibrium(c) 
   global lambda_0_c lambda_1 mu_c nu_c;

   A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];    
   sigma = max(eig(A));
   f_0 = A(2,1)/(sigma-A(1,1)+A(2,1));
end

% a function to calcluate f0'
function dydt = f0_deriv(c, f0)
   global lambda_0_c lambda_1 mu_c nu_c;
   f = lambda_1-lambda_0_c(c);
   g = lambda_1-lambda_0_c(c)+mu_c(c)+nu_c(c);
   
   dydt = f.*f0.^2-g.*f0+nu_c(c);
end