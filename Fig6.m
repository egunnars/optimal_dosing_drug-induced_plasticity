lambda_0 = 0.04;
delta_d_0 = 0.08;
mu = 0.0004;
delta_mu = [0.004,0,0.004];
lambda_1 = 0.001;
nu = 0.004;
delta_nu = [0,0.003,0.003];

max_dose = 10;
T = 1000;

points = 1000;
figure;
tiledlayout(2,3);

for k=1:3
    [times, doses, fs, c_l, f_l] = calculate_optimal_longterm_uniform(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu(k), nu, delta_nu(k), points);
    
    nexttile(k);
    plot(times, doses/max_dose,'LineWidth',3,'Color',[78 107 166]/255)
    hold on;
    times_new = times(end)+(T-times(end))*[0:0.01:1];
    plot(times_new,c_l*ones(1,size(times_new,2)),'LineWidth',3,'Color',[78 107 166]/255);
    yline(c_l,'--black','LineWidth',2.5);
    ylim([0 1.1]);
    set(gca,'fontsize', 14)
    xlabel('Time $t$','Interpreter','Latex','FontSize',19);
    ylabel('Drug proportion $\varphi$','Interpreter','Latex','FontSize',19);
    legend('Optimal treatment','','Optimal equilibrium','Interpreter','Latex','FontSize',19);
    xlim([0 T]);
    
    nexttile(k+3);
    plot(times, fs,'LineWidth',3,'Color',[212 10 0]/255);
    hold on
    plot(times_new,f_l*ones(1,size(times_new,2)),'LineWidth',3,'Color',[212 10 0]/255);
    yline(f_l,'--black','LineWidth',3);
    ylim([0 1.1]);
    set(gca,'fontsize', 14)
    xlabel('Time $t$','Interpreter','Latex','FontSize',19);
    ylabel('Sensitive cell proportion $f_0$','Interpreter','Latex','FontSize',19);
    legend('Optimal treatment','','Optimal equilibrium','Interpreter','Latex','FontSize',19);
    xlim([0 T]);
end