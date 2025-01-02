figure;
tiledlayout(1,3);

lambda_0 = 0.04; delta_d_0 = 0.08;
mu = 0.0004; 
lambda_1 = 0.001;
nu = 0.004; 
max_dose = 10;

delta_mu = 10.^([-2.5:0.005:2]);
delta_nu = 10.^([-2.2:0.005:0]);
results = zeros(size(delta_nu,2),size(delta_mu,2));
for k=1:size(delta_mu,2)
    for j=1:size(delta_nu,2)
        results(j,k) = rapid_uniform_optimize_phi(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu(k)*mu, nu, delta_nu(j)*nu);
    end
end

nexttile(1);
contourf(delta_mu,delta_nu,results, "ShowText", true, 'LineColor', 'k', 'LevelList', 0.5:0.05:1);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'fontsize', 14);
xlabel('$\Delta \mu/\mu_0$','Interpreter','latex','FontSize',19);
ylabel('$\Delta \nu/\nu_0$','Interpreter','latex','FontSize',19)
title('Optimal exposure proportion','Interpreter','latex','FontSize',19)

nexttile(3);
plot(delta_mu,results(1,:),'--','Color',[0 0 0]/255,'LineWidth',2.5);
set(gca,'XScale','log');
ylim([0 1.1]);
set(gca,'fontsize', 14);
xlabel('$\Delta \mu/\mu_0$','Interpreter','latex','FontSize',19);
ylabel('Drug exposure proportion $(\varphi)$','Interpreter','latex','FontSize',19);
legend('Optimal proportion','Interpreter','latex','FontSize',19);

nexttile(2);
plot(delta_nu,results(:,1),'Color',[0 0 0]/255,'LineWidth',2.5);
set(gca,'XScale','log');
ylim([0 1.1]);
set(gca,'fontsize', 14);
xlabel('$\Delta \nu/\nu_0$','Interpreter','latex','FontSize',19);
ylabel('Drug exposure proportion $(\varphi)$','Interpreter','latex','FontSize',19);
legend('Optimal proportion','Interpreter','latex','FontSize',19);