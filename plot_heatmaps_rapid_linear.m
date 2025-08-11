function plot_heatmaps_rapid_linear(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, number)

global lambda_1 lambda_c_phi mu_c_phi nu_c_phi;

c_val = linspace(1, max_dose, 500); % values for c 
phi_val = linspace(0.5, 1, 500); % values for phi 

lambda_c_phi = @(c, phi) lambda_0 - delta_d_0 * c / (c + 1) * phi;
mu_c_phi = @(c, phi) mu + k * c^2/(c^2+(5)^2) * phi;
nu_c_phi = @(c, phi) nu - m * c^2/(c^2+(5)^2) * phi;

rate_c_phi = @(c, phi) 1/2 * (lambda_c_phi(c, phi) - mu_c_phi(c, phi) + lambda_1 - nu_c_phi(c, phi) + sqrt((lambda_c_phi(c, phi) - mu_c_phi(c, phi) - lambda_1 + nu_c_phi(c, phi))^2 + 4 * mu_c_phi(c, phi) * nu_c_phi(c, phi)));

[C, PHI] = meshgrid(c_val, phi_val);

RATE = arrayfun(rate_c_phi, C, PHI);
Y_AVG = arrayfun(@f0_c_phi, C, PHI);

% calculate minimum for each c
min_phi = zeros(size(c_val));
min_phi_v = zeros(size(c_val));
min_c = zeros(size(phi_val));
min_c_v = zeros(size(phi_val));
min_phi_y_avg = zeros(size(c_val));

for i = 1:length(c_val)
   [min_v, min_i] = min(RATE(:,i));
   min_phi(1,i) = phi_val(min_i);
   min_phi_v(1,i) = min_v;
   min_phi_y_avg(1,i) = Y_AVG(min_i, i);
end

for i = 1:length(phi_val)
   [min_v, min_i] = min(RATE(i,:));
   min_c(1,i) = c_val(min_i);
   min_c_v(1,i) = min_v;
end

% plot global minimum
minMatrix = min(RATE(:));
[minrow,mincol] = find(RATE==minMatrix);

figure(1);
nexttile(number);

contourf(C, PHI, RATE, 'LineStyle', 'none');

hold on;
shading interp;

contour(C, PHI, RATE, "ShowText", true, 'LineColor', 'k', 'LevelList', -0.003:0.003:0.03);
contour(C, PHI, RATE, "ShowText", true, 'LineColor', 'k', 'LevelList', -0.003:0.001:0);
contour(C, PHI, RATE, 'LineColor', 'k', 'LineStyle', ':', 'LineWidth', 1, 'LevelList', -0.004:0.0001:-0.001);
%contour(C, PHI, RATE, "ShowText", true, 'LineColor', 'k', 'LineStyle', ':', 'LineWidth', 1, 'LevelList', -0.004:0.0004:-0.002);

%plot(c_val, min_phi, "r");
plot(min_c, phi_val,'Color',[212 10 0]/255,'LineWidth',4);
plot(c_val(mincol), phi_val(minrow),'.r', 'MarkerSize', 35);

set(gca,'fontsize', 14);
xlabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','latex','FontSize',19);
ylabel('Drug proportion $\varphi$','Interpreter','latex','FontSize',19);
set(gca,'Xscale','log');

figure(2);
nexttile(number);
fcn = zeros(1,size(phi_val,2));
for k=1:size(fcn,2)
    a = find(c_val == min_c(k));
    fcn(k) = RATE(k,a);
end
plot(phi_val, fcn,'Color',[78 107 166]/255,'LineWidth',3);
ylim([-3 0]*10^(-3));
set(gca,'fontsize', 14);
xlabel('Drug proportion $\varphi$','Interpreter','latex','FontSize',19);
legend('Long-run growth rate','Interpreter','latex','FontSize',19);
xlim([0.6 1]);

nexttile(3+number);
fcn_2 = zeros(1,size(phi_val,2));
for k=1:size(fcn_2,2)
    fcn_2(k) = min_c(k)*phi_val(k);
end
plot(phi_val, fcn_2/c_val(mincol),'Color',[78 107 166]/255,'LineWidth',3);
ylim([0.9 3]);
set(gca,'fontsize', 14);
xlabel('Drug proportion $\varphi$','Interpreter','latex','FontSize',19);
legend('Average dose (relative to optimal)','Interpreter','latex','FontSize',19);
xlim([0.6 1]);
end

function f0 = f0_c_phi(c, phi) 
   global lambda_1 lambda_c_phi mu_c_phi nu_c_phi;
   a1 = lambda_1 - lambda_c_phi(c, phi);
   a2 = -(lambda_1 - lambda_c_phi(c, phi) + mu_c_phi(c, phi) + nu_c_phi(c, phi)); 
   a3 = nu_c_phi(c, phi);

   rs = roots([a1 a2 a3]);
   if 0 <= rs(1) && rs(1) <= 1
      f0 = rs(1);
   else
      f0 = rs(2);
   end
end