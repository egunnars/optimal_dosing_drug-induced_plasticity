figure(1); tiledlayout(1,3);
figure(2); tiledlayout(2,3);

lambda_0 = 0.04; delta_d_0 = 0.08;
k_0 = 0.0004; m_0 = 0.0004;
mu = 0.0004; nu = 0.004;
lambda_1 = 0.001;

max_dose = 100;

k = k_0;
m = 0;
plot_heatmaps_rapid_linear(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, 1);

k = 0;
m = m_0;
plot_heatmaps_rapid_linear(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, 2);

k = k_0;
m = m_0;
plot_heatmaps_rapid_linear(max_dose, lambda_0, delta_d_0, lambda_1, mu, k, nu, m, 3);