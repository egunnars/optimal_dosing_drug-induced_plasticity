figure(1); tiledlayout(1,3);

lambda_0 = 0.04; delta_d_0 = 0.08;
mu = 0.0004; nu = 0.004;
delta_mu_0 = 0.004; delta_nu_0 = 0.003;
lambda_1 = 0.001;

max_dose = 10;

delta_mu = delta_mu_0;
delta_nu = 0;
plot_heatmaps_rapid_uniform(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu, 1);

delta_mu = 0;
delta_nu = delta_nu_0;
plot_heatmaps_rapid_uniform(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu, 2);

delta_mu = delta_mu_0;
delta_nu = delta_nu_0;
plot_heatmaps_rapid_uniform(max_dose, lambda_0, delta_d_0, lambda_1, mu, delta_mu, nu, delta_nu, 3)