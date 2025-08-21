% Define the parameters
n1 = 1.48;  % Refractive index of core
n2 = 1.45;  % Refractive index of cladding
lambda_0 = 0.65e-6;  % Wavelength in meters
a = 5e-6;  % Core radius in meters
k = 2 * pi / lambda_0;  % Wave number

% Calculate Numerical Aperture (NA)
NA = sqrt(n1^2 - n2^2);

% Calculate V number
V = (2 * pi * a * NA) / lambda_0;

% Calculate propagation constant (β)
b = 1 - (1.1428 / V)^2;
n_eff = sqrt(n2^2 + b * (n1^2 - n2^2));
beta = k * n_eff;

% Calculate cutoff wavelength (λc) for single mode (V=2.405)
lambda_c = (2 * pi * a * NA) / 2.405;

fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('V number: %.4f\n', V);
fprintf('Propagation constant (β): %.4e m^-1\n', beta);
fprintf('Cutoff wavelength (λc): %.4f μm\n', lambda_c * 1e6);

% Sweep wavelength to get a range of V numbers and propagation constants
lambdas = linspace(0.5e-6, 1.5e-6, 100);
Vs = (2 * pi * a * NA) ./ lambdas;
bs = 1 - (1.1428 ./ Vs).^2;
n_effs = sqrt(n2^2 + bs .* (n1^2 - n2^2));
betas = k * n_effs;  % Note: k is fixed at lambda_0

% Create figure with subplots
figure('Position', [100, 100, 1400, 500]);
subplot(1, 2, 1);
plot(Vs, betas, 'b-', 'LineWidth', 2);
xlabel('V Number');
ylabel('Propagation Constant β (m^{-1})');
title('Propagation Constant vs V Number');
grid on;

subplot(1, 2, 2);
plot(lambdas * 1e6, Vs, 'r-', 'LineWidth', 2);
xlabel('Wavelength (μm)');
ylabel('V Number');
title('Wavelength vs V Number');
grid on;