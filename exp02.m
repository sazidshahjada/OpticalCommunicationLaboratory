% Parameters for graded-index fiber
n1 = 1.48;  % Core refractive index
n2 = 1.45;  % Cladding refractive index
a = 5e-6;   % Core radius in meters
lambda_0 = 0.65e-6;  % Wavelength in meters
delta = (n1^2 - n2^2) / (2 * n1^2);  % Relative index difference

% Calculate Numerical Aperture (NA)
NA = sqrt(n1^2 - n2^2);

% Calculate V number
V = (2 * pi * a * NA) / lambda_0;

% Calculate propagation constant (β) for LP_01 mode
n_eff = n1 * sqrt(1 - 2 * delta / sqrt(V));
k = 2 * pi / lambda_0;
beta = k * n_eff;

% Calculate cutoff wavelength (λc) for single mode (V ≈ 3.4 for parabolic profile)
V_c = 3.4;
lambda_c = (2 * pi * a * NA) / V_c;

fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('V number: %.4f\n', V);
fprintf('Propagation constant (β): %.4e m^-1\n', beta);
fprintf('Cutoff wavelength (λc): %.4f μm\n', lambda_c * 1e6);

% Sweep wavelength for plots
lambdas = linspace(0.5e-6, 1.5e-6, 100);
Vs = (2 * pi * a * NA) ./ lambdas;
n_effs = n1 * sqrt(1 - 2 * delta ./ sqrt(Vs));
ks = 2 * pi ./ lambdas;  % Use wavelength-dependent k
betas = ks .* n_effs;

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