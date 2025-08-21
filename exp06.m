% Sellmeier coefficients for pure silica
A1 = 0.6961663; A2 = 0.4079426; A3 = 0.8974794;
B1 = 0.0684043; B2 = 0.1162414; B3 = 9.896161;

% Wavelength range (μm) based on the figure
lambda_um = linspace(0.6, 2.0, 1000);

% Refractive index using Sellmeier equation
function n = refractive_index(lambda_um, A1, A2, A3, B1, B2, B3)
    lambda2 = lambda_um.^2;
    n2 = 1 + (A1 * lambda2 ./ (lambda2 - B1^2)) + ...
             (A2 * lambda2 ./ (lambda2 - B2^2)) + ...
             (A3 * lambda2 ./ (lambda2 - B3^2));
    n = sqrt(n2);
end

% Material dispersion function
function [D_m, lambda_um_disp] = material_dispersion(lambda_um, n)
    c = 3e8;  % Speed of light (m/s)
    lambda_m = lambda_um * 1e-6;  % Convert to meters
    d_lambda = lambda_m(2) - lambda_m(1);
    dn_dlambda = gradient(n, d_lambda);
    d2n_dlambda2 = gradient(dn_dlambda, d_lambda);
    D_m = - (lambda_m / c) .* d2n_dlambda2 * 1e12 / 1e-9;  % ps/(nm·km)
    lambda_um_disp = lambda_um(6:end-5);
    D_m = D_m(6:end-5);
end

% Compute refractive index
n = refractive_index(lambda_um, A1, A2, A3, B1, B2, B3);

% Compute material dispersion
[D_m, lambda_um_disp] = material_dispersion(lambda_um, n);

% Plot
figure('Position', [100, 100, 1200, 500]);
subplot(1, 2, 1);
plot(lambda_um, n, 'b-', 'LineWidth', 2);
xlabel('Wavelength (μm)');
ylabel('Refractive Index (n)');
title('Refractive Index vs. Wavelength for Pure Silica');
grid on;

subplot(1, 2, 2);
plot(lambda_um_disp, D_m, 'r-', 'LineWidth', 2);
xlabel('Wavelength (μm)');
ylabel('Material Dispersion (ps/(nm·km))');
title('Material Dispersion vs. Wavelength for Pure Silica');
yline(0, 'k--', 'LineWidth', 0.5);
grid on;
