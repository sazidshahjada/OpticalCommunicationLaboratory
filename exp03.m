% Parameters for step-index fiber
n1 = 1.46;  % Core refractive index
n2 = 1.45;  % Cladding refractive index
a = 5e-6;   % Core radius (m)
c = 3e8;    % Speed of light (m/s)

% Function to calculate Numerical Aperture (NA)
function na = calculate_na(n1, n2)
    na = sqrt(n1^2 - n2^2);
end

% Function to calculate Acceptance angle (θa) in degrees
function theta_a_deg = calculate_acceptance_angle(na)
    theta_a_deg = rad2deg(asin(na));
end

% Function to calculate V number
function v = calculate_v_number(a, na, wavelength)
    v = (2 * pi * a * na) ./ wavelength;
end

% Function to calculate normalized propagation constant (b)
function b = calculate_b(v)
    b = 1 - (1.1428 ./ v).^2;
end

% Function to calculate effective refractive index
function n_eff = calculate_n_eff(n2, b, n1)
    n_eff = sqrt(n2^2 + b .* (n1^2 - n2^2));
end

% Function to calculate propagation constant (β)
function beta = calculate_beta(k0, n_eff)
    beta = k0 .* n_eff;
end

% Function to calculate waveguide dispersion (D_w) in ps/(nm·km)
function [D_w, lambda_range_disp] = calculate_waveguide_dispersion(lambda_range, beta_range, c)
    d_lambda = lambda_range(2) - lambda_range(1);
    d_beta = gradient(beta_range, d_lambda);  % First derivative
    d2_beta = gradient(d_beta, d_lambda);    % Second derivative
    lambda_range_disp = lambda_range(6:end-5);  % Adjust for length reduction
    d2_beta = d2_beta(6:end-5);             % Trim to match
    D_w = - (lambda_range_disp ./ c) .* d2_beta * 1e6;  % Convert to ps/(nm·km)
end

% Calculate NA
NA = calculate_na(n1, n2);

% Calculate Acceptance angle
theta_a_deg = calculate_acceptance_angle(NA);

% Wavelength range (500 nm to 2000 nm)
lambda_range = linspace(500e-9, 2000e-9, 100);

% Calculate V number, b, n_eff, and beta for the range
V_range = calculate_v_number(a, NA, lambda_range);
b_range = calculate_b(V_range);
n_eff_range = calculate_n_eff(n2, b_range, n1);
k0_range = 2 * pi ./ lambda_range;
beta_range = calculate_beta(k0_range, n_eff_range);

% Calculate waveguide dispersion
[D_w_range, lambda_range_disp] = calculate_waveguide_dispersion(lambda_range, beta_range, c);

% Calculate V number at a reference wavelength (e.g., 650 nm for consistency)
V_0 = calculate_v_number(a, NA, 650e-9);

% Print results
fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('Acceptance angle (θa): %.4f degrees\n', theta_a_deg);
fprintf('V number at 650 nm: %.4f\n', V_0);

% Plot
figure('Position', [100, 100, 800, 500]);
plot(lambda_range_disp*1e6, D_w_range, 'b-', 'LineWidth', 2);
xlabel('Wavelength (μm)');
ylabel('Waveguide Dispersion D_w (ps/nm/km)');
title('Waveguide Dispersion vs Wavelength (Step-Index Fiber)');
grid on;
set(gca, 'FontSize', 12);
axis tight;