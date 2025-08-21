% Parameters for graded-index fiber
n1 = 1.5;  % Core refractive index at center
n2 = 1.48;  % Cladding refractive index
delta = (n1^2 - n2^2) / (2 * n1^2);  % Relative refractive index difference
alpha_values = [1, 2, 4];  % Profile parameters
r_a = linspace(0, 1, 100);  % Normalized radius (r/a) from 0 to 1

% Function to calculate refractive index n(r) for a given alpha
function n_r = refractive_index(r_a, alpha, n1, delta)
    n_r = n1 * sqrt(1 - 2 * delta * (r_a .^ alpha));
end

% Function to calculate acceptance angle theta_a (in degrees)
function theta_a = acceptance_angle(r_a, alpha, n1, n2, delta)
    n_r = refractive_index(r_a, alpha, n1, delta);
    NA = sqrt(n_r.^2 - n2^2);  % Numerical aperture
    theta_a = rad2deg(asin(NA));  % Convert to degrees
    theta_a(n_r <= n2) = 0;  % Set to 0 if n(r) <= n2
end

% Create figure with subplots
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Acceptance Angle vs. Normalized Core Radius
subplot(1, 2, 1);
hold on;
for i = 1:length(alpha_values)
    alpha = alpha_values(i);
    theta_a = acceptance_angle(r_a, alpha, n1, n2, delta);
    plot(r_a, theta_a, 'LineWidth', 2, 'DisplayName', ['α = ' num2str(alpha)]);
end
xlabel('Normalized Core Radius (r/a)');
ylabel('Acceptance Angle (θa, degrees)');
title('Acceptance Angle vs. Core Radius');
legend;
grid on;

% Subplot 2: Refractive Index vs. Normalized Core Radius
subplot(1, 2, 2);
hold on;
for i = 1:length(alpha_values)
    alpha = alpha_values(i);
    n_r = refractive_index(r_a, alpha, n1, delta);
    plot(r_a, n_r, 'LineWidth', 2, 'DisplayName', ['α = ' num2str(alpha)]);
end
xlabel('Normalized Core Radius (r/a)');
ylabel('Refractive Index (n)');
title('Refractive Index vs. Core Radius');
ylim([1.47, 1.51]);  % Set y-axis limits for clarity
legend;
grid on;