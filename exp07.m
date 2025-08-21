% Wavelength range (0.8 to 1.8 μm)
wavelength = linspace(0.8, 1.8, 1000);

% Loss calculations
A_R = 0.9;
rayleigh_loss = A_R ./ (wavelength .^ 4);

A_UV = 0.1;
B_UV = 5.0;
uv_loss = A_UV * exp(-B_UV * wavelength);

A_IR = 0.01;
B_IR = 7.0;
lambda_0 = 1.27;
ir_loss = A_IR * exp(B_IR * (wavelength - lambda_0));

total_loss = rayleigh_loss + uv_loss + ir_loss;

% Find minimum loss point
[min_loss, min_loss_idx] = min(total_loss);
min_wavelength = wavelength(min_loss_idx);

fprintf('Minimum Loss: %.2f dB/km at Wavelength: %.3f nm\n', min_loss, min_wavelength * 1000);

% Create the plot (normal scale)
figure('Position', [100, 100, 1200, 800]);

% Plot individual loss components (linear scale)
plot(wavelength, rayleigh_loss, 'b--', 'LineWidth', 2, 'DisplayName', 'Rayleigh Scattering');
hold on;
plot(wavelength, uv_loss, 'g--', 'LineWidth', 2, 'DisplayName', 'UV Absorption');
plot(wavelength, ir_loss, 'm--', 'LineWidth', 2, 'DisplayName', 'IR Absorption');
plot(wavelength, total_loss, 'r-', 'LineWidth', 3, 'DisplayName', 'Total Loss');

% Mark minimum loss point
plot(min_wavelength, min_loss, 'ro', 'MarkerSize', 8, 'DisplayName', sprintf('Min Loss: %.2f dB/km @ %.3f nm', min_loss, min_wavelength * 1000));

% Vertical line at minimum loss
xline(min_wavelength, 'gray--', 'LineWidth', 2);

% Formatting
xlabel('Wavelength (μm)', 'FontSize', 12);
ylabel('Loss (dB/km)', 'FontSize', 12);
title('Silica Optical Fiber Loss vs Wavelength (Linear Scale)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('FontSize', 10);
xlim([0.8, 1.8]);
ylim([0, max(total_loss) * 1.1]);