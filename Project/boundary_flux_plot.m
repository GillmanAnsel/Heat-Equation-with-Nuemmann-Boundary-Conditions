%% Boundary Flux Over Time
% This script creates a plot showing how the prescribed Neumann boundary conditions
% (heat flux) evolve over time at both boundaries (x=0 and x=L)

% Add parent directory to path to access Project.m
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) x + t;
u0  = @(x) x; 
f_func = @(x,t) 1;
g_left_func  = @(t) 1;       % Left Neumann BC: du/dx at x=0
g_right_func = @(t) 1;      % Right Neumann BC: du/dx at x=L

% --- Time parameters ---
Tf = 8;            % Final time
Nt = 200;          % High resolution for smooth plot
T = linspace(0, Tf, Nt);

% --- Compute boundary fluxes over time ---
flux_left = g_left_func(T);
flux_right = g_right_func(T);

% --- Create figure ---
fig = figure('Position', [100, 100, 1000, 600]);

% Plot boundary fluxes
plot(T, flux_left, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Left BC: $g_{left}(t) = 1$ at $x=0$');
hold on;
plot(T, flux_right, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Right BC: $g_{right}(t) = 1$ at $x=L$');
plot(T, zeros(size(T)), 'k--', 'LineWidth', 1, 'DisplayName', 'Zero flux line');
hold off;

xlabel('Time t', 'FontSize', 12);
ylabel('Heat Flux: $\frac{\partial u}{\partial x}$ (Neumann BC value)', 'FontSize', 12);
title('Boundary Heat Flux Over Time for $u(x,t) = x + t$', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
legend('FontSize', 11, 'Location', 'best', 'Interpreter', 'latex');
xlim([0, Tf]);

% Add annotations for both boundaries at t=0 and t=Tf
% Left boundary (blue)
text(0.05, flux_left(1)*0.85, sprintf('Left init: %.3f', flux_left(1)), ...
    'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold', 'BackgroundColor', 'white');
text(Tf-0.5, flux_left(end)*1.15 +0.25, sprintf('Left final: %.2e', flux_left(end)), ...
    'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold', 'BackgroundColor', 'white');

% Right boundary (red)
text(0.05, flux_right(1)*0.85, sprintf('Right init: %.3f', flux_right(1)), ...
    'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold', 'BackgroundColor', 'white');
text(Tf-0.35, flux_right(end)*0.85 -0.25, sprintf('Right final: %.2e', flux_right(end)), ...
    'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold', 'BackgroundColor', 'white');

% Save figure
saveas(gcf, 'boundary_flux_plot.png');
fprintf('Boundary flux plot saved to: boundary_flux_plot.png\n');

% Create a second figure showing flux magnitude (absolute value)
fig2 = figure('Position', [100, 100, 1000, 600]);

flux_mag_left = abs(flux_left);
flux_mag_right = abs(flux_right);

plot(T, flux_mag_left, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Left BC magnitude: $|g_{left}(t)|$');
hold on;
plot(T, flux_mag_right, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Right BC magnitude: $|g_{right}(t)|$');
hold off;

xlabel('Time t', 'FontSize', 12);
ylabel('Absolute Heat Flux Magnitude', 'FontSize', 12);
title('Boundary Heat Flux Magnitude for $u(x,t) = x + t$ Over Time', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
legend('FontSize', 11, 'Location', 'best', 'Interpreter', 'latex');
xlim([0, Tf]);
ylim([0, max(flux_mag_left)*1.1]);

% Save second figure
saveas(gcf, 'boundary_flux_magnitude_plot.png');
fprintf('Boundary flux magnitude plot saved to: boundary_flux_magnitude_plot.png\n');

% Print summary statistics
fprintf('\n========== Boundary Flux Analysis ==========\n');
fprintf('Left BC (x=0):\n');
fprintf('  Initial flux: g_left(0) = %.6f\n', flux_left(1));
fprintf('  Final flux:   g_left(1) = %.6e\n', flux_left(end));
fprintf('  Decay factor: g_left(1) / g_left(0) = %.6e\n', flux_left(end) / flux_left(1));

fprintf('\nRight BC (x=L):\n');
fprintf('  Initial flux: g_right(0) = %.6f\n', flux_right(1));
fprintf('  Final flux:   g_right(1) = %.6e\n', flux_right(end));
fprintf('  Decay factor: g_right(1) / g_right(0) = %.6e\n', flux_right(end) / flux_right(1));


close all;
