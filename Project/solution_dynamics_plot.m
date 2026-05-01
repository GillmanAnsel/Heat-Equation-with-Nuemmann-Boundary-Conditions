%% Solution Dynamics and Error Evolution Visualization
% This script creates a comprehensive 4-panel visualization showing:
% 1. L2 norm decay over time
% 2. Max norm decay over time
% 3. Solution heatmap in space-time
% 4. Error heatmap in space-time

% Add parent directory to path
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) x + t;
u0  = @(x) x; 
f_func = @(x,t) 1;
g_left_func  = @(t) 1;
g_right_func = @(t) 1;

% --- Compute solution with fine mesh for visualization ---
fprintf('Computing fine mesh solution for visualization...\n');
[X_fine, T_fine, U_final_fine, U_hist_fine] = FiniteElementMethod(...
    'L', 1, ...
    'Nx', 256, ...
    'Tf', 5, ...
    'Nt', 500, ...
    'u0', @(x) u_exact(x, 0), ...
    'a_func', @(x) 1, ...
    'f_func', f_func, ...
    'g_left', g_left_func, ...
    'g_right', g_right_func, ...
    'save_plots', false);

% --- Compute L2 norm, max norm, and error over time ---
fprintf('Computing norms and errors...\n');
L2_norm = zeros(length(T_fine), 1);
max_norm = zeros(length(T_fine), 1);
L2_error = zeros(length(T_fine), 1);

for n = 1:length(T_fine)
    dx = X_fine(2) - X_fine(1);
    L2_norm(n) = sqrt(dx * sum(U_hist_fine(:, n).^2));
    max_norm(n) = max(abs(U_hist_fine(:, n)));
    U_exact_n = u_exact(X_fine, T_fine(n))';
    L2_error(n) = sqrt(dx * sum((U_hist_fine(:, n) - U_exact_n).^2));
end

% --- Create a 2x2 subplot figure ---
fprintf('Creating visualization...\n');
fig = figure('Position', [100, 100, 1400, 1000]);

% Plot 1: L2 norm decay (should be exponential)
subplot(2,2,1);
semilogy(T_fine, L2_norm, 'b-', 'LineWidth', 2.5);
hold on;
semilogy(T_fine, 0.33 + T_fine + 0.5*T_fine.^2, 'r--', 'LineWidth', 2, 'DisplayName', '$\frac{1}{3} + t + \frac{t^2}{2}$ (theoretical)');
xlabel('Time t', 'FontSize', 11);
ylabel('$\|u(\cdot,t)\|_{L^2}$', 'FontSize', 11);
title('L² Norm Evolution for $u(x,t) = x + t$', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
legend('FEM Solution', 'Theoretical decay', 'FontSize', 10, 'Interpreter', 'latex');
hold off;

% Plot 2: Max norm decay
subplot(2,2,2);
semilogy(T_fine, max_norm, 'g-', 'LineWidth', 2.5);
hold on;
semilogy(T_fine, 1 + T_fine, 'r--', 'LineWidth', 2, 'DisplayName', '$1 + t$ (theoretical max at x=1)');
xlabel('Time t', 'FontSize', 11);
ylabel('$\max_x |u(x,t)|$', 'FontSize', 11);
title('Max Norm Evolution for $u(x,t) = x + t$', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
legend('FEM Solution', 'Theoretical decay', 'FontSize', 10, 'Interpreter', 'latex');
hold off;

% Plot 3: Heatmap of solution over space-time
subplot(2,2,3);
contourf(T_fine, X_fine, U_hist_fine, 20, 'LineColor', 'none');
colorbar;
xlabel('Time t', 'FontSize', 11);
ylabel('Position x', 'FontSize', 11);
title('Solution $u(x,t) = x + t$ Heatmap', 'FontSize', 12, 'Interpreter', 'latex');
colormap(gca, 'hot');

% Plot 4: Heatmap of error over space-time
subplot(2,2,4);
U_exact_grid = zeros(size(U_hist_fine));
for n = 1:length(T_fine)
    U_exact_grid(:, n) = u_exact(X_fine, T_fine(n))';
end
error_grid = abs(U_hist_fine - U_exact_grid);
contourf(T_fine, X_fine, error_grid, 20, 'LineColor', 'none');
colorbar;
xlabel('Time t', 'FontSize', 11);
ylabel('Position x', 'FontSize', 11);
title('Absolute Error $|u_{FEM} - u_{exact}|$ Heatmap', 'FontSize', 12, 'Interpreter', 'latex');
colormap(gca, 'jet');

sgtitle('Solution Dynamics and Error Evolution for $u(x,t) = x + t$', 'FontSize', 14, 'Interpreter', 'latex');
saveas(gcf, 'solution_dynamics.png');
fprintf('Solution dynamics plot saved to: solution_dynamics.png\n');
close;
