%% Spatial vs Temporal Error Analysis
% This script decomposes the total error into spatial and temporal contributions
% by solving with different levels of refinement in space and time

% Add parent directory to path
addpath(pwd);

% --- Manufactured Solution ---
u_exact = @(x,t) exp(-t) .* sin(pi * x);
u0  = @(x) sin(pi * x); 
f_func = @(x,t) (pi^2 - 1) * exp(-t) .* sin(pi * x);
g_left_func  = @(t) pi * exp(-t);
g_right_func = @(t) -pi * exp(-t);

% --- Reference solution (fine grid) ---
fprintf('Computing reference solution (fine spatial & temporal grids)...\n');
[X_ref, T_ref, ~, U_ref] = FiniteElementMethod(...
    'L', 1, ...
    'Nx', 128, ...
    'Tf', 5, ...
    'Nt', 1000, ...
    'u0', @(x) u_exact(x, 0), ...
    'a_func', @(x) 1, ...
    'f_func', f_func, ...
    'g_left', g_left_func, ...
    'g_right', g_right_func, ...
    'save_plots', false);

% --- Define coarse and fine grid/time parameters (balanced error: Nt = 5*Nx²) ---
% Coarse (balanced): Nx=7, Nt=5*49=245
% Fine (balanced): Nx=15, Nt=5*225=1125
% Solution 1: Nx_coarse=7, Nt_fine=1125 (tests spatial error)
% Solution 2: Nx_fine=15, Nt_coarse=245 (tests temporal error)
% Solution 3: Nx_coarse=7, Nt_coarse=245 (both balanced)
Nx_coarse = 7;
Nx_fine = 15;
Nt_coarse = 245;      % 5*7^2
Nt_fine = 1125;       % 5*15^2

% --- Solution 1: Coarse spatial, fine temporal (spatial error dominates) ---
fprintf('Computing coarse spatial / fine temporal solution...\n');
[X_coarse_fine_t, T_coarse_fine_t, ~, U_coarse_fine_t] = FiniteElementMethod(...
    'L', 1, ...
    'Nx', Nx_coarse, ...
    'Tf', 5, ...
    'Nt', Nt_fine, ...
    'u0', @(x) u_exact(x, 0), ...
    'a_func', @(x) 1, ...
    'f_func', f_func, ...
    'g_left', g_left_func, ...
    'g_right', g_right_func, ...
    'save_plots', false);

% --- Solution 2: Fine spatial, coarse temporal (temporal error dominates) ---
fprintf('Computing fine spatial / coarse temporal solution...\n');
[X_fine_coarse_t, T_fine_coarse_t, ~, U_fine_coarse_t] = FiniteElementMethod(...
    'L', 1, ...
    'Nx', Nx_fine, ...
    'Tf', 5, ...
    'Nt', Nt_coarse, ...
    'u0', @(x) u_exact(x, 0), ...
    'a_func', @(x) 1, ...
    'f_func', f_func, ...
    'g_left', g_left_func, ...
    'g_right', g_right_func, ...
    'save_plots', false);

% --- Solution 3: Coarse spatial and temporal (both errors) ---
fprintf('Computing coarse spatial / coarse temporal solution...\n');
[X_coarse, T_coarse, ~, U_coarse_coarse] = FiniteElementMethod(...
    'L', 1, ...
    'Nx', Nx_coarse, ...
    'Tf', 5, ...
    'Nt', Nt_coarse, ...
    'u0', @(x) u_exact(x, 0), ...
    'a_func', @(x) 1, ...
    'f_func', f_func, ...
    'g_left', g_left_func, ...
    'g_right', g_right_func, ...
    'save_plots', false);

% --- Interpolate all solutions to reference grid for comparison ---
fprintf('Interpolating solutions to reference grid...\n');

% Create meshgrids for interpolation
[T_coarse_fine_t_grid, X_coarse_fine_t_grid] = meshgrid(T_coarse_fine_t, X_coarse_fine_t);
[T_fine_coarse_t_grid, X_fine_coarse_t_grid] = meshgrid(T_fine_coarse_t, X_fine_coarse_t);
[T_coarse_grid, X_coarse_grid] = meshgrid(T_coarse, X_coarse);
[T_ref_grid, X_ref_grid] = meshgrid(T_ref, X_ref);

% Coarse spatial / fine temporal
U_coarse_fine_t_interp = interp2(T_coarse_fine_t_grid, X_coarse_fine_t_grid, U_coarse_fine_t, T_ref_grid, X_ref_grid, 'linear');

% Fine spatial / coarse temporal
U_fine_coarse_t_interp = interp2(T_fine_coarse_t_grid, X_fine_coarse_t_grid, U_fine_coarse_t, T_ref_grid, X_ref_grid, 'linear');

% Coarse / coarse
U_coarse_coarse_interp = interp2(T_coarse_grid, X_coarse_grid, U_coarse_coarse, T_ref_grid, X_ref_grid, 'linear');

% --- Compute errors ---
fprintf('Computing errors...\n');
dx_ref = X_ref(2) - X_ref(1);

% Spatial error (from coarse spatial / fine temporal vs reference)
spatial_error = abs(U_coarse_fine_t_interp - U_ref);

% Temporal error (from fine spatial / coarse temporal vs reference)
temporal_error = abs(U_fine_coarse_t_interp - U_ref);

% Total error (from coarse / coarse vs reference)
total_error = abs(U_coarse_coarse_interp - U_ref);

% Compute L2 errors over time
spatial_L2 = zeros(length(T_ref), 1);
temporal_L2 = zeros(length(T_ref), 1);
total_L2 = zeros(length(T_ref), 1);

for n = 1:length(T_ref)
    spatial_L2(n) = sqrt(dx_ref * sum(spatial_error(:, n).^2));
    temporal_L2(n) = sqrt(dx_ref * sum(temporal_error(:, n).^2));
    total_L2(n) = sqrt(dx_ref * sum(total_error(:, n).^2));
end

% --- Create visualization ---
fprintf('Creating visualizations...\n');
fig = figure('Position', [100, 100, 1400, 1000]);

% Plot 1: L2 errors over time (semilog)
subplot(2,2,1);
semilogy(T_ref, spatial_L2, 'b-', 'LineWidth', 2.5, 'DisplayName', sprintf('Spatial error (Nx=%d, Nt=%d)', Nx_coarse, Nt_fine));
hold on;
semilogy(T_ref, temporal_L2, 'r-', 'LineWidth', 2.5, 'DisplayName', sprintf('Temporal error (Nx=%d, Nt=%d)', Nx_fine, Nt_coarse));
semilogy(T_ref, total_L2, 'k-', 'LineWidth', 2.5, 'DisplayName', sprintf('Total error (Nx=%d, Nt=%d)', Nx_coarse, Nt_coarse));
xlabel('Time t', 'FontSize', 11);
ylabel('L² Error', 'FontSize', 11);
title('Error Components Over Time', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
legend('FontSize', 10, 'Location', 'best');
hold off;
set(gca, 'YScale', 'log');

% Plot 2: Error contribution comparison at different times
subplot(2,2,2);
time_indices = [1, round(length(T_ref)*0.25), round(length(T_ref)*0.5), round(length(T_ref)*0.75), length(T_ref)];
time_points = T_ref(time_indices);
spatial_vals = spatial_L2(time_indices);
temporal_vals = temporal_L2(time_indices);

x_pos = 1:length(time_points);
width = 0.35;
bar(x_pos - width/2, spatial_vals, width, 'b', 'DisplayName', 'Spatial error');
hold on;
bar(x_pos + width/2, temporal_vals, width, 'r', 'DisplayName', 'Temporal error');
set(gca, 'YScale', 'log');
set(gca, 'XTick', x_pos);
set(gca, 'XTickLabel', arrayfun(@(t) sprintf('t=%.1f', t), time_points, 'UniformOutput', false));
ylabel('L² Error', 'FontSize', 11);
title('Spatial vs Temporal Error at Key Time Points', 'FontSize', 12, 'Interpreter', 'latex');
legend('FontSize', 10);
grid on;
hold off;

% Plot 3: Spatial error heatmap
subplot(2,2,3);
contourf(T_ref, X_ref, spatial_error, 20, 'LineColor', 'none');
colorbar;
xlabel('Time t', 'FontSize', 11);
ylabel('Position x', 'FontSize', 11);
title('Spatial Error: $|u_{coarse-spatial} - u_{ref}|$', 'FontSize', 12, 'Interpreter', 'latex');
colormap(gca, 'hot');

% Plot 4: Temporal error heatmap
subplot(2,2,4);
contourf(T_ref, X_ref, temporal_error, 20, 'LineColor', 'none');
colorbar;
xlabel('Time t', 'FontSize', 11);
ylabel('Position x', 'FontSize', 11);
title('Temporal Error: $|u_{coarse-temporal} - u_{ref}|$', 'FontSize', 12, 'Interpreter', 'latex');
colormap(gca, 'cool');

sgtitle('Spatial vs Temporal Error Decomposition', 'FontSize', 14, 'Interpreter', 'latex');
saveas(gcf, 'spatial_temporal_error_analysis.png');
fprintf('Spatial vs temporal error plot saved to: spatial_temporal_error_analysis.png\n');
close;

% --- Print summary statistics ---
fprintf('\n========== Error Analysis Summary ==========\n');
fprintf('Reference solution: Nx=%d, Nt=%d\n', length(X_ref), length(T_ref));
fprintf('\nBalanced Error Analysis (Nt ≈ 5*Nx² for optimal error balance):\n');
fprintf('\nSpatial error (Nx=%d, Nt=%d - balanced): \n', Nx_coarse, Nt_fine);
fprintf('  h = %.4f, τ = %.6f, ratio τ/h² = %.4f\n', 1/Nx_coarse, 5/Nt_fine, (5/Nt_fine)/((1/Nx_coarse)^2));
fprintf('  Initial: %.6e\n', spatial_L2(1));
fprintf('  Final:   %.6e\n', spatial_L2(end));
fprintf('  Max:     %.6e\n', max(spatial_L2));

fprintf('\nTemporal error (Nx=%d, Nt=%d - coarse time): \n', Nx_fine, Nt_coarse);
fprintf('  h = %.4f, τ = %.6f, ratio τ/h² = %.4f\n', 1/Nx_fine, 5/Nt_coarse, (5/Nt_coarse)/((1/Nx_fine)^2));
fprintf('  Initial: %.6e\n', temporal_L2(1));
fprintf('  Final:   %.6e\n', temporal_L2(end));
fprintf('  Max:     %.6e\n', max(temporal_L2));

fprintf('\nTotal error (Nx=%d, Nt=%d - both coarse):\n', Nx_coarse, Nt_coarse);
fprintf('  Initial: %.6e\n', total_L2(1));
fprintf('  Final:   %.6e\n', total_L2(end));
fprintf('  Max:     %.6e\n', max(total_L2));

fprintf('\nDominant error source at final time:\n');
[~, dominant] = max([spatial_L2(end), temporal_L2(end)]);
if dominant == 1
    fprintf('  SPATIAL ERROR dominates (ratio: %.2f)\n', spatial_L2(end) / (temporal_L2(end) + eps));
    fprintf('  Recommendation: Increase Nx (refine spatial mesh)\n');
else
    fprintf('  TEMPORAL ERROR dominates (ratio: %.2f)\n', temporal_L2(end) / (spatial_L2(end) + eps));
    fprintf('  Recommendation: Increase Nt or use smaller time step\n');
end
