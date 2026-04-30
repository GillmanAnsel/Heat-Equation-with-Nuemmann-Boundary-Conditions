%% Animation of FEM Solutions for Different Mesh Sizes
% This script creates an animation showing how the FEM solution improves
% as the mesh is refined (h decreases from 1/2 to 1/64)

% Add parent directory to path to access Project.m
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) exp(-t) .* sin(pi * x);
u0  = @(x) sin(pi * x); 
f_func = @(x,t) (pi^2 - 1) * exp(-t) .* sin(pi * x);
g_left_func  = @(t) pi * exp(-t);
g_right_func = @(t) -pi * exp(-t);

% --- Animation parameters ---
h_values = [1/2, 1/4, 1/8, 1/16, 1/32, 1/64];
Nx_values = 1 ./ h_values + 1;
tau = 5e-5;  % Fixed time step
Tf = 5;      % Final time
Nt = ceil(Tf / tau) + 1;

% GIF output parameters
gif_filename = 'surface_plot_animation.gif';
delay_time = 1.5;  % Delay in seconds between frames

% --- Create figure for animation ---
fig = figure('Position', [100, 100, 1000, 800]);

% --- Define time grid ---
T = linspace(0, Tf, Nt);

% --- Compute exact solution on a fine grid (constant across all frames) ---
X_fine = linspace(0, 1, 300);  % Very fine spatial grid for smooth exact solution
U_exact_fine = zeros(length(X_fine), length(T));
for n = 1:length(T)
    U_exact_fine(:, n) = u_exact(X_fine, T(n))';
end
[T_mesh_fine, X_mesh_fine] = meshgrid(T, X_fine);

% --- Loop through mesh sizes and create animation frames ---
fprintf('\n========== Creating Animation ==========\n');
fprintf('Animation will show 6 frames as mesh is refined (h: 1/2 → 1/64)\n\n');

for k = 1:length(h_values)
    h = h_values(k);
    Nx = Nx_values(k);
    
    fprintf('Frame %d/%d: Computing solution for h = 1/%d (Nx = %d)...\n', ...
        k, length(h_values), round(1/h), Nx);
    
    % Run FEM solver
    [X_grid, T, U_final, U_hist] = FiniteElementMethod(...
        'L', 1, ...
        'Nx', Nx, ...
        'Tf', Tf, ...
        'Nt', Nt, ...
        'u0', @(x) u_exact(x, 0), ...
        'a_func', @(x) 1, ...
        'f_func', f_func, ...
        'g_left', g_left_func, ...
        'g_right', g_right_func, ...
        'save_plots', false);
    
    % Clear previous plot
    clf;
    
    % Create 3D surface plot for final time

    [T_mesh, X_mesh] = meshgrid(T, X_grid);
    surf(T_mesh, X_mesh, U_hist, 'EdgeColor', 'none', 'FaceAlpha', 0.70);
    colormap(gca, 'hot');
    colorbar;

    % Hold on and plot exact solution as very smooth translucent overlay
    hold on;
    surf(T_mesh_fine, X_mesh_fine, U_exact_fine, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    colormap(gca, 'cool');
    hold off;

    xlabel('Time t', 'FontSize', 11);
    ylabel('Position x', 'FontSize', 11);
    zlabel('u(x, t)', 'FontSize', 11);
    title(sprintf('FEM Solution converging: $u_{exact} = e^{-t}\\sin(\\pi x)$, h = 1/%d (Nx = %d) | Exact (red)', round(1/h), Nx), 'FontSize', 12, 'Interpreter', 'latex');
    view(45, 30);
    zlim([-1, 1]);
    
    % Create legend with visual transparency difference
    p1 = patch(nan, nan, nan, 'FaceColor', [0 0 1], 'FaceAlpha', 0.70);
    p2 = patch(nan, nan, nan, 'FaceColor', [0 0 1], 'FaceAlpha', 0.25);
    legend([p1, p2], 'FEM Solution (opaque)', 'Exact Solution (translucent)', 'FontSize', 10);
    
    % Capture frame and save to GIF
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF (append mode for subsequent frames)
    if k == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        fprintf('  GIF started: %s\n', gif_filename);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    
    % Pause to show frame
    pause(delay_time);
end

fprintf('\nAnimation complete!\n');
fprintf('GIF saved to: %s\n', gif_filename);
close(fig);
