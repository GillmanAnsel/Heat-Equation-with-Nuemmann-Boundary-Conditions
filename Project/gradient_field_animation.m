%% Gradient Field Animation for Heat Flux
% This script creates an animation showing the spatial gradient (heat flux)
% du/dx over time, visualizing where heat is flowing most rapidly

% Add parent directory to path to access Project.m
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) exp(-t) .* sin(pi * x);
u0  = @(x) sin(pi * x); 
f_func = @(x,t) (pi^2 - 1) * exp(-t) .* sin(pi * x);
g_left_func  = @(t) pi * exp(-t);
g_right_func = @(t) -pi * exp(-t);

% Exact gradient (du/dx)
du_exact = @(x,t) pi * exp(-t) .* cos(pi * x);

% --- Animation parameters ---
Nx = 100;          % Fixed mesh size
Tf = 1;            % Final time
Nt = 50;           % Number of time steps
tau = Tf / (Nt - 1);

% Time steps to capture for animation (evenly spaced)
time_frames = linspace(1, Nt, 12);  % 12 frames from start to end
time_frames = round(time_frames);
time_frames = unique(time_frames);  % Remove duplicates

% GIF output parameters
gif_filename = 'gradient_field_animation.gif';
delay_time = 0.5;  % Delay in seconds between frames

% --- Run FEM solver once to get full solution ---
fprintf('\n========== Computing Full Solution ==========\n');
fprintf('Computing FEM solution for Nx = %d, Nt = %d...\n', Nx, Nt);

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

% Compute exact gradient on a fine grid for comparison
X_fine = linspace(0, 1, 300);
dU_exact_fine_all = zeros(length(X_fine), length(T));
for n = 1:length(T)
    dU_exact_fine_all(:, n) = du_exact(X_fine, T(n))';
end

% --- Create figure for animation ---
fig = figure('Position', [100, 100, 1000, 800]);

% --- Create animation frames at selected time steps ---
fprintf('\n========== Creating Gradient Field Animation ==========\n');
fprintf('Creating %d frames showing heat flux evolution...\n\n', length(time_frames));

for frame_idx = 1:length(time_frames)
    n = time_frames(frame_idx);
    
    fprintf('Frame %d/%d: t = %.4f\n', frame_idx, length(time_frames), T(n));
    
    % Clear previous plot
    clf;
    
    % Ensure we have at least 2 columns for surf
    if n == 1
        n_plot = 2;
    else
        n_plot = n;
    end
    
    % Extract data up to current time step
    T_current = T(1:n_plot);
    U_current = U_hist(:, 1:n_plot);
    
    % Compute FEM gradient (spatial derivative) numerically
    dU_current = zeros(size(U_current));
    for j = 1:n_plot
        % Use central differences in space (except at boundaries)
        dU_current(1, j) = (U_current(2, j) - U_current(1, j)) / (X_grid(2) - X_grid(1));
        dU_current(end, j) = (U_current(end, j) - U_current(end-1, j)) / (X_grid(end) - X_grid(end-1));
        for i = 2:(length(X_grid)-1)
            dU_current(i, j) = (U_current(i+1, j) - U_current(i-1, j)) / (X_grid(i+1) - X_grid(i-1));
        end
    end
    
    % Create meshgrid
    X_grid_col = X_grid(:);  % Column vector
    [T_mesh, X_mesh] = meshgrid(T_current, X_grid_col);
    
    % Exact gradient up to current time
    dU_exact_current = zeros(length(X_grid_col), n_plot);
    for j = 1:n_plot
        dU_exact_current(:, j) = du_exact(X_grid_col, T(j))';
    end
    
    % Plot FEM gradient field
    surf(T_mesh, X_mesh, dU_current, 'EdgeColor', 'none', 'FaceAlpha', 0.70);
    colormap(gca, 'cool');
    colorbar('Location', 'southoutside');
    
    % Hold on and plot exact gradient as translucent overlay
    hold on;
    X_fine_col = X_fine(:);  % Column vector
    [T_mesh_fine, X_mesh_fine] = meshgrid(T_current, X_fine_col);
    dU_exact_fine_current = dU_exact_fine_all(:, 1:n_plot);
    surf(T_mesh_fine, X_mesh_fine, dU_exact_fine_current, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    hold off;
    
    xlabel('Time t', 'FontSize', 11);
    ylabel('Position x', 'FontSize', 11);
    zlabel('du/dx (Heat Flux)', 'FontSize', 11);
    title(sprintf('Heat Flux Field: $\\frac{\\partial u}{\\partial x} = \\pi e^{-t}\\cos(\\pi x)$ at t = %.4f (Nx = %d)', T(n), Nx), ...
        'FontSize', 12, 'Interpreter', 'latex');
    view(45, 45);
    zlim([-4, 4]);
    xlim([0, Tf]);
    ylim([0, 1]);
    
    % Add Neumann boundary condition annotations
    t_current = T(n);
    bc_left = g_left_func(t_current);
    bc_right = g_right_func(t_current);
    
    % Add text labels for Neumann BCs at boundaries
    text(t_current, 0, bc_left + 0.5, sprintf('BC left: du/dx = %.2f', bc_left), ...
        'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    text(t_current, 1, bc_right + 0.5, sprintf('BC right: du/dx = %.2f', bc_right), ...
        'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Create legend with visual transparency difference
    p1 = patch(nan, nan, nan, 'FaceColor', [0 0 1], 'FaceAlpha', 0.70);
    p2 = patch(nan, nan, nan, 'FaceColor', [0 0 1], 'FaceAlpha', 0.25);
    legend([p1, p2], 'FEM Gradient (opaque)', 'Exact Gradient (translucent)', 'FontSize', 10);
    
    % Capture frame and save to GIF
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF (append mode for subsequent frames)
    if frame_idx == 1
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

% Save final frame as a high-resolution figure
final_filename = 'gradient_field_final.png';
saveas(gcf, final_filename);
fprintf('Final frame saved to: %s\n', final_filename);

close(fig);
