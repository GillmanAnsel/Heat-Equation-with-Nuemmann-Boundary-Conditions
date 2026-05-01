%% Error Over Time Animation
% This script creates an animation showing the L2 error between FEM and exact solution
% evolving over time, visualizing how well the numerical solution tracks the exact solution

% Add parent directory to path to access Project.m
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) x + t;
u0  = @(x) x; 
f_func = @(x,t) 1;
g_left_func  = @(t) 1;
g_right_func = @(t) 1;

% --- Animation parameters ---
Nx = 100;          % Fixed mesh size
Tf = 20;            % Final time
Nt = 100;           % Number of time steps
tau = Tf / (Nt - 1);

% Time steps to capture for animation (evenly spaced)
time_frames = linspace(1, Nt, 12);  % 12 frames from start to end
time_frames = round(time_frames);
time_frames = unique(time_frames);  % Remove duplicates

% GIF output parameters
gif_filename = 'error_over_time_animation.gif';
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

% --- Create figure for animation ---
fig = figure('Position', [100, 100, 1000, 800]);

% --- Pre-compute maximum error across all time steps for fixed z-axis ---
fprintf('\n========== Computing maximum error for fixed axis scaling ==========\n');
max_error_global = 0;
for t_idx = 1:length(T)
    U_exact_t = u_exact(X_grid, T(t_idx))';
    U_approx_t = U_hist(:, t_idx);
    error_t = abs(U_approx_t - U_exact_t);
    max_error_global = max(max_error_global, max(error_t));
end
fprintf('Maximum error across all time steps: %.6e\n', max_error_global);

% --- Create animation frames at selected time steps ---
fprintf('\n========== Creating Error Over Time Animation ==========\n');
fprintf('Creating %d frames showing error evolution...\n\n', length(time_frames));

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
    
    % Compute exact solution at all points and times
    U_exact_current = zeros(length(X_grid), n_plot);
    for j = 1:n_plot
        U_exact_current(:, j) = u_exact(X_grid, T(j))';
    end
    
    % Compute pointwise error
    error_current = abs(U_current - U_exact_current);
    
    % Create meshgrid
    X_grid_col = X_grid(:);  % Column vector
    [T_mesh, X_mesh] = meshgrid(T_current, X_grid_col);
    
    % Plot error surface
    surf(T_mesh, X_mesh, error_current, 'EdgeColor', 'none', 'FaceAlpha', 0.70);
    colormap(gca, 'jet');
    c = colorbar('Location', 'eastoutside');
    ylabel(c, 'Absolute Error |u_{exact} - u_{FEM}|', 'FontSize', 10);
    
    xlabel('Time t', 'FontSize', 11);
    ylabel('Position x', 'FontSize', 11);
    zlabel('Error', 'FontSize', 11);
    title(sprintf('Pointwise Error Evolution for $u(x,t) = x + t$ at t = %.4f (Nx = %d)', T(n), Nx), ...
        'FontSize', 12, 'Interpreter', 'latex');
    view(45, 30);
    zlim([0, max_error_global*1.1]);
    xlim([0, Tf]);
    ylim([0, 1]);
    
    % Compute and display global L2 error at current time
    dx = X_grid(2) - X_grid(1);
    l2_error_current = sqrt(dx * sum(error_current(:, end).^2));
    
    % Add text with L2 error information
    text(0.02, 0.98, sprintf('L^2 Error at t = %.4f: %.6e', T(n), l2_error_current), ...
        'Units', 'normalized', 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 5);
    
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
final_filename = 'error_over_time_final.png';
saveas(gcf, final_filename);
fprintf('Final frame saved to: %s\n', final_filename);

close(fig);
