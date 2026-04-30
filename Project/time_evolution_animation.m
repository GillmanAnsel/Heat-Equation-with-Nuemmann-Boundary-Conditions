%% Time Evolution Animation for FEM Solution
% This script creates an animation showing how the FEM solution evolves over time
% for a fixed mesh, visualizing the heat diffusion process

% Add parent directory to path to access Project.m
addpath(pwd);

% --- Manufactured Solution (same as in Project.m) ---
u_exact = @(x,t) exp(-t) .* sin(pi * x);
u0  = @(x) sin(pi * x); 
f_func = @(x,t) (pi^2 - 1) * exp(-t) .* sin(pi * x);
g_left_func  = @(t) pi * exp(-t);
g_right_func = @(t) -pi * exp(-t);

% --- Animation parameters ---
Nx = 100;          % Fixed mesh size
Tf = 5;            % Final time
Nt = 50;           % Number of time steps (will create frames at selected time steps)
tau = Tf / (Nt - 1);

% Time steps to capture for animation (evenly spaced)
time_frames = linspace(1, Nt, 12);  % 12 frames from start to end
time_frames = round(time_frames);
time_frames = unique(time_frames);  % Remove duplicates

% GIF output parameters
gif_filename = 'time_evolution_animation.gif';
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

% --- Create animation frames at selected time steps ---
fprintf('\n========== Creating Time Evolution Animation ==========\n');
fprintf('Creating %d frames showing time evolution...\n\n', length(time_frames));

for frame_idx = 1:length(time_frames)
    n = time_frames(frame_idx);
    
    fprintf('Frame %d/%d: t = %.4f\n', frame_idx, length(time_frames), T(n));
    
    % Clear previous plot
    clf;
    
    % surf requires at least 2 time points — pad to 2 columns if on first frame
    if n == 1
        n_plot = 2;
    else
        n_plot = n;
    end

    % Extract data up to current time step
    T_current = T(1:n_plot);
    U_current = U_hist(:, 1:n_plot);
    
    % Ensure X_grid and T_current are properly formatted for meshgrid
    X_grid_col = X_grid(:);  % Column vector
    [T_mesh, X_mesh] = meshgrid(T_current, X_grid_col);
    
  
    
    % Plot FEM solution
    surf(T_mesh, X_mesh, U_current, 'EdgeColor', 'none', 'FaceAlpha', 0.70);
    colormap(gca, 'hot');
    colorbar;
    
    
    xlabel('Time t', 'FontSize', 11);
    ylabel('Position x', 'FontSize', 11);
    zlabel('u(x, t)', 'FontSize', 11);
    title(sprintf('Time Evolution: $u_{exact} = e^{-t}\\sin(\\pi x)$ at t = %.4f (Nx = %d)', T(n), Nx), ...
        'FontSize', 12, 'Interpreter', 'latex');
    view(45, 30);
    zlim([-1, 1]);
    xlim([0, Tf]);
    ylim([0, 1]);
    
    
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
close(fig);