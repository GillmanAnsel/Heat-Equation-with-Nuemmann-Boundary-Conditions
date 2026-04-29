% filepath: /workspaces/Heat-Equation-with-Nuemmann-Boundary-Conditions/Project/Project_Copy.m
%% ----------APDE Project----------
%%-----Task 1-----

%% STEP 1: Main FEM Solver Function (Heat Equation with Neumann BCs)
function [X_grid, T, U_final, U_history] = FiniteElementMethod(varargin)
 
    % --- Parse input parameters ---
    p = inputParser;
    addParameter(p, 'L', 1);           % Domain length
    addParameter(p, 'Nx', 20);         % Number of spatial nodes
    addParameter(p, 'Tf', 1);          % Final time
    addParameter(p, 'Nt', 20);         % Number of time steps
    addParameter(p, 'a_func', @(x) 1); % Diffusion coefficient
    addParameter(p, 'f_func', @(x,t) 0); % Source term
    addParameter(p, 'g_left', @(t) 0);   % Left Neumann BC
    addParameter(p, 'g_right', @(t) 0);  % Right Neumann BC
    addParameter(p, 'u0', @(x) 0);       % Initial condition
    addParameter(p, 'save_plots', false);% Save plots flag
    parse(p, varargin{:});
    params = p.Results;

    % --- Extract parameters for use ---
    L = params.L;
    Nx = params.Nx;
    Tf = params.Tf;
    Nt = params.Nt;
    a_func = params.a_func;
    f_func = params.f_func;
    g_left = params.g_left;
    g_right = params.g_right;
    u0 = params.u0;
    save_plots = params.save_plots;

    % --- Set up spatial and temporal grids ---
    X_grid = linspace(0, L, Nx);    % Spatial grid
    T = linspace(0, Tf, Nt);        % Time grid
    dt = T(2) - T(1);               % Time step size

    % --- Assemble FEM matrices ---
    K = StiffnessAssembler1D(X_grid, a_func, [0,0]); % Stiffness matrix
    M = MassAssembler1D(X_grid);                     % Mass matrix
    A_system = (1/dt)*M + K;                         % System matrix

    % --- Initial condition ---
    U_prev = u0(X_grid)';            % Initial solution vector
    U_history = zeros(Nx, Nt);       % Store solution at all time steps
    U_history(:, 1) = U_prev;

    % --- Time-stepping loop (implicit Euler) ---
    for n = 2:length(T)
        t_n = T(n);
        F = LoadAssembler1D(X_grid, f_func, t_n);    % Load vector

        % --- Neumann BCs: add boundary fluxes to RHS ---
        B = zeros(Nx, 1);
        h = X_grid(2) - X_grid(1);
        B(1) = -a_func(0) * g_left(t_n);     % Left boundary
        B(Nx) = a_func(L) * g_right(t_n);    % Right boundary
        % BC Interpretation:
        % g_left:  Heat flux at x=0. Positive value means heat leaves domain at left.
        % g_right: Heat flux at x=L. Positive value means heat leaves domain at right.

        % --- Solve linear system for new time step ---
        b = (1/dt)*M*U_prev + F + B;
        U_current = A_system \ b;
        U_history(:, n) = U_current;
        U_prev = U_current;
    end

    % --- Output final solution and optionally plot ---
    U_final = U_history(:, end);
    if save_plots
        plot_fem_results(X_grid, T, U_history, U_final);
    end
end

%% STEP 2: Plotting Functions
function plot_fem_results(X_grid, T, U_history, U_final)
    % Create output directory if needed
    if ~exist('Project', 'dir')
        mkdir('Project');
    end
    % Plot final solution at last time
    figure;
    plot(X_grid, U_final, 'b-', 'LineWidth', 2);
    xlabel('x');
    ylabel('u(x, T)');
    title(sprintf('FEM Solution of Heat Equation at t = %.2f', T(end)));
    grid on;
    saveas(gcf, fullfile('Project', 'solution.png'));
    close;

    % Plot solution evolution as a heatmap
    figure;
    imagesc(T, X_grid, U_history);
    colorbar;
    xlabel('Time t');
    ylabel('Position x');
    title('Heat Equation Solution Over Time');
    saveas(gcf, fullfile('Project', 'solution_heatmap.png'));
    close;
end

%% STEP 3: Basis Function (Hat/Tent Function)
function phi = hat_function(X_grid, k)
    % Returns the k-th hat function evaluated on X_grid
    h = X_grid(2) - X_grid(1);
    phi = zeros(size(X_grid));
    if k == 1
        xk = X_grid(k);
        xk_next = X_grid(k+1);
        idx1 = (X_grid >= xk & X_grid <= xk_next);
        phi(idx1) = (xk_next - X_grid(idx1)) / h;
    elseif k == length(X_grid)
        xk = X_grid(k);
        xk_prev = X_grid(k-1);
        idx2 = (X_grid >= xk_prev & X_grid <= xk);
        phi(idx2) = (X_grid(idx2) - xk_prev) / h;
    else
        xk_prev = X_grid(k-1);
        xk = X_grid(k);
        xk_next = X_grid(k+1);
        idx1 = (X_grid >= xk_prev & X_grid < xk);
        phi(idx1) = (X_grid(idx1) - xk_prev) / h;
        idx2 = (X_grid >= xk & X_grid <= xk_next);
        phi(idx2) = (xk_next - X_grid(idx2)) / h;
    end
end
% Plotting some hat functions for visualization
figure;
hold on;
for k = 1:length(X_grid)
    phi_k = hat_function(X_grid, k);
    plot(X_grid, phi_k + (k-1)*0.5, 'LineWidth', 1.5); % Shift up for visibility
end
xlabel('x');
ylabel('Hat Functions');
title('Hat (Tent) Functions for 1D FEM');
grid on;
saveas(gcf, fullfile('Project', 'hat_functions.png'));
close;  
%% STEP 4: Stiffness Matrix Assembly
function A = StiffnessAssembler1D(x, a, kappa)
    % Assemble stiffness matrix for 1D FEM
    n = length(x)-1;
    A = zeros(n+1, n+1);
    for i = 1:n
        h = x(i+1) - x(i);
        xmid = (x(i+1) + x(i))/2;
        amid = a(xmid);
        A(i,i) = A(i,i) + amid/h;
        A(i,i+1) = A(i,i+1) - amid/h;
        A(i+1,i) = A(i+1,i) - amid/h;
        A(i+1,i+1) = A(i+1,i+1) + amid/h;
    end
    A(1,1) = A(1,1) + kappa(1);
    A(n+1,n+1) = A(n+1,n+1) + kappa(2);
end

%% STEP 5: Mass Matrix Assembly
function M = MassAssembler1D(x)
    % Assemble mass matrix for 1D FEM
    n = length(x)-1;
    M = zeros(n+1, n+1);
    h = x(2) - x(1);
    M(1,1) = h/3;
    for i = 2:n
        M(i,i) = 2*h/3;
    end
    M(n+1,n+1) = h/3;
    for i = 1:n
        M(i,i+1) = h/6;
        M(i+1,i) = h/6;
    end
end

%% STEP 6: Load Vector Assembly
function F = LoadAssembler1D(x, f_func, t)
    % Assemble load vector for 1D FEM
    n = length(x)-1;
    F = zeros(n+1, 1);
    h = x(2) - x(1);
    for i = 1:n
        xmid = (x(i) + x(i+1)) / 2;
        f_mid = f_func(xmid, t);
        F(i) = F(i) + f_mid * h / 2;
        F(i+1) = F(i+1) + f_mid * h / 2;
    end
end

%% STEP 7: Example Run (Test the solver)

u_exact = @(x,t) x + t;
f_func = @(x,t) 1; 
g_left_func = @(t) 1;  
g_right_func = @(t) 1; 

X_grid = linspace(0, 1, 21);
[X, T, U_final, U_hist] = FiniteElementMethod(...
    'Nx', 50, ... % Number of spatial nodes
    'Nt', 100, ... % Number of time steps
    'Tf', 1, ... % Final time
    'u0', @(x) u_exact(x, 0), ... % Initial condition
    'a_func', @(x) 1, ... % Diffusion coefficient
    'f_func', f_func, ... % f(x,t)
    'g_left', g_left_func, ... % Left Neumann BC
    'g_right', g_right_func, ...  % Right Neumann BC
    'save_plots', true);


disp('Task 1 complete!!!!!!');



%%-----Task 2-----
%% Method of Manufactured Solutions

disp('Task 2 complete!!!!!!');

%%-----Task 3-----
%% Spatial Convergence Study using Method of Manufactured Solutions



% Fixed parameters
tau = 5e-3;  % Fixed time step
Tf = 1;    % Final time
Nt = ceil(Tf / tau) + 1;  % Number of time steps

% Mesh sizes (h values)
h_values = [1/2, 1/4, 1/8, 1/16, 1/32, 1/64];
Nx_values = 1 ./ h_values + 1;  % Number of nodes corresponding to each h

% Storage for errors and convergence rates
E_h = zeros(size(h_values));
convergence_rates = zeros(size(h_values));

% Convergence study loop
fprintf('\n========== TASK 3: Spatial Convergence Study ==========\n');
fprintf('Fixed time step: tau = %.2e\n', tau);
fprintf('Final time: Tf = %.2f\n\n', Tf);

for k = 1:length(h_values)
    h = h_values(k);
    Nx = Nx_values(k);
    
    fprintf('Computing solution for h = 1/%d (Nx = %d)...\n', round(1/h), Nx);
    
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
    
    % Compute L2 error at each time step
    max_error = 0;
    for n = 1:length(T)
        U_exact_n = u_exact(X_grid, T(n))';
        U_approx_n = U_hist(:, n);
        
        % Compute L2 norm of error
        dx = X_grid(2) - X_grid(1);
        error_n = sqrt(dx * sum((U_approx_n - U_exact_n).^2));
        
        % Track maximum error
        max_error = max(max_error, error_n);
    end
    
    E_h(k) = max_error;
    fprintf('  E_h = %.6e\n', E_h(k));
end

% Compute convergence rates
fprintf('\n========== Convergence Rates ==========\n');
fprintf('h\t\tE_h\t\t\tln(E_2h/E_h)/ln(2)\n');
fprintf('%.4f\t\t%.6e\tN/A\n', h_values(1), E_h(1));

for k = 2:length(h_values)
    convergence_rates(k) = log(E_h(k-1) / E_h(k)) / log(2);
    fprintf('%.4f\t\t%.6e\t%.6f\n', h_values(k), E_h(k), convergence_rates(k));
end

% Create Table 1: Convergence Table
fprintf('\n========== Table 1: h-rate Convergence Table ==========\n');
fprintf('h\t\tE_h = max|u(t_n) - u_h^n|_L2\tln(E_2h/E_h)/ln(2)\n');
fprintf('%-10s\t%-30s\t%-20s\n', '1/2', sprintf('%.6e', E_h(1)), 'N/A');
for k = 2:length(h_values)
    fprintf('%-10s\t%-30s\t%-20.6f\n', sprintf('1/%d', round(1/h_values(k))), sprintf('%.6e', E_h(k)), convergence_rates(k));
end

% Create log-log plot
figure;
log_h = -log(h_values);
log_E = -log(E_h);

plot(log_h, log_E, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('-log(h)', 'FontSize', 12);
ylabel('-log(E_h)', 'FontSize', 12);
title('Spatial Convergence: Log-Log Plot', 'FontSize', 14);
grid on;
hold on;

% Fit a line to estimate convergence order
p = polyfit(log_h, log_E, 1);
log_h_fit = linspace(min(log_h), max(log_h), 100);
log_E_fit = polyval(p, log_h_fit);
plot(log_h_fit, log_E_fit, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Linear fit (slope=%.2f)', p(1)));
legend('Numerical results', 'Linear fit', 'FontSize', 11);

hold off;
saveas(gcf, fullfile('Project', 'convergence_loglog_plot.png'));
fprintf('\nLog-log plot saved to Project/convergence_loglog_plot.png\n');
close;

disp('Task 3 complete!!!!!!');
disp(sprintf('Estimated convergence order: %.2f', p(1)));