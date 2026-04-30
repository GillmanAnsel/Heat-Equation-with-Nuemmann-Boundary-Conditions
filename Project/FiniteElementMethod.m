%% STEP 1: Main FEM Solver Function (Heat Equation with Neumann BCs)
function [X_grid, T, U_final, U_history] = FiniteElementMethod(varargin)
 
    % --- Parse input parameters ---
    p = inputParser;
    addParameter(p, 'L', 1);           % Domain length
    addParameter(p, 'Nx', 20);         % Number of spatial nodes
    addParameter(p, 'Tf', 5);          % Final time
    addParameter(p, 'Nt', 20);         % Number of time steps
    addParameter(p, 'a_func', @(x) 1); % Diffusion coefficient
    addParameter(p, 'f_func', @(x,t) 0); % Source term
    addParameter(p, 'g_left', @(t) 0);   % Left Neumann BC
    addParameter(p, 'g_right', @(t) 0);  % Right Neumann BC
    addParameter(p, 'u0', @(x) 0);       % Initial condition
    addParameter(p, 'save_plots', false);% Save plots flag
    addParameter(p, 'plot_title', 'Heat Equation Solution');    % Plot title
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
    plot_title = params.plot_title;

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
        B(1) = -a_func(0) * g_left(t_n);      % Left boundary (outward normal is -x)
        B(Nx) = a_func(L)  * g_right(t_n);   % Right boundary (outward normal is +x)
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
        plot_fem_results(X_grid, T, U_history, U_final, plot_title);
    end
end

%% STEP 2: Plotting Functions
function plot_fem_results(X_grid, T, U_history, U_final, plot_title)
    % Plot solution evolution as a 3D surface
    [T_mesh, X_mesh] = meshgrid(T, X_grid);
    figure;
    surf(T_mesh, X_mesh, U_history, 'EdgeColor', 'none');
    colorbar;
    xlabel('Time t');
    ylabel('Position x');
    zlabel('u(x, t)');
    title(plot_title, 'Interpreter', 'latex');
    view(45, 30);  % Set viewing angle
    saveas(gcf, 'solution_3d_surface.png');
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

%% Plot Hat Functions
function plot_hat_functions(X_grid)
    % Plotting some hat functions for visualization
    figure;
    hold on;
    for k = 1:10:length(X_grid)
        phi_k = hat_function(X_grid, k);
        plot(X_grid, phi_k, 'LineWidth', 1.5);
    end
    xlabel('x');
    ylabel('Hat Functions');
    title('Plotting (Some) Hat Functions');
    grid on;
    saveas(gcf, 'hat_functions.png');
    close;
end

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
