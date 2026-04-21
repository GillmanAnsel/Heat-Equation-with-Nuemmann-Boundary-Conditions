%% ----------APDE Project----------
%%-----Task 1-----

%% STEP 1: Time Discretization using Implicit Euler Method
% =========================================================
% Set up uniform time partition: t_n = n*tau, where tau = Tf/N
Tf = 1; % final time
nt = 20; % number of time steps
T = linspace(0, Tf, nt); % time vector
dt = T(2) - T(1); % time step size

%% STEP 2: Weak Formulation
% =========================================================
% Starting from the heat equation:
%   du/dt - d/dx(a(x) du/dx) = f(x,t)    in (0,L)
%   du/dx(0,t) = a,  du/dx(L,t) = b      Neumann BCs
%   u(x,0) = u0(x)                       Initial condition
%
% Multiply by test function v, integrate by parts to get weak form:
%   (1/dt)(u^n, v) + (a du^n/dx, dv/dx) = (1/dt)(u^{n-1}, v) + (f^n, v) + boundary_terms
%
% This leads to the system:  [(1/dt)M + K] U^n = (1/dt)M U^{n-1} + F^n + B^n

%% STEP 3: Finite Element Method with Piecewise Linear Basis Functions
% =========================================================
% Set up spatial domain
L = 1; % domain length
Nx = 20; % number of spatial nodes
X_grid = linspace(0, L, Nx); % uniform spatial mesh
h = X_grid(2) - X_grid(1); % mesh spacing

% Define problem parameters
a_func = @(x) 1; % diffusion coefficient a(x)
f_func = @(x,t) 0; % source term f(x,t)
g_left = @(t) 0; % Neumann BC at x=0: du/dx(0,t) = a
g_right = @(t) 0; % Neumann BC at x=L: du/dx(L,t) = b
u0 = @(x) 0; % initial condition u(x,0) = u0(x)

% plotting some hat functions
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;
for k = 1:5:Nx
    phi_k = hat_function(X_grid, k);
    plot(X_grid, phi_k, 'b');
end
xlabel('x');
ylabel('\phi_k(x)');
title('Piecewise Linear Hat Functions (Basis Functions)');
saveas(gcf, fullfile('Project', 'hat_functions.png'));
close;

%% STEP 4: Fully Discrete Scheme and System of Linear Equations
% =========================================================
% Assemble matrices:
%   K_ij = integral a(x)*phi_i'(x)*phi_j'(x) dx  (stiffness matrix)
%   M_ij = integral phi_i(x)*phi_j(x) dx  (mass matrix)
K = StiffnessAssembler1D(X_grid, a_func, [0,0]);
M = MassAssembler1D(X_grid);

% System matrix: A = (1/dt)*M + K
A_system = (1/dt)*M + K;

% Initial condition (expand in basis functions)
U_prev = u0(X_grid)';
U_current = U_prev;

% Solving the System [(1/dt)M + K] U^n = (1/dt)M U^{n-1} + F^n + B^n
for n = 2:length(T)
    t_n = T(n);
    
    % Assemble load vector: F_i = integral f(x,t_n)*phi_i(x) dx
    F = LoadAssembler1D(X_grid, f_func, t_n);
    
    % Boundary contributions (Neumann conditions on RHS)
    B = zeros(Nx, 1);
    B(1) = g_left(t_n);
    B(Nx) = g_right(t_n);
    
    % Right-hand side: b = (1/dt)*M*U^{n-1} + F + B
    b = (1/dt)*M*U_prev + F + B;
    
    % Solve the linear system
    U_current = A_system \ b;
    
    % Update for next time step
    U_prev = U_current;
end

% Visualize final solution
figure;
plot(X_grid, U_current, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('u(x, T)');
title(sprintf('FEM Solution of Heat Equation at t = %.2f', T(end)));
grid on;
saveas(gcf, fullfile('Project', 'solution.png'));
close;

disp('Task 1 complete!!!!!');

%% ========== FUNCTIONS ==========

function phi = hat_function(X_grid, k)
    % Piecewise linear hat (tent) basis function
    % phi_k is 1 at node k, 0 elsewhere, linear between nodes
    h = X_grid(2) - X_grid(1);
    phi = zeros(size(X_grid));
    
    if k == 1
        % Left boundary: linear from 1 at x_1 to 0 at x_2
        xk = X_grid(k);
        xk_next = X_grid(k+1);
        idx1 = (X_grid >= xk & X_grid <= xk_next);
        phi(idx1) = (xk_next - X_grid(idx1)) / h;
    elseif k == length(X_grid)
        % Right boundary: linear from 0 at x_{n-1} to 1 at x_n
        xk = X_grid(k);
        xk_prev = X_grid(k-1);
        idx2 = (X_grid >= xk_prev & X_grid <= xk);
        phi(idx2) = (X_grid(idx2) - xk_prev) / h;
    else
        % Interior: tent with peak at x_k
        xk_prev = X_grid(k-1);
        xk = X_grid(k);
        xk_next = X_grid(k+1);
        
        % Left slope
        idx1 = (X_grid >= xk_prev & X_grid < xk);
        phi(idx1) = (X_grid(idx1) - xk_prev) / h;
        
        % Right slope
        idx2 = (X_grid >= xk & X_grid <= xk_next);
        phi(idx2) = (xk_next - X_grid(idx2)) / h;
    end
end

function A = StiffnessAssembler1D(x, a, kappa)
    % Assemble 1D stiffness matrix K_ij = integral a(x)*phi_i'(x)*phi_j'(x) dx
    % using piecewise linear hat functions on uniform mesh
    n = length(x)-1;
    A = zeros(n+1, n+1);
    for i = 1:n
        h = x(i+1) - x(i);
        xmid = (x(i+1) + x(i))/2; % interval mid-point
        amid = a(xmid); % value of a(x) at mid-point
        A(i,i) = A(i,i) + amid/h;
        A(i,i+1) = A(i,i+1) - amid/h;
        A(i+1,i) = A(i+1,i) - amid/h;
        A(i+1,i+1) = A(i+1,i+1) + amid/h;
    end
    % Add boundary terms (Robin/Dirichlet penalties if kappa > 0)
    A(1,1) = A(1,1) + kappa(1);
    A(n+1,n+1) = A(n+1,n+1) + kappa(2);
end

function M = MassAssembler1D(x)
    % Assemble 1D mass matrix M_ij = integral phi_i(x)*phi_j(x) dx
    % For uniform mesh with piecewise linear hat functions:
    % M_ii = 2*h/3 (diagonal), except endpoints M_00 = M_NN = h/3
    % M_{i,i+1} = M_{i+1,i} = h/6 (off-diagonal)
    n = length(x)-1;
    M = zeros(n+1, n+1);
    h = x(2) - x(1); % uniform mesh spacing
    
    % Diagonal entries
    M(1,1) = h/3;
    for i = 2:n
        M(i,i) = 2*h/3;
    end
    M(n+1,n+1) = h/3;
    
    % Off-diagonal entries (neighbors)
    for i = 1:n
        M(i,i+1) = h/6;
        M(i+1,i) = h/6;
    end
end

function F = LoadAssembler1D(x, f_func, t)
    % Assemble 1D load vector F_i = integral f(x,t)*phi_i(x) dx
    % Using midpoint quadrature on each element
    n = length(x)-1;
    F = zeros(n+1, 1);
    h = x(2) - x(1); % uniform mesh spacing
    
    % Contribution from each interval [x_i, x_{i+1}]
    for i = 1:n
        xmid = (x(i) + x(i+1)) / 2; % midpoint of interval
        f_mid = f_func(xmid, t); % evaluate f at midpoint
        
        % Each node receives h/2 of the contribution per adjacent interval
        F(i) = F(i) + f_mid * h / 2;
        F(i+1) = F(i+1) + f_mid * h / 2;
    end
end


