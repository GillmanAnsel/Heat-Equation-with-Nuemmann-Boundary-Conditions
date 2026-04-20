%% ----------APDE Project----------
%%-----Task 1-----

%% Step 1 - using implicit Euler method for time discretization
% Grid setup 
Tf = 1; % final time
nt = 100; % number of time steps
T = linspace(0, Tf, nt); % time vector


% Step 2 - derive the weak form 
% Step 3 - implement the finite element method
L = 1;
Nx = 100; % number of spatial steps
X_grid= linspace(0, L, Nx); % spatial vector
phi = zeros(Nx, 1); % piecewise linear basis functions
function phi = hat_function(x, X_grid, k)
    % x: vector of input points
    % x_coords: the full mesh/grid
    % k: the index of the node
    h = X_grid(2) - X_grid(1); % spatial step size
    % Initialize phi to zeros

    % Making half hats for end nodes
    phi = zeros(size(x));
    if k == 1
        xk = X_grid(k);
        xk_next = X_grid(k+1);
        idx1 = (x >= xk & x <= xk_next);
        phi(idx1) = (xk_next - x(idx1)) / h;
        return;
    elseif k == length(X_grid)
        xk = X_grid(k);
        xk_prev = X_grid(k-1);
        idx2 = (x >= xk_prev & x <= xk);
        phi(idx2) = (x(idx2) - xk_prev) / h;
        return;
    end

    % full hats for interior nodes
    % Step 1: 
    % Define points and intervals
    xk_prev = X_grid(k-1);
    xk      = X_grid(k);
    xk_next = X_grid(k+1);
    % half hats
    % Left slope (x_{k-1} to x_k)
    idx1 = (x >= xk_prev & x < xk);
    phi(idx1) = (x(idx1) - xk_prev) / h ;
    
    % Step 2: Right slope (x_k to x_{k+1})
    idx2 = (x >= xk & x <= xk_next);
    phi(idx2) = (xk_next - x(idx2)) / h;
    
    % Everything else remains 0
end

% plotting some hat functions
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;
for k = [1, 25, 50, 75, 100]
    phi_k = hat_function(X_grid, X_grid, k);
    plot(X_grid, phi_k, 'b');
end
saveas(gcf, 'Project\hat_functions.png');

% Step 4- formulate the fully discrete scheme and the system of equations
disp('Task 1 complete!!!!!');