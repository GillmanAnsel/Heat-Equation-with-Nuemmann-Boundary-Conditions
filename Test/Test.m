% Simple 1D Heat Diffusion Test (No Toolboxes)
% Grid setup
L = 1;              % Length of rod
Nx = 50;            % Spatial steps
dx = L/(Nx-1);      
x = linspace(0, L, Nx);

% Initial Condition: A "spike" of heat in the center
u = zeros(1, Nx);
u(round(Nx/2)) = 100; 

% Simulation parameters
alpha = 0.01;       % Thermal diffusivity
dt = 0.01;          % Time step
steps = 500;

% Finite Difference Loop
for t = 1:steps
    u_new = u;
    for i = 2:Nx-1
        % Standard second-order central difference
        u_new(i) = u(i) + alpha * dt/dx^2 * (u(i+1) - 2*u(i) + u(i-1));
    end
    u = u_new;
end

% Plot final result
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
plot(x, u, 'r');
title(['Heat Diffusion at step ', num2str(steps)]);
xlabel('Position (x)');
ylabel('Temperature (u)');
grid on;
drawnow;

% Save single output graph
frame = getframe(gcf);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
saveas(gcf, 'Test\test_graph.png');
disp('Test complete!!!!!');