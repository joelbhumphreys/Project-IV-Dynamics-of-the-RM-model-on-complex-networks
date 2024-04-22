rng(1);

% Plot Poincare section

beta = 10;    
K = 10;
alpha = 16;
d = 1; 

currentYLim = ylim;
currentZLim = zlim;

% Generate a grid of points in 3D
[y, z] = meshgrid(currentYLim(1):0.1:currentYLim(2), currentZLim(1):0.1:currentZLim(2));

x = alpha / (alpha - beta); % define the constant coordinate of the plane

% Plot the plane 
figure;
surf(x, y, z, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlay the plane with a plot of a trajectory in a 3D subspace 

num_patches = 7; % number of patches

% Set up appropriate adjacency matrix

% A = find_ad_matrix(num_patches);
    
A = [0,1,1,0,0,1,1;
     1,0,1,1,0,0,1;
     1,1,0,1,1,0,0;
     0,1,1,0,1,1,0;
     0,0,1,1,0,1,1;
     1,0,0,1,1,0,1;
     1,1,0,0,1,1,0];  

L = A - diag(sum(A, 1)); % define Laplacian matrix

% Define ODE system
f = @(u, v) u.*((1-1/K.*u)-v./(1+u));
g = @(u, v) v.*(alpha./(1+u).*u-beta);

ui = 1:num_patches;
vi = num_patches + 1:2*num_patches;

odeSystem = @(t, U) [f(U(ui), U(vi)); g(U(ui), U(vi)) + d*L*U(vi)];

initialConditions = 1 + abs(1e1 * randn(2*size(A,1),1))';

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5); % settings for ode45

[t, U] = ode45(odeSystem, linspace(0, 1e4, 1e3), initialConditions, options); % disgarded burn-in period
[t, U] = ode45(odeSystem, linspace(0, 2*1e4, 2*1e4), U(end,:), options);

% Plot the 3D trajectory
plot3(U(1*end/2:end, 1), U(1*end/2:end, 7), U(end/2:end, 12), 'LineWidth', 1);

xlabel('$u_1$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$u_7$', 'Interpreter', 'latex', 'FontSize', 16);
zlabel('$v_5$', 'Interpreter', 'latex', 'FontSize', 16);

hold off;



