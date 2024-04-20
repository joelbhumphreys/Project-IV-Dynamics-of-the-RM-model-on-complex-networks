rng(1); % set random seed

num_patches = 7; % set the number of patches

% Find the adjacency matrix

A = find_ad_matrix(num_patches,num_connec);

% Find the Laplacian matrix
L = A - diag(sum(A, 1));

% Set the parameter value
alpha = 6;
beta = 4;
K = 7;
d = 1;
   
% Define the ODE system
f = @(u, v) u.*((1-1/K.*u)-v./(1+u));
g = @(u, v) v.*(alpha./(1+u).*u-beta);

ui = 1:num_patches;
vi = num_patches + 1:2*num_patches;

odeSystem = @(t, U) [f(U(ui), U(vi)); g(U(ui), U(vi)) + d*L*U(vi)];

initialConditions = 1+abs(1e1*randn(2*size(A,1),1))'; % set the initial conditions

% Solve the ODE system
I = eye(num_patches);
JacSparse = sparse([I, I; I, L]);

options = odeset('reltol', 1e-4, 'abstol', 1e-5, 'JPattern', JacSparse);

[t, U] = ode45(odeSystem, linspace(0, 1e3, 1e4), initialConditions, options); % discarded burn-in period
[t, U] = ode45(odeSystem, linspace(0, 1e3, 1e4), U(end,:), options);

figure;

num_colors = num_patches * 2;
colors = jet(num_colors);

% Plot u values with different colors
for i = 1:num_patches
    plot(t, U(:, i), 'LineWidth', 2, 'Color', colors(i, :));
    hold on;
end

% Plot v values with different colors
for i = num_patches + 1:2*num_patches
    plot(t, U(:, i), 'LineWidth', 2, 'Color', colors(i, :));
    hold on;
end

hold off;

% Generate legend entries 
legend_entries = cell(1, 2*num_patches);
for i = 1:num_patches
    legend_entries{i} = ['u_', num2str(i)];
end
for i = num_patches + 1:2*num_patches
    legend_entries{i} = ['v_', num2str(i - num_patches)];
end

legend(legend_entries, 'FontSize', 14);

xlabel('Time');
ylabel('Population');

grid on;



