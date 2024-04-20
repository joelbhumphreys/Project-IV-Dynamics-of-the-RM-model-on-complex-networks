rng(1);

num_patches = 7; % set number of patches

% Set up appropriate adjacency matrix

A = find_ad_matrix(num_patches); 

% A = [0,1,1,0,0,1,1;
%      1,0,1,1,0,0,1;
%      1,1,0,1,1,0,0;
%      0,1,1,0,1,1,0;
%      0,0,1,1,0,1,1;
%      1,0,0,1,1,0,1;
%      1,1,0,0,1,1,0];

% Set parameter values
beta = 10;    
K = 10;
alpha = 16;
d = 1;  

L = A - diag(sum(A, 1)); % find Laplacian matrix

% Define ODE system
f = @(u, v) u.*((1-1/K.*u)-v./(1+u));
g = @(u, v) v.*(alpha./(1+u).*u-beta);

ui = 1:num_patches;
vi = num_patches + 1:2*num_patches;

odeSystem = @(t, U) [f(U(ui), U(vi)); g(U(ui), U(vi)) + d*L*U(vi)];

initialconditions = 1 + abs(1e1 * randn(2*size(A,1),1))';

% Solve the ODE system
I = eye(num_patches);

JacSparse = sparse([I,I;I,L]);

options = odeset('reltol', 1e-11, 'abstol', 1e-11, 'JPattern', JacSparse);

[t, U] = ode45(odeSystem, linspace(0, 1e4, 1e3), initialconditions, options);
[t, U] = ode45(odeSystem, linspace(0, 2*1e4, 2*1e4), U(end,:), options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time plot

figure;

num_colours = num_patches * 2;
colours = jet(num_colours);

% Plot u values with different colours
for i = 1:num_patches
    plot(t, U(:, i), 'LineWidth', 2, 'Color', colours(i, :));
    hold on;
end

% Plot v values with different colours
for i = num_patches + 1:2*num_patches
    plot(t, U(:, i), 'LineWidth', 2, 'Color', colours(i, :));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D plot
figure;
plot3(U(:, 1), U(:, 7), U(:, 12))
xlabel('$u_1$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$v_1$', 'Interpreter', 'latex', 'FontSize', 16)
zlabel('$v_5$', 'Interpreter', 'latex', 'FontSize', 16)

xlim([0, 13]);
ylim([0, 17]);
zlim([0, 16]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous?
figure
plot(U(end/2:5:end,1),'*')

% xlim([0, 10000]);
% ylim([0, 12]);


