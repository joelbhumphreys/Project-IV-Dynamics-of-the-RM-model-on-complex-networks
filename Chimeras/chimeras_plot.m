rng(1); % Set random seed

num_patches = 100; % Set number of patches

% Set parameter values
d = 1e-3;
K = 4;
beta = 1;

% Iterate over alpha
for i = [1]
    alpha = i; % set value of alpha
        
    A = find_ad_matrix2(num_patches); % set up appropriate adjacency matrix
    L = A - diag(sum(A, 1)); % set Laplacian matrix
    
    x0 = 1 + abs(1e1 * randn(2 * size(A, 1), 1)); % initial conditions

    % Define ODE system
    odeSystem = @(t, u) [u(1:num_patches).*((1-1/K.*u(1:num_patches))-u(num_patches+1:end)./(1+u(1:num_patches)));
                 u(num_patches+1:end).*(alpha./(1+u(1:num_patches)).*u(1:num_patches)-beta) + d*L*u(num_patches+1:end)];
        
    % Run ODE system
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);  % options for ode45
    
    [t, xx] = ode45(odeSystem, linspace(0, 1e4, 1e5), x0, options); % discarded burn-in period
    [t, xx] = ode45(odeSystem, linspace(0, 1e3, 1e4), xx(end,:), options);  
    
    % Generate heat map 
    figure;
    imagesc(xx);
    colorbar;
    caxis([min(xx(:)), max(xx(:))/10]); % Adjust colorbar range
    
    xlabel('Node Index');
    ylabel('Time');
    

end    
    

   

