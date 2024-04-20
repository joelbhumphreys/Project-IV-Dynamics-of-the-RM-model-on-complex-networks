for i = [1,2,3,4,5,6,7,8,9,10]
    rng(i); % set the random seed
    
    num_patches = 7; % set the number of patches
    
    A = find_ad_matrix(num_patches); % adjacency matrix
    
    beta = 1;
    K = 2.51;
    alpha = 9.96; 
    d = 0.5;
    
    L = A - diag(sum(A, 1)); % Laplacian matrix
    
    % Define ODE system
    f = @(u, v) u.*((1-1/K.*u)-v./(1+u));
    g = @(u, v) v.*(alpha./(1+u).*u-beta);
    
    ui = 1:num_patches;
    vi = num_patches + 1:2*num_patches;
    
    odeSystem = @(t, U) [f(U(ui), U(vi)); g(U(ui), U(vi)) + d*L*U(vi)];
    
    initialConditions = 1 + abs(1e1 * randn(2*size(A,1),1))';
    
    I = eye(num_patches);
    JacSparse = sparse([I, I; I, L]);
    options = odeset('reltol', 1e-4, 'abstol', 1e-5, 'JPattern', JacSparse); % settings for ode45
    % options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11, 'JPattern', JacSparse, 'MaxStep',0.1);
    
    [t, U] = ode45(odeSystem, linspace(0, 1e4, 1e3), initialConditions, options);
    [t, U] = ode45(odeSystem, linspace(0, 7500000, 1000), U(end,:), options);
    
    % Extract the x-component of the system
    x_component = U(:, 1);
    
    % Apply the z1test function
    kmedian = z1test(x_component);
    
    % Display the result
    disp(['Result of the 0-1 test for chaos: ', num2str(kmedian)]);
    
end



