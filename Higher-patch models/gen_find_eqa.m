function equilibria = gen_find_eqa(time_limit,num_patches, alpha, beta, K, d)
    
    % Set up the appropriate adjacency matrix

    %   A = [0,1,1,0,0,1,1;
    %       1,0,1,1,0,0,1;
    %       1,1,0,1,1,0,0;
    %       0,1,1,0,1,1,0;
    %       0,0,1,1,0,1,1;
    %       1,0,0,1,1,0,1;
    %       1,1,0,0,1,1,0];
   
    %   A = [0,1,0,0,0,0,1;
    %        1,0,1,0,0,0,0;
    %        0,1,0,1,0,0,0;
    %        0,0,1,0,1,0,0;
    %        0,0,0,1,0,1,0;
    %        0,0,0,0,1,0,1;
    %        1,0,0,0,0,1,0];

    %   A = [0,1,1,1,1,1,1;
    %      1,0,1,1,1,1,1;
    %      1,1,0,1,1,1,1;
    %      1,1,1,0,1,1,1;
    %      1,1,1,1,0,1,1;
    %      1,1,1,1,1,0,1;
    %      1,1,1,1,1,1,0];   
    

    A = find_ad_matrix(num_patches,num_connec);

    % Find the Laplacian matrix
    L = A - diag(sum(A, 1));
    
    % Define the ODE system
    f = @(u, v) u.*((1-1/K.*u)-v./(1+u));
    g = @(u, v) v.*(alpha./(1+u).*u-beta);
    
    ui = 1:num_patches;
    vi = num_patches + 1:2*num_patches;

    odeSystem = @(U) [f(U(ui), U(vi)); g(U(ui), U(vi)) + d*L*U(vi)];
    
    % Set options for fsolve
    options = optimset('Display', 'off', 'Algorithm', 'trust-region-dogleg');

    % Initialize matrix to store equilibrium solutions
    equilibria = [];

    % Set a time limit for searching equilibrium solutions
    timeLimit = time_limit;

    % Use fsolve in a loop with different initial guesses
    tic; % Start the timer

    while true
        initialGuess = (24/24*K) * rand(2 * num_patches, 1);  % Random initial guess
        equilibrium = fsolve(odeSystem, initialGuess, options);
        
        % Check if the equilibrium satisfies the ODE system
        odeCheck = odeSystem(equilibrium);
        if all(abs(odeCheck) < 1e-6)

            % Remove rows with negative values
            if all(equilibrium >= -1e-6)

                % Check if the equilibrium is not already in the list
                isUnique = true;
                for j = 1:size(equilibria, 1)
                    if max(abs(transpose(equilibrium) - equilibria(j, :))) < 1e-2
                        isUnique = false;
                        break;
                    end
                end
                
                % Check if unique
                if isUnique
                    equilibria = [equilibria; transpose(equilibrium)];
                end

            end
        end

        % Check the elapsed time
        if toc > timeLimit
            break;
        end

    end

    % Remove rows with infinite values    
    rowsWithInf = any(~isfinite(equilibria), 2);
    equilibria(rowsWithInf, :) = [];

end


