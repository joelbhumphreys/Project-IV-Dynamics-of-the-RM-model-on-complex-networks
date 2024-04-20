function equilibria = find_eqa(alpha, beta, K,d)
    % Define the ODE system
    odeSystem = @(y) [
        y(1) * ((1 - y(1)/K) - y(3)/(1 + y(1)));
        y(2) * ((1 - y(2)/K) - y(4)/(1 + y(2)));
        y(3) * ((alpha * y(1))/(1 + y(1)) - beta) + d * (y(4) - y(3));
        y(4) * ((alpha * y(2))/(1 + y(2)) - beta) + d * (y(3) - y(4))
    ];

    % Set options for fsolve
    options = optimset('Display', 'off', 'Algorithm', 'trust-region-dogleg');

    % Initialize matrix to store equilibrium solutions
    equilibria = [];

    % Set a time limit for searching equilibrium solutions 
    timeLimit = 0.1;

    % Use fsolve in a loop with different initial guesses
    tic; % Start the timer

    while true
        initialGuess = K * rand(4, 1);  % Random initial guess
        equilibrium = fsolve(odeSystem, initialGuess, options);

        % Remove rows with negative values
        isNegative = true;
        if any(equilibrium > -1e-3, 2)
            isNegative = false;
            equilibria = abs(equilibria);
        end

        % Check if the equilibrium is already in the list
        isUnique = true;
        for j = 1:size(equilibria, 1)
            if max(abs(transpose(equilibrium) - equilibria(j, :))) < 1e-2
                isUnique = false;
                break;
            end
        end

        if isUnique && isNegative == false
            equilibria = [equilibria; transpose(equilibrium)];
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


