function eigenvalues = find_eig(u10,u20,v10,v20,alpha,beta,K,d)

    % Evaluate Jacobian at equilibrium point
    M = [1 - 2/K*u10 - v10/(1+u10)^2, 0, -u10/(1+u10), 0;
        0, 1 - 2/K*u20 - v20/(1+u20)^2, 0, -u20/(1+u20);
        (alpha*v10)/(1+u10)^2, 0, (alpha*u10)/(1+u10) - beta - d, d;
        0, (alpha*v20)/(1+u20)^2, d, (alpha*u20)/(1+u20) - beta - d];

    % Find eigenvalues
    eigenvalues = eig(M);
    
end
