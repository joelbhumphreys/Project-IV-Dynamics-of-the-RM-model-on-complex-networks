% Set parameter values
alpha = 11;
betaList = 0.01:0.05:20.01;
KList = 0.01:0.05:20.01;

matrixWithCells = cell(length(betaList), length(KList));

for i = 1:length(betaList)
    for j = 1:length(KList)
        matrixWithCells{i, j} = [betaList(i), KList(j)];
    end
end

% Construct matrix to hold the colour outputs
colourMatrix = cell(length(betaList), length(KList));

% Iterate over the parameter space 
beta_index = 0;

for beta = 0.01:0.05:20.01
    K_index = 0;
    beta_index = beta_index + 1;
    for K = 0.01:0.05:20.01
        K_index = K_index + 1;

        % Find eigenvalues of the Jacobian at the coexistence equilibrium         
        A = [(beta*(1+1/K))/alpha - (2*beta)/(K*(alpha-beta)), -beta/alpha;
            alpha - beta - beta/K, 0];
        
        eigenvalues = eig(A);

        % Identify the type of dynamics         
        isComplexRealNumber = ~isreal(eigenvalues(1));
        if isComplexRealNumber == 0 

            if sign(eigenvalues(1)) < 0 && sign(eigenvalues(2)) < 0
                colourMatrix{beta_index, K_index} = [1, 0, 0];
            elseif sign(eigenvalues(1)) > 0 && sign(eigenvalues(2)) > 0
                colourMatrix{beta_index, K_index} = [0, 1, 0];
            elseif sign(eigenvalues(1)) ~= sign(eigenvalues(2))
                colourMatrix{beta_index, K_index} = [0, 0, 1];
            end
        
        else

            if real(eigenvalues(1)) < 0
                colourMatrix{beta_index, K_index} = [255, 255, 0];
            elseif real(eigenvalues(1)) == 0
                colourMatrix{beta_index, K_index} = [255, 165, 0];
            elseif real(eigenvalues(1)) > 0
                colourMatrix{beta_index, K_index} = [128, 0, 128];
            end
            
        end

    end
end

% Extract coordinates from the parameter matrix
coordinates = cell2mat(matrixWithCells(:));

% Extract colours from the colour matrix
colours = cell2mat(colourMatrix(:));

% Create a scatter plot of the colours
scatter(coordinates(:, 1), coordinates(:, 2), 50, colours, 'filled');
xlabel('\beta', 'FontSize', 25);
ylabel('K', 'FontSize', 19);

xlim([0.01, 20]); 
ylim([0.01, 20]); 




