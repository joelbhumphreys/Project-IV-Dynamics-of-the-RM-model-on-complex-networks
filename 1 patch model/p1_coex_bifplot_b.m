% Set parameter values
beta = 11;
alphaList = 0.01:0.05:20.01;
KList = 0.01:0.05:20.01;

parameterMatrix = cell(length(alphaList), length(KList));

for i = 1:length(alphaList)
    for j = 1:length(KList)
        parameterMatrix{i, j} = [alphaList(i), KList(j)];
    end
end

% Construct matrix to hold the colour outputs
colourMatrix = cell(length(alphaList), length(KList));

% Iterate over the parameter space 
alpha_index = 0;

for alpha = 0.01:0.05:20.01
    K_index = 0;
    alpha_index = alpha_index + 1;
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
                colourMatrix{alpha_index, K_index} = [1, 0, 0];
            elseif sign(eigenvalues(1)) > 0 && sign(eigenvalues(2)) > 0
                colourMatrix{alpha_index, K_index} = [0, 1, 0];
            elseif sign(eigenvalues(1)) ~= sign(eigenvalues(2))
                colourMatrix{alpha_index, K_index} = [0, 0, 1];
    
            end
         
        else

            if real(eigenvalues(1)) < 0
                colourMatrix{alpha_index, K_index} = [255, 255, 0];
            
            elseif real(eigenvalues(1)) == 0
                colourMatrix{alpha_index, K_index} = [255, 165, 0];

            elseif real(eigenvalues(1)) > 0
                colourMatrix{alpha_index, K_index} = [128, 0, 128];
            
            end
        end

    end
end

% Extract coordinates from the parameter matrix
coordinates = cell2mat(parameterMatrix(:));

% Extract colours from the colour matrix
colours = cell2mat(colourMatrix(:));

% Create a scatter plot of the colours
scatter(coordinates(:, 1), coordinates(:, 2), 50, colours, 'filled');
xlabel('\alpha', 'FontSize', 25);
ylabel('K', 'FontSize', 19);

xlim([0, 20]); 
ylim([0, 20]); 



