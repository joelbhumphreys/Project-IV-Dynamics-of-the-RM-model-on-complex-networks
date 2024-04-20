% Iterate over values of K
for K = [6, 11, 16]
    % Set the ranges for the parameter space        
    alphaList = 22.01:0.01:47.01;
    betaList = 15.01:0.01:35.01;
    
    matrixWithCells = cell(length(alphaList), length(betaList));
    
    for i = 1:length(alphaList)
        for j = 1:length(betaList)
            matrixWithCells{i, j} = [alphaList(i), betaList(j)];
        end
    end
    
    % Construct matrix to hold the colour outputs            
    colourMatrix = cell(length(alphaList), length(betaList));
    
    % Iterate over the parameter space             
    alpha_index = 0;
    
    for alpha = 22.01:0.01:47.01
        beta_index = 0;
        alpha_index = alpha_index + 1;
        for beta = 15.01:0.01:35.01
            beta_index = beta_index + 1;
    
            % Find eigenvalues of the Jacobian at the predator-free equilibrium and determine its stability                                           
            A = [-1, -K/(1+K);
                0, alpha*K/(1+K) - beta];
            eigenvaluesA = eig(A);        
            
    
            if real(eigenvaluesA(1)) < 0 && real(eigenvaluesA(2)) < 0
                A_stab = true;
            else
                A_stab = false;
            end
            
            % Find eigenvalues of the Jacobian at the coexistence equilibrium and determine its stability                                                       
            B = [(beta*(1+1/K))/alpha - (2*beta)/(K*(alpha-beta)), -beta/alpha;
                alpha - beta - beta/K, 0];
            eigenvaluesB= eig(B);        
    
            if real(eigenvaluesB(1)) < 0 && real(eigenvaluesB(2)) < 0
                B_stab = true;
            else
                B_stab = false;
            end

            % Determining the overall dynamics of the system                        
            if alpha <= beta && B_stab == true
                colourMatrix{alpha_index, beta_index} = [1,1,1];  % check for feasibility
            
            else
                % Produce colour plots
                if A_stab == true && B_stab == true
                    colourMatrix{alpha_index, beta_index} = [0, 1, 0]; % both are stable 
                
                elseif A_stab == true && B_stab == false
                    colourMatrix{alpha_index, beta_index} = [255, 255, 0]; % only predator-free is stable
        
                elseif A_stab == false && B_stab == true
                    colourMatrix{alpha_index, beta_index} = [0, 0, 1]; % only coexistence is stable
        
                elseif A_stab == false && B_stab == false
                    colourMatrix{alpha_index, beta_index} = [1, 0, 0]; % neither eqm point is stable
                    
                end
            end

        end
    end
    
    % Extract coordinates from the parameter matrix
    coordinates = cell2mat(matrixWithCells(:));
    
    % Extract colours from the colour matrix
    colors = cell2mat(colourMatrix(:));
    
    % Create a scatter plot of the colours
    figure;
    scatter(coordinates(:, 1), coordinates(:, 2), 50, colors, 'filled');

    title_str = sprintf('$K = %g$', K);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 25);    
    xlabel('\alpha', 'FontSize', 25);
    ylabel('\beta', 'FontSize', 25);
    
    % Adjust plot limits to avoid dots going over the axis
    axis([min(alphaList) max(alphaList) min(betaList) max(betaList)]);
    
end


