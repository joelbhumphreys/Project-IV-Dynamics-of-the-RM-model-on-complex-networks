% Iterate over d
for d = [1, 5, 20]
    % Set the ranges for the parameter space        
    betaList = 5.01:0.13:15.01;
    KList = 0.01:0.13:12.01;
    alpha = 16;   

    parameterMatrix = cell(length(betaList), length(KList));
    
    for i = 1:length(betaList)
        for j = 1:length(KList)
            parameterMatrix{i, j} = [betaList(i), KList(j)];
        end
    end
    
    % Construct matrix to hold the colour outputs            
    colourMatrix = cell(length(betaList), length(KList));
    
    % Iterate over the parameter space                 
    beta_index = 0;
    
    for beta = 5.01:0.13:15.01
        K_index = 0;
        beta_index = beta_index + 1;

        for K = 0.01:0.13:12.01
            K_index = K_index + 1;  
                
            % Analytical forms of most of the equilibrium points
            eqa = [K, K, 0, 0;
            0, (d + beta - (d^2)/(d+beta))/(alpha - d - beta + d^2/(d+beta)), d/(d+beta)*(1 - 1/K*(d + beta - d^2/(d+beta))/(alpha - d - beta + (d^2)/(d+beta))) * ( alpha/(alpha - d - beta + d^2/(d+beta))), ( 1-1/K * (d + beta - d^2/(d+beta))/(alpha - d - beta + d^2/(d+beta))) * (alpha/(alpha - d - beta + d^2/(d+beta)));
            (d + beta - (d^2)/(d+beta))/(alpha - d - beta + d^2/(d+beta)), 0, ( 1-1/K * (d + beta - d^2/(d+beta))/(alpha - d - beta + d^2/(d+beta))) * (alpha/(alpha - d - beta + d^2/(d+beta))), d/(d+beta)*(1 - 1/K*(d + beta - d^2/(d+beta))/(alpha - d - beta + (d^2)/(d+beta))) * ( alpha/(alpha - d - beta + d^2/(d+beta)));
            alpha/(alpha-beta),alpha/(alpha-beta),alpha/(alpha-beta)*(1-1/K*(beta/(alpha-beta))),alpha/(alpha-beta)*(1-1/K*(beta/(alpha-beta)))];

            % Finding the non-symmetric coexistence equilibrium
            extra_eqa = find_eqa(alpha,beta,K,d);
            num_rows = size(extra_eqa,1);

            eqa = [eqa; extra_eqa];

            unstable = false(1,4+num_rows);

            % Determining whether any of the eigenvalues are negative
            for i = 1:(4+num_rows) % for each equilibrium
                infinite = false;
                eigenvalues = find_eig(eqa(i,1),eqa(i,2),eqa(i,3),eqa(i,4),alpha,beta,K,d);

                if any(~isfinite(eqa(i)))
                    unstable(i) = true;

                else 
                    if any(real(eigenvalues) > 1e-6)
                        unstable(i) = true;
                    end

                end

            end

            % Produce colour plots
            if unstable(1) == false && all(unstable(2:4+num_rows)) == true % predator-free
                colourMatrix{beta_index, K_index} = [1, 1, 0]; % yellow 

            elseif unstable(1) == true && unstable(2) == true && unstable(3) == true && any(~unstable(4:4+num_rows)) == true % coexistence 
                colourMatrix{beta_index, K_index} = [0, 0, 1]; % blue               

            elseif unstable(1) == true && (unstable(2) == false || unstable(3) == false) && all(unstable(4:4+num_rows)) == true % partial-prey-death 
                colourMatrix{beta_index, K_index} = [0, 1, 0]; % green

            elseif unstable(1) == true && unstable(2) == true && unstable(3) == true && all(unstable(4:4+num_rows)) == true % non-equilibrium
                colourMatrix{beta_index, K_index} = [1, 0, 0]; % red

            else
                colourMatrix{beta_index, K_index} = [1, 0, 1]; % magenta

            end
            
        end
    end
    
    
    % Extract coordinates from the parameter matrix
    coordinates = cell2mat(parameterMatrix(:));
    
    % Extract colours from the colour matrix
    colours = cell2mat(colourMatrix(:));
    
    % Create a scatter plot of the colours
    figure;
    scatter(coordinates(:, 1), coordinates(:, 2), 50, colours, 'filled');

    title_str = sprintf('$d = %g, \\alpha = 16$', d);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 25);  

    xlabel('\beta', 'FontSize', 25);
    ylabel('K', 'FontSize', 19);
    
    % Adjust plot limits to avoid dots going over the axis
    axis([min(betaList) max(betaList) min(KList) max(KList)]);
    
end

