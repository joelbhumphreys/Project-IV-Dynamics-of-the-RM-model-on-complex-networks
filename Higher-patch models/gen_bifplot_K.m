function gen_bif_plot_K()

    num_patches = 2;
    timeLimit = 0.1;
    
    % Iterate over d        
    for d = [1, 5, 20]
        % Set the ranges for the parameter space                                
        alphaList = 22.01:0.13:47.01;
        betaList = 15.01:0.13:35.01;
        K = 16;  
    
        parameterMatrix = cell(length(alphaList), length(betaList));
        
        for i = 1:length(alphaList)
            for j = 1:length(betaList)
                parameterMatrix{i, j} = [alphaList(i), betaList(j)];
            end
        end
        
        % Construct matrix to hold the colour outputs                                    
        colourMatrix = cell(length(alphaList), length(betaList));
        
        % Iterate over the parameter space                                         
        alpha_index = 0;
        
        for alpha = 22.01:0.13:47.01
            beta_index = 0;
            alpha_index = alpha_index + 1;
    
            for beta = 15.01:0.13:35.01
                beta_index = beta_index + 1;
    
                % Find the equilibrium points                
                eqa = gen_find_eqa(timeLimit,num_patches,alpha,beta,K,d);
                coex_eqa(1:num_patches) = beta/(alpha-beta);
                coex_eqa(num_patches+1:2*num_patches) = alpha/(alpha-beta)*(1-1/K*beta/(alpha-beta));
                pred_free_eqa(1:num_patches) = K;
                pred_free_eqa(num_patches+1:2*num_patches) = 0;                
                combined_eqa = vertcat(eqa, coex_eqa, pred_free_eqa);

                num_rows = size(combined_eqa,1);
                ordered_eqa = zeros(num_rows,num_patches*2);
    
                filled = 0;          
    
                % Classify full extinction
                num_removed = 0; % Initialize the counter for removed elements
                for j = num_rows:-1:1  % Backward loop
                    if all(combined_eqa(j,:) < 1e-3)
                        filled = filled + 1;                        
                        ordered_eqa(filled,:) = combined_eqa(j,:); 
                        combined_eqa(j,:) = [];  % Remove the element from eqa
                        num_removed = num_removed + 1; % Increment the counter
                    end
                end
                endmark1 = filled;
                num_rows_temp = num_rows - num_removed; % Update num_rows
                
                % Classify predator-free
                num_removed = 0; % Reset the counter
                for j = num_rows_temp:-1:1  % Backward loop
                    if all(combined_eqa(j, num_patches+1:2*num_patches) < 1e-3)
                        filled = filled + 1;                        
                        ordered_eqa(filled,:) = combined_eqa(j,:);                 
                        combined_eqa(j,:) = [];  
                        num_removed = num_removed + 1; 
                    end
                end
                endmark2 = filled;
                num_rows_temp = num_rows_temp - num_removed; 
                
                % Classify full coexistence
                num_removed = 0; % Reset the counter
                for j = num_rows_temp:-1:1  % Backward loop
                    if all(combined_eqa(j,:) > 1e-3)
                        filled = filled + 1;                        
                        ordered_eqa(filled,:) = combined_eqa(j,:);                       
                        combined_eqa(j,:) = [];  
                        num_removed = num_removed + 1; 
                    end
                end
                endmark3 = filled;
                num_rows_temp = num_rows_temp - num_removed; 
                                  
    
                % Classify partial-prey-death
                for j = 1:num_rows_temp
                    filled = filled + 1;                        
                    ordered_eqa(filled,:) = combined_eqa(j,:); 
                end
    
    
                % Determine the stability 
                unstable = false(1,num_rows);
    
                for i = 1:(num_rows) % for each equilibrium
                    eigenvalues = ad_jac_finder(num_patches,ordered_eqa(i,:),alpha, beta, K, d);
    
                    if any(isAlways(real(eigenvalues) > 1e-6))
                        unstable(i) = true;
                    end
    
                end
    
                % Produce colour plots
                if any(unstable(1:endmark1) == false) && all(unstable(endmark1+1:num_rows) == true) % full extinction 
                    colourMatrix{beta_index, K_index} = [255, 165, 0]; % orange
    
                elseif any(unstable(endmark1+1:endmark2) == false) && all(unstable(1:endmark1) == true) && all(unstable(endmark2+1:num_rows) == true) % predator-free
                    colourMatrix{beta_index, K_index} = [1, 1, 0]; % yellow
    
                elseif any(unstable(endmark2+1:endmark3) == false) && all(unstable(1:endmark2) == true) && all(unstable(endmark3+1:num_rows) == true) % coexistence
                    colourMatrix{beta_index, K_index} = [0, 0, 1]; % blue             
    
                elseif any(unstable(endmark3+1:num_rows) == false) && all(unstable(1:endmark3) == true) % partial-prey-death
                    colourMatrix{beta_index, K_index} = [0, 1, 0]; % green   
    
                elseif all(unstable == true) % non-equilibrium                    
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
        title_str = sprintf('$d = %g, K = 16$', d);    
        title(title_str, 'Interpreter', 'latex', 'FontSize', 25); 
        xlabel('\alpha', 'FontSize', 25);
        ylabel('\beta', 'FontSize', 19);
            
        % Adjust plot limits to avoid dots going over the axis
        axis([min(alphaList) max(alphaList) min(betaList) max(betaList)]);
        
    end
    
end


