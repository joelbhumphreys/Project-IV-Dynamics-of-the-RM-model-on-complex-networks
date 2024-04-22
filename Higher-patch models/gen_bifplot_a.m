function gen_bif_plot_a()

    num_patches = 7;
    timeLimit = 8;
    
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

                % Find the equilibrium points                
                eqa = gen_find_eqa(timeLimit,num_patches,alpha,beta,K,d);
                coex_eqa(1:num_patches) = beta/(alpha-beta);
                coex_eqa(num_patches+1:2*num_patches) = alpha/(alpha-beta)*(1-1/K*beta/(alpha-beta));
                pred_free_eqa = find_pred_death(num_patches, K);
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
                
                % Classify coexistence
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
    
        title_str = sprintf('$d = %g, \\alpha = 16$', d);    
        title(title_str, 'Interpreter', 'latex', 'FontSize', 25); 
        xlabel('\beta', 'FontSize', 25);
        ylabel('K', 'FontSize', 19);
    
        % Adjust plot limits to avoid dots going over the axis
        axis([min(betaList) max(betaList) min(KList) max(KList)]);
        
        
    end
    
end


