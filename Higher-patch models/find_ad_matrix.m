function [ A ] = find_ad_matrix(num_patches,num_connec)

    % Generates adjacency matrices for ring network models with a specified
    % number of patches and edges

    for i = 1:num_patches % iterate over patches
        for k = 1:num_connec % iterate over edges
            
            if i-k < 1
                A(i,i+k) = 1;                
                A(i,num_patches + i - k) = 1;

            elseif i+k > num_patches
                A(i,i-k) = 1;                                
                A(i, i + k - num_patches) = 1;

            else
                A(i,i-k) = 1;
                A(i,i+k) = 1; 

            end
        end
    end
    
end

