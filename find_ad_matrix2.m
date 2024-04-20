function [ A ] = find_ad_matrix2(num_patches)

    % Generates the fully-connected adjacency matrix

    A = ones(num_patches);
    
    for i = 1:num_patches
        A(i,i) = 0;
    end
    
end






