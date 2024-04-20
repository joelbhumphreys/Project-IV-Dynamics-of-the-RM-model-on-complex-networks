function [ A ] = find_ad_matrix(num_patches)

    % Defines an adjacency matrix with specified number of patches,
    % connecting each patch to its closest three neighbours

    for i = 1:num_patches
        for j = 1:num_patches
            if ((j == i-1) || (j == i+1) || (j == i-2) || (j == i+2) || (j == i-3) || (j == i+3)) 
                A(i,j) = 1;
            else
                A(i,j) = 0;
            end
        end
    end

    A(1,num_patches) = 1;
    A(1,num_patches-1) = 1;    
    A(1,num_patches-2) = 1;        
    A(num_patches,1) = 1;
    A(num_patches-1,1) = 1;
    A(num_patches-2,1) = 1;

    A(2,num_patches) = 1;
    A(2,num_patches-1) = 1;    
    A(num_patches,2) = 1;
    A(num_patches-1,2) = 1;  

    A(3,num_patches) = 1;    
    A(num_patches,3) = 1;  

end






