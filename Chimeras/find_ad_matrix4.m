function [ A ] = find_ad_matrix3(num_patches)

    % Defines an adjacency matrix with specified number of patches, 
    % connecting each patch to all the others, except the three furthest away 

    % Only works for models with an even number of patches
    
    A = ones(num_patches);

    for i = 1:num_patches
        for j = 1:num_patches
            A(i,i) = 0;
        end
    end

    for i = 1:num_patches
        for j = 1:(num_patches)

            if (j == i-(num_patches/2)) || (j == i+(num_patches/2) || j == i-(num_patches/2)+1) || (j == i+(num_patches/2)+1 || j == i-(num_patches/2)-1) || (j == i+(num_patches/2)-1)
                A(i,j) = 0;
            end

        end
    end  

     
end




