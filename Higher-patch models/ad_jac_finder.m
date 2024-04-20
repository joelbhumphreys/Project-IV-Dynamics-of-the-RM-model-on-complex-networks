function eigenvalues = ad_jac_finder(num_patches,eqa,alpha, beta, K, d)
    eqa_row = eqa;

    % generate a matrix of zeros of the correct size
    Jacobian = zeros(num_patches*2, num_patches*2); 

    % fill the matrix with the components of the matrix by partioning it
    % into five components
    for i = 1:num_patches
        
        Jacobian(i,i) = computeC1(K, eqa_row(i), eqa_row(i+num_patches));
        Jacobian(i,i+num_patches) = computeC2(eqa_row(i));
        Jacobian(i+num_patches,i) = computeC3(alpha,eqa_row(i),eqa_row(i+num_patches));
        Jacobian(i+num_patches,i+num_patches) = computeC4(alpha,beta,eqa_row(i)); 
        Jacobian(i+num_patches,i+num_patches) = Jacobian(i+num_patches,i+num_patches) - d * indicator(i, num_patches);
        
        for j = 1:num_patches
            if (i == j + 1 || i == j - 1)
                Jacobian(i+num_patches,j+num_patches) = d;
            end
        end

    end

    % find the eigenvalues
    eigenvalues = eig(Jacobian);

end

