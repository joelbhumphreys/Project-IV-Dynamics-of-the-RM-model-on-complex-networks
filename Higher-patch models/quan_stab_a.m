function num_stab_pnts = quan_stab_a(num_patches,d)
    % set the fixed parameter
    alpha = 16;  
        
    % Iterate over the parameter space                                 
    beta_index = 0;
    
    for beta = 4.01:0.33:16.01
        K_index = 0;
        beta_index = beta_index + 1;

        for K = 0.01:0.33:14.01
            K_index = K_index + 1;

            eqa = gen_find_eqa(num_patches,alpha,beta,K,d);
            eqa2 = find_pred_death(num_patches, K);
            combined_eqa = vertcat(eqa, eqa2);
            
            ppd_eqa = [];
            num_ppd_eqa = 0;

            % Identify partial-prey-death        
            for j = 1:size(combined_eqa,1) % for each equilibrium
                if all(combined_eqa(j, num_patches+1:2*num_patches) > 1e-6) && any(combined_eqa(j, 1:num_patches) <= 1e-6)
                    num_ppd_eqa = num_ppd_eqa + 1;
                    ppd_eqa(num_ppd_eqa,:) = combined_eqa(j,:); 
                    disp('pop');
                end
            end

            % Determine the stability 
            num_stab_pnts = 0;

            for i = 1:num_ppd_eqa % for each equilibrium
                eigenvalues = ad_jac_finder(num_patches,ppd_eqa(i,:),alpha, beta, K, d);

                if all(isAlways(real(eigenvalues) < 0))
                    num_stab_pnts = num_stab_pnts + 1; % count instances of partial-prey-death
                end

            end

        end
    end
  
end



