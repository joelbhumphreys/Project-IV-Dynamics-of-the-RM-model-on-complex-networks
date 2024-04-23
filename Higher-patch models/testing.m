% This function was used to test individual parameter regimes

timeLimit = 8; % set time limit for finding eqa
initial_cond = 14/24; % parameter used in finding eqa function

num_patches = 7; % set the number of patches

% Set the parameter values
alpha = 16;
beta = 8;
K = 9;
d = 20;

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
unstable = ones(1,num_rows);

for i = 1:(num_rows) % for each equilibrium
    eigenvalues = ad_jac_finder(num_patches,ordered_eqa(i,:),alpha, beta, K, d);

    if all(isAlways(real(eigenvalues) < 0))
        unstable(i) = false;
    end

end

% Produce colour plots
if any(unstable(1:endmark1) == false) && all(unstable(endmark1+1:num_rows) == true) % full extinction 
    colour = [255, 165, 0]; % orange

elseif any(unstable(endmark1+1:endmark2) == false) && all(unstable(1:endmark1) == true) && all(unstable(endmark2+1:num_rows) == true) % predator-free
    colour = [1, 1, 0]; % yellow

elseif any(unstable(endmark2+1:endmark3) == false) && all(unstable(1:endmark2) == true) && all(unstable(endmark3+1:num_rows) == true) % coexistence
    colour = [0, 0, 1]; % blue             

elseif any(unstable(endmark3+1:num_rows) == false) && all(unstable(1:endmark3) == true) % partial-prey-death
    colour = [0, 1, 0]; % green 

elseif all(unstable == true) % non-equilibrium                    
    colour = [1, 0, 0]; % red
    
else
    colour = [1, 0, 1]; % magenta

end

disp(ordered_eqa);
disp(colour);



