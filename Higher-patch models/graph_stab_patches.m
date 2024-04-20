d = 1; % set diffusion rate
num_patches = 1; % set number of parameters

num_stab_pnts = zeros(1,num_patches); % initialise vector to record results

% Iterate quan_stab_a function for graphs with one to seven patches  
for patches = 1:num_patches
    num_stab_pnts(patches) = quan_stab_a(patches,d);
end

% Plot the line graph 
figure;
plot(num_patches, num_stab_pnts, '-o', 'LineWidth', 2, 'MarkerSize', 8);

xlabel('Num patches', 'FontSize', 14);
ylabel('Num points where partial-prey-death eqm is stable', 'FontSize', 14);
title('$\alpha = 1, \beta = 1, K = 1, d = 1$', 'Interpreter', 'latex', 'FontSize', 16);

grid on;


