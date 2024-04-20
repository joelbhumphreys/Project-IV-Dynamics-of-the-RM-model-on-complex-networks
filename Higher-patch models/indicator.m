function edges = indicator(i, num_patches)
    % Check if i is equal to 1 or num_patches
    if (i == 1) || (i == num_patches)
        edges = 1;
    else
        edges = 2;
    end
end