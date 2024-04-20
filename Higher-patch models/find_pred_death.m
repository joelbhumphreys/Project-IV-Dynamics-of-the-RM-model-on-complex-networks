function resultMatrix = find_pred_death(num_patches, K)

    possibilities = [0; K];
    combinations = generateCombinations(possibilities, num_patches);

    combinations = combinations(~all(combinations == 0, 2), :);
    resultMatrix = [combinations, zeros(size(combinations, 1), num_patches)];

end

function combinations = generateCombinations(possibilities, num_patches)
    % Recursive function to generate all combinations of 0 and K

    if num_patches == 1
        combinations = possibilities;
    else
        smallerCombinations = generateCombinations(possibilities, num_patches-1);
        combinations = [];
        for i = 1:numel(possibilities)
            combinations = [combinations; possibilities(i)*ones(size(smallerCombinations, 1), 1), smallerCombinations];
        end
    end
    
end


