function C1 = computeC1(K, u0, v0)
    C1 = 1 - 2/K * u0 - v0/(1 + u0)^2;
end

