function rho = calc_rho_stirring_tank(x, k, M)
    % Extract the states
    ck = x(1);
    vk = x(2);
    
    % Avoid division by zero
    vk = max(vk, 0.1); 
    ck = max(ck, 0.0);
    ck = min(ck, 1.0);
    
    % Compute the schedulling variables
    rho1 = k * exp(-M/vk);
    rho2 = k * ck * M * exp(-M/vk) / (vk^2);
    rho3 = vk;
    
    % Column vector
    rho = [rho1; rho2; rho3];
end