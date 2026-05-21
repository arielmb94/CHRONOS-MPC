function sigma_k = calc_jacob_stirring_tank(x, k_cte, M)
    % Extract the states
    ck = x(1);
    vk = max(x(2), 0.1); % To avoid division by zero
    
    % Diff of rho1 = k * exp(-M/vk)
    d_rho1_dc = 0;
    d_rho1_dv = k_cte * M / (vk^2) * exp(-M/vk);
    
    % Diff of rho2 = k * ck * M * exp(-M/vk) / (vk^2)
    d_rho2_dc = k_cte * M / (vk^2) * exp(-M/vk);
    d_rho2_dv = k_cte * ck * M * exp(-M/vk) * (M/(vk^4) - 2/(vk^3));
    
    % Diff of rho3 = vk
    d_rho3_dc = 0;
    d_rho3_dv = 1;
    
    % Compute the Jacobian sigma_k (n_rho x nx) -> (3 x 2)
    sigma_k = [d_rho1_dc, d_rho1_dv;
               d_rho2_dc, d_rho2_dv;
               d_rho3_dc, d_rho3_dv];
end