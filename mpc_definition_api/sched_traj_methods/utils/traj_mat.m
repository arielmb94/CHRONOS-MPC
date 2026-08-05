%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   M_full = traj_mat(compute_M_fun, P_traj, n_rho, N)
%
% Converts a concatenated scheduling parameter trajectory into a 3D matrix 
% array along the prediction horizon.
%
% In:
%   - compute_M_fun: function handle (@(rho)), computes the 2D matrix for a 
%   single time step given the scheduling parameters rho.
%   - P_traj: (N * n_rho) column vector, the scheduling trajectory.
%   - n_rho: integer scalar, number of scheduling variables.
%   - N: integer scalar, prediction horizon length.
%
% Out:
%   - M_full: 3D array of size [nx, nx, N] containing the evaluated matrices.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M_full = traj_mat(compute_M_fun, P_traj, n_rho, N)
    
    % Evaluate the first element to preallocate the 3D array correctly
    rho_initial = P_traj(1:n_rho);
    M_initial = compute_M_fun(rho_initial);
    
    [rows, cols] = size(M_initial);
    M_full = zeros(rows, cols, N);
    
    M_full(:,:,1) = M_initial;
    
    % Loop through the rest of the horizon
    for j = 2:N
        rho_j = P_traj((j-1)*n_rho + 1 : j*n_rho);
        M_full(:,:,j) = compute_M_fun(rho_j);
    end
end