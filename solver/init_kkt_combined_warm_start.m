function mpc_maps = init_kkt_combined_warm_start(phi, n)
    m = size(phi,1) - n;
    H  = phi(1:n, 1:n);
    AT = phi(1:n, n+1:end);
    
    L_H = chol(H, 'lower');
    opts_LT.LT = true;
    U   = linsolve(L_H, AT, opts_LT);
    S   = U' * U;
    L_S = chol(S, 'lower');
    
    [mpc_maps.L_H_idx, mpc_maps.L_H_len] = createMapMat(L_H);
    [mpc_maps.L_S_idx, mpc_maps.L_S_len] = createMapMat(L_S);
    
    [mpc_maps.fwd_H_idx, mpc_maps.fwd_H_len] = createRowMapMat(L_H);
    [mpc_maps.bwd_H_idx, mpc_maps.bwd_H_len] = createRowMapMat(L_H'); 
end

function [map_idx, map_len] = createMapMat(L)
    n = size(L, 1);
    map_idx = zeros(n, n);
    map_len = zeros(n, 1);
    for j = 1:n
        idx = find(L(j+1:end, j) ~= 0) + j;
        map_len(j) = length(idx);
        if map_len(j) > 0
            map_idx(1:map_len(j), j) = idx;
        end
    end
end

function [map_idx, map_len] = createRowMapMat(L)
    n = size(L, 1);
    map_idx = zeros(n, n);
    map_len = zeros(n, 1);
    for i = 1:n
        idx = find(L(i, :) ~= 0);
        idx(idx == i) = [];
        map_len(i) = length(idx);
        if map_len(i) > 0
            map_idx(1:map_len(i), i) = idx(:); 
        end
    end
end