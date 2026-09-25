function M_full = fill_mat(M_full, M, fill_style)

step_m = size(M, 3);
step_full = size(M_full, 3);

% Clamp stage copying to the target buffer length
step_copy = min(step_m, step_full);

% Copy available stages up to the target limit
M_full(:,:, 1:step_copy) = M(:,:, 1:step_copy);

% Fill remaining stages if input is shorter than target horizon
if step_copy < step_full
    if fill_style == 0
        % Zero fill remaining stages
        M_full(:,:, step_copy+1:step_full) = 0;
    else
        % Repeat last available stage
        for k = step_copy+1:step_full
            M_full(:,:, k) = M(:,:, step_copy);
        end
    end
end
end