function v_full = fill_vec(v_full, v, fill_style)
step_v = size(v,2);
step_full = size(v_full,2);

step_copy = min(step_v, step_full);

% Copy available stages
v_full(:,1:step_copy) = v(:,1:step_copy);

% Fill remaining stages
if step_copy < step_full
    if fill_style == 0
        % Zero fill
        v_full(:,step_copy+1:step_full) = 0;
    else
        % Repeat last stage
        for i = step_copy+1:step_full
            v_full(:,i) = v(:,step_copy);
        end
    end
end
end