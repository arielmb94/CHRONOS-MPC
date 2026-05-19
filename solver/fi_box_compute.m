% Computes fi (fi<0) for box constraints
function cnstr = fi_box_compute(cnstr,x,N,n)

% Let H be the number of time steps (horizon)
H = N / n;

if cnstr.min_limit
    % Compute f_i for ALL time instances simultaneously
    cnstr.fi_min_x0 = repmat(cnstr.min, H, 1) - x;
end

if cnstr.max_limit
    % Compute f_i for ALL time instances simultaneously
    cnstr.fi_max_x0 = x - repmat(cnstr.max, H, 1);
end

end