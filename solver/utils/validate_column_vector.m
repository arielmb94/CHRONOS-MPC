function validate_column_vector(val, expected_length, var_name, allow_stage_scalar)
% VALIDATE_COLUMN_VECTOR Validate the per-stage shape of a staged vector.

    if nargin < 4
        allow_stage_scalar = false;
    end

    % Empty is allowed (means constraint is disabled)
    if isempty(val)
        return;
    end
    
    % Scalars are allowed (will be expanded later)
    if isscalar(val)
        return;
    end

    % Validate only the first supplied stage. The remaining stages are handled
    % by the existing horizon fill logic.
    stage_val = val(:,1);

    if allow_stage_scalar && isscalar(stage_val)
        return;
    end

    % If the per-stage shape is wrong, throw a custom error.
    if size(stage_val,1) ~= expected_length || size(stage_val,2) ~= 1
        error('CHRONOS:DimensionMismatch', ...
            'Input "%s" must have per-stage dimensions %d x 1; supplied per-stage dimensions are %d x %d.', ...
            var_name, expected_length, size(stage_val,1), size(stage_val,2));
    end
end
