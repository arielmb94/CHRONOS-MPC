function validate_column_vector(val, expected_length, var_name)
% VALIDATE_COLUMN_VECTOR Throws a clear error if the input is not a scalar,
% empty, or a column vector of the exact expected length.

    % Empty is allowed (means constraint is disabled)
    if isempty(val)
        return;
    end
    
    % Scalars are allowed (will be expanded later)
    if isscalar(val)
        return;
    end
    
    % If it's not a column vector, or the length is wrong, throw a custom error
    if ~iscolumn(val) || length(val) ~= expected_length
        error('CHRONOS:DimensionMismatch', ...
            'Error in constraint setup: "%s" must be a scalar or a column vector of size %d x 1. You provided a size of %d x %d.', ...
            var_name, expected_length, size(val,1), size(val,2));
    end
end