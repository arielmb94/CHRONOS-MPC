function validate_matrix(val, expected_rows, expected_columns, var_name, allow_stage_scalar)
% VALIDATE_MATRIX Validate the per-stage shape of a staged matrix.

    if nargin < 5
        allow_stage_scalar = false;
    end

    % Empty is allowed for optional matrix inputs.
    if isempty(val)
        return;
    end

    % Validate only the first supplied page. The remaining pages are handled
    % by the existing horizon fill logic.
    stage_val = val(:,:,1);

    if allow_stage_scalar && isscalar(stage_val)
        return;
    end

    % A scalar is accepted by the default policy only for a 1-by-1 matrix,
    % because the scalar stage has exactly those per-stage dimensions.
    if size(stage_val,1) ~= expected_rows || size(stage_val,2) ~= expected_columns
        error('CHRONOS:DimensionMismatch', ...
            'Input "%s" must have per-stage dimensions %d x %d; supplied per-stage dimensions are %d x %d.', ...
            var_name, expected_rows, expected_columns, ...
            size(stage_val,1), size(stage_val,2));
    end
end
