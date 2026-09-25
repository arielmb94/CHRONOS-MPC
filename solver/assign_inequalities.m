function mpc = assign_inequalities(mpc)

nx = mpc.nx;
nu = mpc.nu;
ny_0 = mpc.ny_0;
ny = mpc.ny;
ny_ter = mpc.ny_ter;
nh_0 = mpc.nh_0;
nh = mpc.nh;
nh_ter = mpc.nh_ter;

%% k = 0
v_rows_0 = [];

start_index = 1;

if ~isempty(mpc.has_u_cnstr)
    if ~isempty(mpc.u_cnstr.min_limit)
        % Lower-bound inequality rows.
        row = start_index:start_index+nu-1;

        mpc.u_cnstr.min_ineqRow_0 = row;

        start_index = start_index + nu;
    end
    if ~isempty(mpc.u_cnstr.max_limit)
        % Upper-bound inequality rows.
        row = start_index:start_index+nu-1;

        mpc.u_cnstr.max_ineqRow_0 = row;

        start_index = start_index + nu;
    end
end

if ~isempty(mpc.has_du_cnstr)
    if ~isempty(mpc.du_cnstr.min_limit)
        % Lower-bound inequality rows.
        row = start_index:start_index+nu-1;

        mpc.du_cnstr.min_ineqRow_0 = row;

        start_index = start_index + nu;
    end
    if ~isempty(mpc.du_cnstr.max_limit)
        % Upper-bound inequality rows.
        row = start_index:start_index+nu-1;

        mpc.du_cnstr.max_ineqRow_0 = row;

        start_index = start_index + nu;
    end
end

start_index_v = 1;
if ~isempty(mpc.has_y_cnstr)
    if ~isempty(mpc.y_cnstr.min_limit) && ~isempty(mpc.y_cnstr.use_k0)
        % Lower-bound inequality rows.
        row = start_index:start_index+ny_0-1;
        row_v = start_index_v:start_index_v+ny_0-1;

        mpc.y_cnstr.min_ineqRow_0 = row;
        mpc.y_cnstr.min_row_v_0 = row_v;

        v_rows_0 = [v_rows_0; row'];

        start_index = start_index + ny_0;
        start_index_v = start_index_v + ny_0;
    end
    if ~isempty(mpc.y_cnstr.max_limit) && ~isempty(mpc.y_cnstr.use_k0)
        % Upper-bound inequality rows.
        row = start_index:start_index+ny_0-1;
        row_v = start_index_v:start_index_v+ny_0-1;

        mpc.y_cnstr.max_ineqRow_0 = row;
        mpc.y_cnstr.max_row_v_0 = row_v;

        v_rows_0 = [v_rows_0; row'];

        start_index = start_index + ny_0;
        start_index_v = start_index_v + ny_0;
    end
end

if ~isempty(mpc.has_h_cnstr)
    if ~isempty(mpc.h_cnstr.min_limit) && ~isempty(mpc.h_cnstr.use_k0)
        % Lower-bound inequality rows.
        row = start_index:start_index+nh_0-1;
        row_v = start_index_v:start_index_v+nh_0-1;

        mpc.h_cnstr.min_ineqRow_0 = row;
        mpc.h_cnstr.min_row_v_0 = row_v;

        v_rows_0 = [v_rows_0; row'];

        start_index = start_index + nh_0;
        start_index_v = start_index_v + nh_0;
    end
    if ~isempty(mpc.h_cnstr.max_limit) && ~isempty(mpc.h_cnstr.use_k0)
        % Upper-bound inequality rows.
        row = start_index:start_index+nh_0-1;
        row_v = start_index_v:start_index_v+nh_0-1;

        mpc.h_cnstr.max_ineqRow_0 = row;
        mpc.h_cnstr.max_row_v_0 = row_v;

        v_rows_0 = [v_rows_0; row'];

        start_index = start_index + nh_0;
        start_index_v = start_index_v + nh_0;
    end
end
mpc.v_rows_0 = v_rows_0;

%% k = 1:N-1
v_rows_k = [];
start_index = 1;
start_index_v = 1;

    if ~isempty(mpc.has_s_cnstr)
        if ~isempty(mpc.s_cnstr.min_limit)
            % Lower-bound inequality rows.
            row = start_index:start_index+nx-1;
            row_v = start_index_v:start_index_v+nx-1;

            mpc.s_cnstr.min_ineqRow_k = row;
            mpc.s_cnstr.min_row_v_k = row_v;
           
            v_rows_k = [v_rows_k; row'];

            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
        if ~isempty(mpc.s_cnstr.max_limit)
            % Upper-bound inequality rows.
            row = start_index:start_index+nx-1;
            row_v = start_index_v:start_index_v+nx-1;

            mpc.s_cnstr.max_ineqRow_k = row;
            mpc.s_cnstr.max_row_v_k = row_v;

            v_rows_k = [v_rows_k; row'];

            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
    end

    if ~isempty(mpc.has_u_cnstr)
        if ~isempty(mpc.u_cnstr.min_limit)
            % Lower-bound inequality rows.
            row = start_index:start_index+nu-1;

            mpc.u_cnstr.min_ineqRow_k = row;

            start_index = start_index + nu;
        end
        if ~isempty(mpc.u_cnstr.max_limit)
            % Upper-bound inequality rows.
            row = start_index:start_index+nu-1;

            mpc.u_cnstr.max_ineqRow_k = row;

            start_index = start_index + nu;
        end
    end

    if ~isempty(mpc.has_du_cnstr)
        if ~isempty(mpc.du_cnstr.min_limit)
            % Lower-bound inequality rows.
            row = start_index:start_index+nu-1;

            mpc.du_cnstr.min_ineqRow_k = row;

            start_index = start_index + nu;
        end
        if ~isempty(mpc.du_cnstr.max_limit)
            % Upper-bound inequality rows.
            row = start_index:start_index+nu-1;

            mpc.du_cnstr.max_ineqRow_k = row;

            start_index = start_index + nu;
        end
    end

    if ~isempty(mpc.has_y_cnstr)
        if ~isempty(mpc.y_cnstr.min_limit)
            % Lower-bound inequality rows.
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+ny-1;

            mpc.y_cnstr.min_ineqRow_k = row;
            mpc.y_cnstr.min_row_v_k = row_v;

            v_rows_k = [v_rows_k; row'];

            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
        if ~isempty(mpc.y_cnstr.max_limit)
            % Upper-bound inequality rows.
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+ny-1;

            mpc.y_cnstr.max_ineqRow_k = row;
            mpc.y_cnstr.max_row_v_k = row_v;

            v_rows_k = [v_rows_k; row'];

            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
    end

    if ~isempty(mpc.has_h_cnstr)
        if ~isempty(mpc.h_cnstr.min_limit)
            % Lower-bound inequality rows.
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.min_ineqRow_k = row;
            mpc.h_cnstr.min_row_v_k = row_v;

            v_rows_k = [v_rows_k; row'];

            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
        if ~isempty(mpc.h_cnstr.max_limit)
            % Upper-bound inequality rows.
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.max_ineqRow_k = row;
            mpc.h_cnstr.max_row_v_k = row_v;

            v_rows_k = [v_rows_k; row'];

            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
    end

mpc.v_rows_k = v_rows_k;

%% k = N
start_index = 1;

if ~isempty(mpc.has_s_cnstr)
    if ~isempty(mpc.s_cnstr.min_limit)
        % Lower-bound inequality rows.
        row = start_index:start_index+nx-1;

        mpc.s_cnstr.min_ineqRow_ter = row;

        start_index = start_index + nx;
    end
    if ~isempty(mpc.s_cnstr.max_limit)
        % Upper-bound inequality rows.
        row = start_index:start_index+nx-1;

        mpc.s_cnstr.max_ineqRow_ter = row;

        start_index = start_index + nx;
    end
end

if ~isempty(mpc.has_y_cnstr)
    if ~isempty(mpc.y_cnstr.min_limit) && ~isempty(mpc.y_cnstr.use_ter)
        % Lower-bound inequality rows.
        row = start_index:start_index+ny_ter-1;

        mpc.y_cnstr.min_ineqRow_ter = row;

        start_index = start_index + ny_ter;
    end
    if ~isempty(mpc.y_cnstr.max_limit) && ~isempty(mpc.y_cnstr.use_ter)
        % Upper-bound inequality rows.
        row = start_index:start_index+ny_ter-1;

        mpc.y_cnstr.max_ineqRow_ter = row;

        start_index = start_index + ny_ter;
    end
end

if ~isempty(mpc.has_h_cnstr)
    if ~isempty(mpc.h_cnstr.min_limit) && ~isempty(mpc.h_cnstr.use_ter)
        % Lower-bound inequality rows.
        row = start_index:start_index+nh_ter-1;

        mpc.h_cnstr.min_ineqRow_ter = row;

        start_index = start_index + nh_ter;
    end
    if ~isempty(mpc.h_cnstr.max_limit) && ~isempty(mpc.h_cnstr.use_ter)
        % Upper-bound inequality rows.
        row = start_index:start_index+nh_ter-1;

        mpc.h_cnstr.max_ineqRow_ter = row;

        start_index = start_index + nh_ter;
    end
end

end
