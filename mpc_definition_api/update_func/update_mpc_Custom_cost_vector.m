% UPDATE_MPC_CUSTOM_COST_VECTOR Update the custom-cost signal model.
%
%   mpc = UPDATE_MPC_CUSTOM_COST_VECTOR(mpc, Cz, Dz, Dsuz, Ddz) updates
%   selected coefficients in
%
%       z_k = Cz_k*s_k + Dz_k*u_k + Dsuz_k*su_k + Ddz_k*dz_k.
%
%   Use [] to leave a coefficient unchanged. The signal must first be
%   defined with INIT_MPC_CUSTOM_COST. This updater cannot change the number
%   of custom costs or add a coefficient term that was not enabled during
%   initialization.
%
%   Coefficients may be constant matrices or contain L horizon stages in
%   their third dimension. If L < N, the last supplied stage is reused for
%   the remaining stages. Here, su_k is the control action preceding u_k,
%   and dz_k is the dedicated fixed known input supplied to MPC_SOLVE.
% 
%   The existing Qz and qz weights are applied to the updated
%   signal model.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Cz      - Optional state coefficient, size nz-by-nx or nz-by-nx-by-L.
%     Dz      - Optional control coefficient, size nz-by-nu or
%               nz-by-nu-by-L.
%     Dsuz    - Optional previous-control coefficient, size nz-by-nu or
%               nz-by-nu-by-L.
%     Ddz     - Optional fixed-known-input coefficient, size nz-by-ndz or
%               nz-by-ndz-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only Cz and Dz:
%
%       mpc = update_mpc_Custom_cost_vector(mpc, Cz, Dz, [], []);
function mpc = update_mpc_Custom_cost_vector(mpc,Cz,Dz,Dsuz,Ddz)

if ~isempty(Cz) && ~isempty(mpc.z_use_s)
    if ~isempty(mpc.quad_custom_cost)
        mpc.update_customcost_quad = true;
    end
    if ~isempty(mpc.lin_custom_cost)
        mpc.update_customcost_lin = true;
    end

<<<<<<< HEAD
    len_Cz = size(Cz,3);
    mpc.Cz(:,:,:) = fill_mat(mpc.Cz, Cz, 1);
    if ~isempty(mpc.z_use_ter)
        ter_stage = len_Cz;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.Cz_ter(:,:) = Cz(mpc.z_rows_ter,:,ter_stage);
    end
    if ~isempty(mpc.z_use_k0), mpc.Cz_0(:,:) = Cz(mpc.z_rows_k0,:,1); end

end

if ~isempty(Dz) && ~isempty(mpc.z_use_u)
    if ~isempty(mpc.quad_custom_cost)
        mpc.update_customcost_quad = 1;
    end
    if ~isempty(mpc.lin_custom_cost)
        mpc.update_customcost_lin = 1;
    end
    mpc.Dz(:,:,:) = fill_mat(mpc.Dz, Dz, 1);
    if ~isempty(mpc.z_use_k0)
        mpc.Dz_0(:,:) = Dz(mpc.z_rows_k0,:,1);
    end
end

if ~isempty(Dsuz) && ~isempty(mpc.z_use_su)
    if ~isempty(mpc.quad_custom_cost)
        mpc.update_customcost_quad = true;
=======
if ~isempty(Cz)
    update_grad = 1;

    len_Cz = size(Cz,3);
    if len_Cz < mpc.N
        mpc.Cz(:,:,:) = fill_mat(mpc.Cz, Cz, 1);
        if mpc.z_use_ter, mpc.Cz_ter(:,:) = mpc.Cz(mpc.z_rows_ter,:,mpc.N-1); end
    else
        mpc.Cz(:,:,:) = Cz(:,:,1:mpc.N-1);
        if mpc.z_use_ter, mpc.Cz_ter(:,:) = Cz(mpc.z_rows_ter,:,mpc.N); end
    end
    if mpc.z_use_k0, mpc.Cz_0(:,:) = mpc.Cz(mpc.z_rows_k0,:,1); end

    if mpc.z_use_ter
        mpc.grad_z_ter(:,:) = mpc.Cz_ter';
    end
end

if ~isempty(Dz)
    update_grad = 1;

    len_Dz = size(Dz,3);
    if len_Dz < mpc.N-1
        mpc.Dz(:,:,:) = fill_mat(mpc.Dz, Dz, 1);
    else
        mpc.Dz(:,:,:) = Dz(:,:,1:mpc.N-1);
    end
    % if there is Dz it means there is k0
    mpc.Dz_0(:,:) = mpc.Dz(mpc.z_rows_k0,:,1);

    mpc.grad_z_0(:,:) = mpc.Dz_0';
end

if ~isempty(Dsuz)
    update_grad = 1;
    
    len_Dsuz = size(Dsuz,3);
    if len_Dsuz < mpc.N-1
        mpc.Dsuz(:,:,:) = fill_mat(mpc.Dsuz, Dsuz, 1);
    else
        mpc.Dsuz(:,:,:) = Dsuz(:,:,1:mpc.N-1);
    end
    if mpc.z_use_k0, mpc.Dsuz_0(:,:) = mpc.Dsuz(mpc.z_rows_k0,:,1); end
end

if ~isempty(Ddz)   

    len_Ddz = size(Ddz,3);
    if len_Ddz < mpc.N-1
        mpc.Ddz(:,:,:) = fill_mat(mpc.Ddz, Ddz, 1);
    else
        mpc.Ddz(:,:,:) = Ddz(:,:,1:mpc.N-1);
    end
    if mpc.z_use_k0, mpc.Ddz_0(:,:) = mpc.Ddz(mpc.z_rows_k0,:,1); end
end

if update_grad

    for k = 1:mpc.N-1
        if mpc.z_use_s && mpc.z_use_su && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_s && mpc.z_use_su
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'];
        elseif mpc.z_use_s && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_su && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_s
            mpc.grad_z(:,:,k) = mpc.Cz(:,:,k)';
        elseif mpc.z_use_su
            mpc.grad_z(:,:,k) = mpc.Dsuz(:,:,k)';
        elseif mpc.z_use_u
            mpc.grad_z(:,:,k) = mpc.Dz(:,:,k)';
        end
>>>>>>> d5a9fce (adapt update api to 3D time staged matrices)
    end
    if ~isempty(mpc.lin_custom_cost)
        mpc.update_customcost_lin = true;
    end
    mpc.Dsuz(:,:,:) = fill_mat(mpc.Dsuz, Dsuz, 1);
    if ~isempty(mpc.z_use_k0), mpc.Dsuz_0(:,:) = Dsuz(mpc.z_rows_k0,:,1); end
end

if ~isempty(Ddz) && ~isempty(mpc.z_use_d)
    mpc.Ddz(:,:,:) = fill_mat(mpc.Ddz, Ddz, 1);
    if ~isempty(mpc.z_use_k0), mpc.Ddz_0(:,:) = Ddz(mpc.z_rows_k0,:,1); end
end

end
