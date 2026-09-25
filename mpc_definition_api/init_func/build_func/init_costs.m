function mpc = init_costs(mpc)

mpc.R_f0_0 = zeros(mpc.nu);
mpc.Q_f0_k = zeros(mpc.nse,mpc.nse,mpc.N-1);
mpc.R_f0_k = zeros(mpc.nu,mpc.nu,mpc.N-1);
mpc.Y_f0_k = zeros(mpc.nu,mpc.nse,mpc.N-1);
mpc.Q_f0_ter = zeros(mpc.nse);

mpc.R_0 = zeros(mpc.nu);
mpc.Q_k = zeros(mpc.nse,mpc.nse,mpc.N-1);
mpc.R_k = zeros(mpc.nu,mpc.nu,mpc.N-1);
mpc.Y_k = zeros(mpc.nu,mpc.nse,mpc.N-1);
mpc.Q_ter = zeros(mpc.nse);

mpc.grad_u_f0_0 = zeros(mpc.nu,1);
mpc.grad_se_f0_k = zeros(mpc.nse,mpc.N-1);
mpc.grad_u_f0_k = zeros(mpc.nu,mpc.N-1);
mpc.grad_se_f0_ter = zeros(mpc.nse,1);

if ~isempty(mpc.tracking_cost)
    mpc = genTrackingCost(mpc);   
end

if ~isempty(mpc.quad_control_cost) || ~isempty(mpc.lin_control_cost)
    mpc = genControlCost(mpc);
end

if ~isempty(mpc.controlrate_cost)
    mpc = genDiffControlCost(mpc);
end

if ~isempty(mpc.quad_custom_cost) || ~isempty(mpc.lin_custom_cost)
    mpc = genCustomCost(mpc);
end

if any(mpc.nv_k)
   mpc = genSoftSlacksCost(mpc);  
end

if ~isempty(mpc.ter_ingredients)
    mpc.P2 = 2*mpc.P;
    mpc.Q_f0_ter(mpc.s_col,mpc.s_col) = mpc.Q_f0_ter(mpc.s_col,mpc.s_col) + mpc.P2;
end

end

function mpc = genCustomCost(mpc)

s_index = mpc.s_col;
su_index = mpc.su_col;

% gradient/Hessian of quadratic cost on z
if ~isempty(mpc.quad_custom_cost)

    % k = 0
    if ~isempty(mpc.z_use_k0)
        mpc.grad_u_Z_0 = mpc.Dz_0'*mpc.Qz_0;
        mpc.R_Z_0 = mpc.Dz_0'*mpc.Qz_0*mpc.Dz_0;

        mpc.R_f0_0 = mpc.R_f0_0 + mpc.R_Z_0;
    else
        mpc.grad_u_Z_0 = [];
        mpc.R_Z_0 = [];
    end

    mpc.grad_s_Z = [];
    mpc.grad_su_Z = [];
    mpc.grad_u_Z = [];
    mpc.Q_Z = [];
    mpc.R_Z = [];
    mpc.Y_Z = [];

    if ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)
        mpc.grad_s_Z = zeros(mpc.nx,mpc.nz,mpc.N-1);
        mpc.grad_su_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.grad_u_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nse,mpc.nse,mpc.N-1);
        mpc.R_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
        mpc.Y_Z = zeros(mpc.nu,mpc.nse,mpc.N-1);
    elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su)
        mpc.grad_s_Z = zeros(mpc.nx,mpc.nz,mpc.N-1);
        mpc.grad_su_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nse,mpc.nse,mpc.N-1);
    elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_u)
        mpc.grad_s_Z = zeros(mpc.nx,mpc.nz,mpc.N-1);
        mpc.grad_u_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nx,mpc.nx,mpc.N-1);
        mpc.R_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
        mpc.Y_Z = zeros(mpc.nu,mpc.nx,mpc.N-1);
    elseif ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)
        mpc.grad_su_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.grad_u_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
        mpc.R_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
        mpc.Y_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
    elseif ~isempty(mpc.z_use_s)
        mpc.grad_s_Z = zeros(mpc.nx,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nx,mpc.nx,mpc.N-1);
    elseif ~isempty(mpc.z_use_su)
        mpc.grad_su_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.Q_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
    elseif ~isempty(mpc.z_use_u)
        mpc.grad_u_Z = zeros(mpc.nu,mpc.nz,mpc.N-1);
        mpc.R_Z = zeros(mpc.nu,mpc.nu,mpc.N-1);
    end

    for k = 1:mpc.N-1

        if ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)

            mpc.grad_s_Z(:,:,k) = mpc.Cz(:,:,k)'*mpc.Qz(:,:,k);
            mpc.grad_su_Z(:,:,k) = mpc.Dsuz(:,:,k)'*mpc.Qz(:,:,k);
            mpc.grad_u_Z(:,:,k) = mpc.Dz(:,:,k)'*mpc.Qz(:,:,k);


            mpc.Q_Z(s_index,s_index,k) = mpc.grad_s_Z(:,:,k)*mpc.Cz(:,:,k);
            mpc.Q_Z(s_index,su_index,k) = mpc.grad_s_Z(:,:,k)*mpc.Dsuz(:,:,k);
            mpc.Q_Z(su_index,s_index,k) = mpc.grad_su_Z(:,:,k)*mpc.Cz(:,:,k);
            mpc.Q_Z(su_index,su_index,k) = mpc.grad_su_Z(:,:,k)*mpc.Dsuz(:,:,k);

            mpc.Y_Z(:,s_index,k) = mpc.grad_u_Z(:,:,k)*mpc.Cz(:,:,k);
            mpc.Y_Z(:,su_index,k) = mpc.grad_u_Z(:,:,k)*mpc.Dsuz(:,:,k);

            mpc.R_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Dz(:,:,k);

        elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su)

            mpc.grad_s_Z(:,:,k) = mpc.Cz(:,:,k)'*mpc.Qz(:,:,k);
            mpc.grad_su_Z(:,:,k) = mpc.Dsuz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.Q_Z(s_index,s_index,k) = mpc.grad_s_Z(:,:,k)*mpc.Cz(:,:,k);
            mpc.Q_Z(s_index,su_index,k) = mpc.grad_s_Z(:,:,k)*mpc.Dsuz(:,:,k);
            mpc.Q_Z(su_index,s_index,k) = mpc.grad_su_Z(:,:,k)*mpc.Cz(:,:,k);
            mpc.Q_Z(su_index,su_index,k) = mpc.grad_su_Z(:,:,k)*mpc.Dsuz(:,:,k);

        elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_u)

            mpc.grad_s_Z(:,:,k) = mpc.Cz(:,:,k)'*mpc.Qz(:,:,k);
            mpc.grad_u_Z(:,:,k) = mpc.Dz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.Q_Z(:,:,k) = mpc.grad_s_Z(:,:,k)*mpc.Cz(:,:,k);

            mpc.Y_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Cz(:,:,k);

            mpc.R_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Dz(:,:,k);

        elseif ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)

            mpc.grad_su_Z(:,:,k) = mpc.Dsuz(:,:,k)'*mpc.Qz(:,:,k);
            mpc.grad_u_Z(:,:,k) = mpc.Dz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.Q_Z(:,:,k) = mpc.grad_su_Z(:,:,k)*mpc.Dsuz(:,:,k);

            mpc.Y_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Dsuz(:,:,k);

            mpc.R_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Dz(:,:,k);

        elseif ~isempty(mpc.z_use_s)

            mpc.grad_s_Z(:,:,k) = mpc.Cz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.Q_Z(:,:,k) = mpc.grad_s_Z(:,:,k)*mpc.Cz(:,:,k);

        elseif ~isempty(mpc.z_use_su)

            mpc.grad_su_Z(:,:,k) = mpc.Dsuz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.Q_Z(:,:,k) = mpc.grad_su_Z(:,:,k)*mpc.Dsuz(:,:,k);

        elseif ~isempty(mpc.z_use_u)

            mpc.grad_u_Z(:,:,k) = mpc.Dz(:,:,k)'*mpc.Qz(:,:,k);

            mpc.R_Z(:,:,k) = mpc.grad_u_Z(:,:,k)*mpc.Dz(:,:,k);

        end
    end

    if ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)
        mpc.Q_f0_k = mpc.Q_f0_k + mpc.Q_Z;
        mpc.R_f0_k = mpc.R_f0_k + mpc.R_Z;
        mpc.Y_f0_k = mpc.Y_f0_k + mpc.Y_Z;

    elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su)
        mpc.Q_f0_k = mpc.Q_f0_k + mpc.Q_Z;

    elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_u)
        mpc.Q_f0_k(s_index,s_index,:) = mpc.Q_f0_k(s_index,s_index,:) + mpc.Q_Z;
        mpc.R_f0_k = mpc.R_f0_k + mpc.R_Z;
        mpc.Y_f0_k(:,s_index,:) = mpc.Y_f0_k(:,s_index,:) + mpc.Y_Z;

    elseif ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)
        mpc.Q_f0_k(su_index,su_index,:) = mpc.Q_f0_k(su_index,su_index,:) + mpc.Q_Z;
        mpc.R_f0_k = mpc.R_f0_k + mpc.R_Z;
        mpc.Y_f0_k(:,su_index,:) = mpc.Y_f0_k(:,su_index,:) + mpc.Y_Z;

    elseif ~isempty(mpc.z_use_s)
        mpc.Q_f0_k(s_index,s_index,:) = mpc.Q_f0_k(s_index,s_index,:) + mpc.Q_Z;

    elseif ~isempty(mpc.z_use_su)
        mpc.Q_f0_k(su_index,su_index,:) = mpc.Q_f0_k(su_index,su_index,:) + mpc.Q_Z;

    elseif ~isempty(mpc.z_use_u)
        mpc.R_f0_k = mpc.R_f0_k + mpc.R_Z;
    end

    if ~isempty(mpc.z_use_ter)
        mpc.grad_s_Z_ter = mpc.Cz_ter'*mpc.Qz_ter;

        mpc.Q_Z_ter = mpc.grad_s_Z_ter*mpc.Cz_ter;
        mpc.Q_f0_ter(s_index,s_index) = mpc.Q_f0_ter(s_index,s_index) + mpc.Q_Z_ter;
    else
        mpc.grad_s_Z_ter = [];
        mpc.Q_Z_ter = [];
    end
end

% gradient of linear cost on z
if ~isempty(mpc.lin_custom_cost)

    % k = 0
    if ~isempty(mpc.z_use_k0)
        mpc.grad_u_Zlin_0 = mpc.Dz_0'*mpc.qz_0;
    else
        mpc.grad_u_Zlin_0 = [];
    end

    mpc.grad_s_Zlin = [];
    mpc.grad_su_Zlin = [];
    mpc.grad_u_Zlin = [];

    if ~isempty(mpc.z_use_s)
        mpc.grad_s_Zlin = zeros(mpc.nx,mpc.N-1);
    end
    if ~isempty(mpc.z_use_su)
        mpc.grad_su_Zlin = zeros(mpc.nu,mpc.N-1);
    end
    if ~isempty(mpc.z_use_u)
        mpc.grad_u_Zlin = zeros(mpc.nu,mpc.N-1);
    end

    for k = 1:mpc.N-1

        if ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)

            mpc.grad_s_Zlin(:,k) = mpc.Cz(:,:,k)'*mpc.qz(:,k);
            mpc.grad_su_Zlin(:,k) = mpc.Dsuz(:,:,k)'*mpc.qz(:,k);
            mpc.grad_u_Zlin(:,k) = mpc.Dz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_su)

            mpc.grad_s_Zlin(:,k) = mpc.Cz(:,:,k)'*mpc.qz(:,k);
            mpc.grad_su_Zlin(:,k) = mpc.Dsuz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_s) && ~isempty(mpc.z_use_u)

            mpc.grad_s_Zlin(:,k) = mpc.Cz(:,:,k)'*mpc.qz(:,k);
            mpc.grad_u_Zlin(:,k) = mpc.Dz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_su) && ~isempty(mpc.z_use_u)

            mpc.grad_su_Zlin(:,k) = mpc.Dsuz(:,:,k)'*mpc.qz(:,k);
            mpc.grad_u_Zlin(:,k) = mpc.Dz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_s)

            mpc.grad_s_Zlin(:,k) = mpc.Cz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_su)

            mpc.grad_su_Zlin(:,k) = mpc.Dsuz(:,:,k)'*mpc.qz(:,k);
        elseif ~isempty(mpc.z_use_u)

            mpc.grad_u_Zlin(:,k) = mpc.Dz(:,:,k)'*mpc.qz(:,k);
        end
    end

    if ~isempty(mpc.z_use_ter)
        mpc.grad_s_Zlin_ter = mpc.Cz_ter'*mpc.qz_ter;
    else
        mpc.grad_s_Zlin_ter = [];
    end
end

end

function mpc = genTrackingCost(mpc)
% gradient of err = ref-y
% gradients are stored positive, but substracted in grad_f0_MPC to account
% for negative sign
s_index = mpc.s_col;

if ~isempty(mpc.y_use_k0)
    mpc.grad_u_E_0 = mpc.D_0'*mpc.Qe_0;
    mpc.R_E_0 =  mpc.D_0'*mpc.Qe_0*mpc.D_0;

    mpc.R_f0_0 = mpc.R_f0_0 + mpc.R_E_0;
else
    mpc.grad_u_E_0 = [];
    mpc.R_E_0 = []; 
end


mpc.grad_s_E = [];
mpc.grad_u_E = [];
mpc.Q_E = [];
mpc.R_E = [];
mpc.Y_E = [];

if ~isempty(mpc.y_use_s) && ~isempty(mpc.y_use_u)
    mpc.grad_s_E = zeros(mpc.nx,mpc.ny,mpc.N-1);
    mpc.grad_u_E = zeros(mpc.nu,mpc.ny,mpc.N-1);
    mpc.Q_E = zeros(mpc.nx,mpc.nx,mpc.N-1);
    mpc.R_E = zeros(mpc.nu,mpc.nu,mpc.N-1);
    mpc.Y_E = zeros(mpc.nu,mpc.nx,mpc.N-1);
elseif ~isempty(mpc.y_use_s)
    mpc.grad_s_E = zeros(mpc.nx,mpc.ny,mpc.N-1);
    mpc.Q_E = zeros(mpc.nx,mpc.nx,mpc.N-1);
elseif ~isempty(mpc.y_use_u)
    mpc.grad_u_E = zeros(mpc.nu,mpc.ny,mpc.N-1);
    mpc.R_E = zeros(mpc.nu,mpc.nu,mpc.N-1);
end

for k = 1:mpc.N-1
    
    if ~isempty(mpc.y_use_s) && ~isempty(mpc.y_use_u)
        mpc.grad_s_E(:,:,k) = mpc.C(:,:,k)'*mpc.Qe(:,:,k);
        mpc.grad_u_E(:,:,k) = mpc.D(:,:,k)'*mpc.Qe(:,:,k);

        mpc.Q_E(:,:,k) = mpc.C(:,:,k)'*mpc.Qe(:,:,k)*mpc.C(:,:,k);
        mpc.R_E(:,:,k) = mpc.D(:,:,k)'*mpc.Qe(:,:,k)*mpc.D(:,:,k);
        mpc.Y_E(:,:,k) = mpc.D(:,:,k)'*mpc.Qe(:,:,k)*mpc.C(:,:,k);
    elseif ~isempty(mpc.y_use_s)
        mpc.grad_s_E(:,:,k) = mpc.C(:,:,k)'*mpc.Qe(:,:,k);

        mpc.Q_E(:,:,k) = mpc.C(:,:,k)'*mpc.Qe(:,:,k)*mpc.C(:,:,k);
    elseif ~isempty(mpc.y_use_u)
        mpc.grad_u_E(:,:,k) = mpc.D(:,:,k)'*mpc.Qe(:,:,k);

        mpc.R_E(:,:,k) = mpc.D(:,:,k)'*mpc.Qe(:,:,k)*mpc.D(:,:,k);
    end
end

if ~isempty(mpc.y_use_s) && ~isempty(mpc.y_use_u)
    mpc.Q_f0_k(s_index,s_index,:) = mpc.Q_f0_k(s_index,s_index,:) + mpc.Q_E;
    mpc.R_f0_k = mpc.R_f0_k + mpc.R_E;
    mpc.Y_f0_k(:,s_index,:) = mpc.Y_f0_k(:,s_index,:) + mpc.Y_E;
elseif ~isempty(mpc.y_use_s)
    mpc.Q_f0_k(s_index,s_index,:) = mpc.Q_f0_k(s_index,s_index,:) + mpc.Q_E;
elseif ~isempty(mpc.y_use_u)
    mpc.R_f0_k = mpc.R_f0_k + mpc.R_E;
end

if ~isempty(mpc.y_use_ter)
    mpc.grad_s_E_ter = mpc.C_ter'*mpc.Qe_ter;
    mpc.Q_E_ter = mpc.C_ter'*mpc.Qe_ter*mpc.C_ter;

    mpc.Q_f0_ter(s_index,s_index) = mpc.Q_f0_ter(s_index,s_index) + mpc.Q_E_ter;
else
    mpc.grad_s_E_ter = [];
    mpc.Q_E_ter = []; 
end

end

function mpc = genControlCost(mpc)

if ~isempty(mpc.quad_control_cost)

    mpc.R_f0_0 = mpc.R_f0_0 + mpc.Ru(:,:,1);

    for k = 1:mpc.N-1
        mpc.R_f0_k(:,:,k) = mpc.R_f0_k(:,:,k) + mpc.Ru(:,:,k+1);
    end
end

end

function mpc = genDiffControlCost(mpc)

su_index = mpc.su_col;

% k = 0
mpc.R_f0_0 = mpc.R_f0_0 + mpc.Rdu(:,:,1);

for k = 1:mpc.N-1

    Rdu_k = mpc.Rdu(:,:,k+1);

    mpc.Q_f0_k(su_index,su_index,k) = mpc.Q_f0_k(su_index,su_index,k) + Rdu_k;
    mpc.R_f0_k(:,:,k) = mpc.R_f0_k(:,:,k) + Rdu_k;
    mpc.Y_f0_k(:,su_index,k) = mpc.Y_f0_k(:,su_index,k) - Rdu_k;
end

end

function mpc = genSoftSlacksCost(mpc)

grad_qv_0 = zeros(mpc.nv_k(1),1);
grad_qv_k = zeros(mpc.nv_k(2),mpc.N-1);
grad_qv_ter = zeros(mpc.nv_k(3),1);

if mpc.nv_k(1)
    if ~isempty(mpc.has_y_cnstr)
        if ~isempty(mpc.y_cnstr.min_limit) && ~isempty(mpc.y_use_k0)
            row = mpc.y_cnstr.min_row_v_0;
            grad_qv_0(row) = mpc.y_cnstr.qv_min_0;
        end
        if ~isempty(mpc.y_cnstr.max_limit) && ~isempty(mpc.y_use_k0)
            row = mpc.y_cnstr.max_row_v_0;
            grad_qv_0(row) = mpc.y_cnstr.qv_max_0;
        end
    end

    if ~isempty(mpc.has_h_cnstr)
        if ~isempty(mpc.h_cnstr.min_limit) && ~isempty(mpc.h_cnstr.use_k0)
            row = mpc.h_cnstr.min_row_v_0;
            grad_qv_0(row) = mpc.h_cnstr.qv_min_0;
        end
        if ~isempty(mpc.h_cnstr.max_limit) && ~isempty(mpc.h_cnstr.use_k0)
            row = mpc.h_cnstr.max_row_v_0;
            grad_qv_0(row) = mpc.h_cnstr.qv_max_0;
        end
    end
end

if mpc.nv_k(2)
for k = 1:mpc.N-1
    if ~isempty(mpc.has_s_cnstr)
        if ~isempty(mpc.s_cnstr.min_limit)
            row = mpc.s_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.s_cnstr.qv_min(:,k);
        end
        if ~isempty(mpc.s_cnstr.max_limit)
            row = mpc.s_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.s_cnstr.qv_max(:,k);
        end
    end

    if ~isempty(mpc.has_y_cnstr)
        if ~isempty(mpc.y_cnstr.min_limit)
            row = mpc.y_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.y_cnstr.qv_min(:,k);
        end
        if ~isempty(mpc.y_cnstr.max_limit)
            row = mpc.y_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.y_cnstr.qv_max(:,k);
        end
    end

    if ~isempty(mpc.has_h_cnstr)
        if ~isempty(mpc.h_cnstr.min_limit)
            row = mpc.h_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.h_cnstr.qv_min(:,k);
        end
        if ~isempty(mpc.h_cnstr.max_limit)
            row = mpc.h_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.h_cnstr.qv_max(:,k);
        end
    end
    
end
end

% k = N

if mpc.ng_k(3)
    
if ~isempty(mpc.has_s_cnstr)
    if ~isempty(mpc.s_cnstr.min_limit)
        row = mpc.s_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_min_ter;
    end
    if ~isempty(mpc.s_cnstr.max_limit)
        row = mpc.s_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_max_ter;
    end
end

if ~isempty(mpc.has_y_cnstr)
    if ~isempty(mpc.y_cnstr.min_limit) && ~isempty(mpc.y_use_ter)
        row = mpc.y_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.y_cnstr.qv_min_ter;
    end
    if ~isempty(mpc.y_cnstr.max_limit) && ~isempty(mpc.y_use_ter)
        row = mpc.y_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.y_cnstr.qv_max_ter;
    end
end

if ~isempty(mpc.has_h_cnstr)
    if ~isempty(mpc.h_cnstr.min_limit) && ~isempty(mpc.h_cnstr.use_ter)
        row = mpc.h_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.h_cnstr.qv_min_ter;
    end
    if ~isempty(mpc.h_cnstr.max_limit) && ~isempty(mpc.h_cnstr.use_ter)
        row = mpc.h_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.h_cnstr.qv_max_ter;
    end
end

end

mpc.grad_qv_0 = grad_qv_0;
mpc.grad_qv_k = grad_qv_k;
mpc.grad_qv_ter = grad_qv_ter;

end
