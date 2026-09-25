% INIT_MPC_CONTROL_COST Add quadratic and/or linear control penalties.
%
%   mpc = INIT_MPC_CONTROL_COST(mpc, Ru) adds the quadratic term
%
%       J_control += 0.5 * u_k' * Ru_k * u_k
%
%   mpc = INIT_MPC_CONTROL_COST(mpc, [], ru) adds the linear term
%
%       J_control += ru_k' * u_k
%
%   mpc = INIT_MPC_CONTROL_COST(mpc, Ru, ru) adds both terms. Leave either
%   weight empty when that term is not needed.
%
%   Ru may contain one matrix or L horizon stages, and ru may contain one
%   column or L horizon columns. Here, L is the number of supplied horizon
%   stages or columns. If L < mpc.N, CHRONOS reuses the last supplied value
%   for the remaining stages; if L >= mpc.N, only the first mpc.N values
%   are used.
%
%   Call this function after INIT_MPC_DYNAMICS and before
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     Ru      - Optional quadratic weight, size nu-by-nu or nu-by-nu-by-L.
%               Each active stage must be symmetric positive semidefinite.
%     ru      - Optional linear weight, size nu-by-1 or nu-by-L. Use only
%               the control action is strictly positive.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - penalize both control effort and a signed control bias:
%
%       Ru = 0.1 * eye(nu);
%       ru = zeros(nu, 1);
%       ru(1) = 0.02;
%       mpc = init_mpc_Control_cost(mpc, Ru, ru);
function mpc = init_mpc_Control_cost(mpc,Ru,ru)
arguments
    mpc
    Ru = []
    ru = [];
end

if any(Ru(:))

    validate_matrix(Ru, mpc.nu, mpc.nu, 'Ru', true);

    if isscalar(Ru)
        Ru = Ru * eye(mpc.nu);
    elseif size(Ru,1) == 1 && size(Ru,2) == 1
        Ru_staged = Ru;
        Ru = zeros(mpc.nu, mpc.nu, size(Ru_staged,3));
        Ru_eye = eye(mpc.nu);
        for k = 1:size(Ru_staged,3)
            Ru(:,:,k) = Ru_staged(1,1,k) * Ru_eye;
        end
    end

    mpc.quad_control_cost = 1;

    mpc.Ru = zeros(mpc.nu,mpc.nu,mpc.N);
    len_Ru = size(Ru,3);
    if len_Ru < mpc.N
        mpc.Ru = fill_mat(mpc.Ru, Ru, 1);
    else
        mpc.Ru = Ru(:,:,1:mpc.N);
    end
end

if any(ru(:))

    validate_column_vector(ru, mpc.nu, 'ru', true);

    if ~isempty(ru) && size(ru,1) == 1
        ru = ones(mpc.nu, 1) * ru;
    end

    mpc.lin_control_cost = 1;

    mpc.ru = zeros(mpc.nu,mpc.N);
    len_ru = size(ru,2);
    if len_ru < mpc.N
        mpc.ru = fill_vec(mpc.ru, ru, 1);
    else
        mpc.ru = ru(:,1:mpc.N);
    end    
end

end
