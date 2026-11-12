function mpc = riccati_KKT(mpc)
[mpc.delta_u,mpc.delta_xe,mpc.Q_hat,mpc.rx_hat,mpc.rp_hat_0,mpc.rp_hat, ...
 mpc.QA_ric,mpc.QB_ric,mpc.LDL_ric,mpc.solve_rhs_ric,mpc.scaled_rhs_ric] = ...
    riccati_KKT_local(mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter, ...
        mpc.ru_0,mpc.rxe_k,mpc.ru_k,mpc.rxe_ter,mpc.rp_0,mpc.rp_k, ...
        mpc.A,mpc.B,mpc.eps_thknv,mpc.xe_col, ...
        mpc.delta_u,mpc.delta_xe,mpc.Q_hat,mpc.rx_hat,mpc.rp_hat_0,mpc.rp_hat, ...
        mpc.QA_ric,mpc.QB_ric,mpc.LDL_ric,mpc.solve_rhs_ric,mpc.scaled_rhs_ric);
end

function [delta_u,delta_xe,Q_hat,rx_hat,rp_hat_0,rp_hat, ...
          QA_ric,QB_ric,LDL_ric,solve_rhs_ric,scaled_rhs_ric] = ...
    riccati_KKT_local(R_0,Q,R,Y,Q_ter,ru_0,rx,ru,rx_ter,rp_0,rp, ...
        A,B,eps_thknv,xe_col,delta_u,delta_xe,Q_hat,rx_hat, ...
        rp_hat_0,rp_hat,QA_ric,QB_ric,LDL_ric,solve_rhs_ric,scaled_rhs_ric)

N = size(delta_u,2);
nu = size(delta_u,1);
residual_col = size(solve_rhs_ric,2);

Q_hat(:,:,N) = Q_ter;
rx_hat(:,N) = rx_ter;

for k = N-1:-1:1
    ku = k+1;

    % Cache Q_hat*A and Q_hat*B once for this stage.
    QA_ric(:,:) = Q_hat(:,:,ku)*A(:,:,ku);
    QB_ric(:,:) = Q_hat(:,:,ku)*B(:,:,ku);

    LDL_ric(:,:,ku) = R(:,:,k) + B(:,:,ku)'*QB_ric;
    for i = 1:nu
        LDL_ric(i,i,ku) = LDL_ric(i,i,ku) + eps_thknv;
    end

    rp_hat(:,k) = rx_hat(:,ku) + Q_hat(:,:,ku)*rp(:,k);
    solve_rhs_ric(:,xe_col) = Y(:,:,k) + B(:,:,ku)'*QA_ric;
    solve_rhs_ric(:,residual_col) = ru(:,k) + B(:,:,ku)'*rp_hat(:,k);

    % Rbar*[K,d] = [Y_hat,ru_hat]; retain T = L^(-1)*RHS and Z = D^(-1)*T.
    [LDL_ric,solve_rhs_ric] = factorizeForwardRiccatiLDL( ...
        LDL_ric,solve_rhs_ric,ku,false);
    scaled_rhs_ric = scaleRiccatiLDLRHS( ...
        scaled_rhs_ric,solve_rhs_ric,LDL_ric,k);

    % T_Y'*Z_Y = Y_hat'*K, without backward substitution of the state RHS.
    Q_hat(:,:,k) = Q(:,:,k) + A(:,:,ku)'*QA_ric - ...
        solve_rhs_ric(:,xe_col)'*scaled_rhs_ric(:,xe_col,k);
    rx_hat(:,k) = rx(:,k) + A(:,:,ku)'*rp_hat(:,k) - ...
        solve_rhs_ric(:,xe_col)'*scaled_rhs_ric(:,residual_col,k);
end

QB_ric(:,:) = Q_hat(:,:,1)*B(:,:,1);
LDL_ric(:,:,1) = R_0 + B(:,:,1)'*QB_ric;
for i = 1:nu
    LDL_ric(i,i,1) = LDL_ric(i,i,1) + eps_thknv;
end

rp_hat_0(:) = rx_hat(:,1) + Q_hat(:,:,1)*rp_0;
solve_rhs_ric(:,residual_col) = ru_0 + B(:,:,1)'*rp_hat_0;

if nu == 1
    % One scalar RHS needs only a division, not reciprocal scaling.
    delta_u(1,1) = -solve_rhs_ric(1,residual_col)/LDL_ric(1,1,1);
else
    [LDL_ric,solve_rhs_ric] = factorizeForwardRiccatiLDL( ...
        LDL_ric,solve_rhs_ric,1,true);
    for i = 1:nu
        delta_u(i,1) = -LDL_ric(i,i,1)*solve_rhs_ric(i,residual_col);
    end
    delta_u = backSubstituteRiccatiLDL(delta_u,LDL_ric,1);
end

% x_1 = B*u_0 + r_p0.
delta_xe(:,1) = B(:,:,1)*delta_u(:,1) + rp_0;

for k = 1:N-1
    ku = k+1;

    % Recover -d_k-K_k*delta_xe with one L' solve after the matrix product.
    delta_u(:,ku) = -scaled_rhs_ric(:,residual_col,k) - ...
        scaled_rhs_ric(:,xe_col,k)*delta_xe(:,k);
    delta_u = backSubstituteRiccatiLDL(delta_u,LDL_ric,ku);

    delta_xe(:,ku) = A(:,:,ku)*delta_xe(:,k) + ...
        B(:,:,ku)*delta_u(:,ku) + rp(:,k);
end
end
