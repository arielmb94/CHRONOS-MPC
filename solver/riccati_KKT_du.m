function mpc = riccati_KKT_du(mpc)
[mpc.delta_u,mpc.delta_xe,mpc.Q_hat,mpc.rx_hat,mpc.rp_hat_0,mpc.rp_hat, ...
 mpc.QA_ric,mpc.G_ric,mpc.LDL_ric,mpc.solve_rhs_ric,mpc.scaled_rhs_ric] = ...
    riccati_KKT_du_local(mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter, ...
        mpc.ru_0,mpc.rxe_k,mpc.ru_k,mpc.rxe_ter,mpc.rp_0,mpc.rp_k, ...
        mpc.A,mpc.B,mpc.eps_thknv,mpc.x_col,mpc.xu_col,mpc.xe_col, ...
        mpc.delta_u,mpc.delta_xe,mpc.Q_hat,mpc.rx_hat,mpc.rp_hat_0,mpc.rp_hat, ...
        mpc.QA_ric,mpc.G_ric,mpc.LDL_ric,mpc.solve_rhs_ric,mpc.scaled_rhs_ric);
end

function [delta_u,delta_xe,Q_hat,rx_hat,rp_hat_0,rp_hat, ...
          QA_ric,G_ric,LDL_ric,solve_rhs_ric,scaled_rhs_ric] = ...
    riccati_KKT_du_local(R_0,Q,R,Y,Q_ter,ru_0,rx,ru,rx_ter,rp_0,rp, ...
        A,B,eps_thknv,x_col,xu_col,xe_col,delta_u,delta_xe,Q_hat,rx_hat, ...
        rp_hat_0,rp_hat,QA_ric,G_ric,LDL_ric,solve_rhs_ric,scaled_rhs_ric)

N = size(delta_u,2);
nu = size(delta_u,1);
residual_col = size(solve_rhs_ric,2);

Q_hat(:,:,N) = Q_ter;
rx_hat(:,N) = rx_ter;

for k = N-1:-1:1
    ku = k+1;

    % Cache the physical-state columns of [B' I]*Q_hat.
    G_ric(:,:) = B(:,:,ku)'*Q_hat(x_col,x_col,ku) + Q_hat(xu_col,x_col,ku);

    LDL_ric(:,:,ku) = R(:,:,k) + G_ric*B(:,:,ku) + ...
        B(:,:,ku)'*Q_hat(x_col,xu_col,ku) + Q_hat(xu_col,xu_col,ku);
    for i = 1:nu
        LDL_ric(i,i,ku) = LDL_ric(i,i,ku) + eps_thknv;
    end

    rp_hat(:,k) = rx_hat(:,ku) + Q_hat(:,:,ku)*rp(:,k);
    solve_rhs_ric(:,x_col) = Y(:,x_col,k) + G_ric*A(:,:,ku);
    solve_rhs_ric(:,xu_col) = Y(:,xu_col,k);
    solve_rhs_ric(:,residual_col) = ru(:,k) + B(:,:,ku)'*rp_hat(x_col,k) + ...
        rp_hat(xu_col,k);

    % Rbar*[K,d] = [Y_hat,ru_hat]; retain T = L^(-1)*RHS and Z = D^(-1)*T.
    [LDL_ric,solve_rhs_ric] = factorizeForwardRiccatiLDL( ...
        LDL_ric,solve_rhs_ric,ku,false);
    scaled_rhs_ric = scaleRiccatiLDLRHS( ...
        scaled_rhs_ric,solve_rhs_ric,LDL_ric,k);

    % Keep the augmented dynamics in their physical blocks.
    QA_ric(:,:) = Q_hat(x_col,x_col,ku)*A(:,:,ku);
    Q_hat(:,:,k) = Q(:,:,k);
    Q_hat(x_col,x_col,k) = Q_hat(x_col,x_col,k) + A(:,:,ku)'*QA_ric;
    % T_Y'*Z_Y = Y_hat'*K, without backward substitution of the state RHS.
    Q_hat(:,:,k) = Q_hat(:,:,k) - ...
        solve_rhs_ric(:,xe_col)'*scaled_rhs_ric(:,xe_col,k);

    rx_hat(:,k) = rx(:,k);
    rx_hat(x_col,k) = rx_hat(x_col,k) + A(:,:,ku)'*rp_hat(x_col,k);
    rx_hat(:,k) = rx_hat(:,k) - ...
        solve_rhs_ric(:,xe_col)'*scaled_rhs_ric(:,residual_col,k);
end

G_ric(:,:) = B(:,:,1)'*Q_hat(x_col,x_col,1) + Q_hat(xu_col,x_col,1);
LDL_ric(:,:,1) = R_0 + G_ric*B(:,:,1) + ...
    B(:,:,1)'*Q_hat(x_col,xu_col,1) + Q_hat(xu_col,xu_col,1);
for i = 1:nu
    LDL_ric(i,i,1) = LDL_ric(i,i,1) + eps_thknv;
end

rp_hat_0(:) = rx_hat(:,1) + Q_hat(:,:,1)*rp_0;
solve_rhs_ric(:,residual_col) = ru_0 + B(:,:,1)'*rp_hat_0(x_col) + rp_hat_0(xu_col);

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

% xe_1 = [B;I]*u_0 + r_p0.
delta_xe(x_col,1) = B(:,:,1)*delta_u(:,1) + rp_0(x_col);
delta_xe(xu_col,1) = delta_u(:,1) + rp_0(xu_col);

for k = 1:N-1
    ku = k+1;

    % Recover -d_k-K_k*delta_xe with one L' solve after the matrix product.
    delta_u(:,ku) = -scaled_rhs_ric(:,residual_col,k) - ...
        scaled_rhs_ric(:,xe_col,k)*delta_xe(:,k);
    delta_u = backSubstituteRiccatiLDL(delta_u,LDL_ric,ku);

    % xe_next = [A 0;0 0]*xe + [B;I]*u + r_p.
    delta_xe(x_col,ku) = A(:,:,ku)*delta_xe(x_col,k) + ...
        B(:,:,ku)*delta_u(:,ku) + rp(x_col,k);
    delta_xe(xu_col,ku) = delta_u(:,ku) + rp(xu_col,k);
end
end
