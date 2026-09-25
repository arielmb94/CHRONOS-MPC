function mpc = riccati_KKT(mpc)
[mpc.delta_u,mpc.delta_se,mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,mpc.rp_hat_0,mpc.rp_hat,mpc.K_ric,mpc.d_ric,mpc.d_ric_0,mpc.QA_ric,mpc.QB_ric,mpc.solve_rhs_ric,mpc.solve_result_ric] = riccati_KKT_local(mpc.N,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,mpc.rp_0,mpc.rp_k,mpc.rp_hat_0,mpc.rp_hat,mpc.K_ric,mpc.d_ric,mpc.d_ric_0,mpc.QA_ric,mpc.QB_ric,mpc.solve_rhs_ric,mpc.solve_result_ric,mpc.A,mpc.B,mpc.delta_u,mpc.delta_se,mpc.eps_thknv,mpc.nu,mpc.nse,mpc.se_col);
end

function [delta_u,delta_se,Q_hat,R_hat_0,R_hat,Y_hat,rs_hat,ru_hat,ru_hat_0,rp_hat_0,rp_hat,K_ric,d_ric,d_ric_0,QA_ric,QB_ric,solve_rhs_ric,solve_result_ric] = riccati_KKT_local(N,R_0,Q,R,Y,Q_ter,Q_hat,R_hat_0,R_hat,Y_hat,ru_0,rs,ru,rs_ter,rs_hat,ru_hat,ru_hat_0,rp_0,rp,rp_hat_0,rp_hat,K_ric,d_ric,d_ric_0,QA_ric,QB_ric,solve_rhs_ric,solve_result_ric,A,B,delta_u,delta_se,eps_thknv,nu,nse,se_col)

linsolve_opts = struct('SYM',true);

Q_hat(:,:,N) = Q_ter;
rs_hat(:,N) = rs_ter;

for k = N-1:-1:1

    % Cache Q_hat*A and Q_hat*B once for this stage.
    QA_ric(:,:) = Q_hat(:,:,k+1)*A(:,:,k+1);
    QB_ric(:,:) = Q_hat(:,:,k+1)*B(:,:,k+1);

    % \hat{R}*k = R_k + B^T\hat{Q}_{k+1}B
    R_hat(:,:,k) = R(:,:,k) + B(:,:,k+1)'*QB_ric;

    for j = 1:nu
        R_hat(j,j,k) = R_hat(j,j,k) + eps_thknv;
    end

    % \hat{Y}_k = Y_k + B^T\hat{Q}_{k+1}A
    Y_hat(:,:,k) = Y(:,:,k) + B(:,:,k+1)'*QA_ric;

    % \hat{r}_{pk} = \hat{r}_{sk+1} + \hat{Q}_{k+1}r_{pk}
    rp_hat(:,k) = rs_hat(:,k+1) + ...
                      Q_hat(:,:,k+1)*rp(:,k);

    % \hat{r}_{uk} = r_{uk} + B^T\hat{r}_{pk}
    ru_hat(:,k) = ru(:,k) + B(:,:,k+1)'*rp_hat(:,k);

    % Solve [K_k, d_k] in one factorization of the regularized block.
    solve_rhs_ric(:,se_col) = Y_hat(:,:,k);
    solve_rhs_ric(:,nse+1) = ru_hat(:,k);

    solve_result_ric(:,:) = linsolve(R_hat(:,:,k),...
                                     solve_rhs_ric,linsolve_opts);

    K_ric(:,:,k) = solve_result_ric(:,se_col);
    d_ric(:,k) = solve_result_ric(:,nse+1);

    % \hat{Q}_k = Q_k + A^T\hat{Q}_{k+1}A - \hat{Y}_k^TK_k
    Q_hat(:,:,k) = Q(:,:,k) + A(:,:,k+1)'*QA_ric - ...
                       Y_hat(:,:,k)'*K_ric(:,:,k);

    % \hat{r}_{sk} = r_{sk} + A^T\hat{r}_{pk} - \hat{Y}_k^Td_k
    rs_hat(:,k) = rs(:,k) + A(:,:,k+1)'*rp_hat(:,k) - ...
                      Y_hat(:,:,k)'*d_ric(:,k);

end

QB_ric(:,:) = Q_hat(:,:,1)*B(:,:,1);

R_hat_0(:,:) = R_0+B(:,:,1)'*QB_ric;

for j = 1:nu
    R_hat_0(j,j) = R_hat_0(j,j) + eps_thknv;
end

rp_hat_0(:) = rs_hat(:,1) + ...
                  Q_hat(:,:,1)*rp_0;

ru_hat_0(:) = ru_0 + B(:,:,1)'*rp_hat_0;

d_ric_0(:) = linsolve(R_hat_0,ru_hat_0,...
                      linsolve_opts);

% u_0 = -d_0
delta_u(:,1) = -d_ric_0;

% s_1 = B u_0 + r_{p0}
delta_se(:,1) = B(:,:,1)*delta_u(:,1) + rp_0;

% \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
%mu(:,1) = rs_hat_ric(:,1) + Q_hat_ric(:,:,1)*delta_se(:,1);

for k = 1:N-1

    ku = k+1; % u mu vectors shifted by 1 stage

    % u_k = -d_k - K_ks_k
    delta_u(:,ku) = -d_ric(:,k) - ...
                    K_ric(:,:,k)*delta_se(:,k);

    % s_{k+1} = As_k + Bu_k + r_{pk}
    delta_se(:,k+1) = A(:,:,k+1)*delta_se(:,k) + ...
                      B(:,:,k+1)*delta_u(:,ku) + rp(:,k);

    % \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
    %mu(:,ku) = rs_hat_ric(:,k+1) + ...
    %           Q_hat_ric(:,:,k+1)*delta_se(:,k+1);

end

end

