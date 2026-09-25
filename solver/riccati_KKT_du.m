function mpc = riccati_KKT_du(mpc)
[mpc.delta_u,mpc.delta_se,mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,mpc.rp_hat_0,mpc.rp_hat,mpc.K_ric,mpc.d_ric,mpc.d_ric_0,mpc.QA_ric,mpc.G_ric,mpc.solve_rhs_ric,mpc.solve_result_ric] = riccati_KKT_du_local(mpc.N,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,mpc.rp_0,mpc.rp_k,mpc.rp_hat_0,mpc.rp_hat,mpc.K_ric,mpc.d_ric,mpc.d_ric_0,mpc.QA_ric,mpc.G_ric,mpc.solve_rhs_ric,mpc.solve_result_ric,mpc.A,mpc.B,mpc.delta_u,mpc.delta_se,mpc.eps_thknv,mpc.nu,mpc.nse,mpc.s_col,mpc.su_col,mpc.se_col);
end

function [delta_u,delta_se,Q_hat,R_hat_0,R_hat,Y_hat,rs_hat,ru_hat,ru_hat_0,rp_hat_0,rp_hat,K_ric,d_ric,d_ric_0,QA_ric,G_ric,solve_rhs_ric,solve_result_ric] = riccati_KKT_du_local(N,R_0,Q,R,Y,Q_ter,Q_hat,R_hat_0,R_hat,Y_hat,ru_0,rs,ru,rs_ter,rs_hat,ru_hat,ru_hat_0,rp_0,rp,rp_hat_0,rp_hat,K_ric,d_ric,d_ric_0,QA_ric,G_ric,solve_rhs_ric,solve_result_ric,A,B,delta_u,delta_se,eps_thknv,nu,nse,s_col,su_col,se_col)

linsolve_opts = struct('SYM',true);

Q_hat(:,:,N) = Q_ter;
rs_hat(:,N) = rs_ter;

for k = N-1:-1:1

    % Cache the state-column block of [B';I]*Q_hat once for this stage.
    G_ric(:,:) = B(:,:,k+1)'*Q_hat(s_col,s_col,k+1) + ...
                 Q_hat(su_col,s_col,k+1);

    % \hat{R}*k = R_k + B_kkt^T\hat{Q}_{k+1}B_kkt
    R_hat(:,:,k) = R(:,:,k) + G_ric*B(:,:,k+1) + ...
                       B(:,:,k+1)'*Q_hat(s_col,su_col,k+1) + ...
                       Q_hat(su_col,su_col,k+1);

    for j = 1:nu
        R_hat(j,j,k) = R_hat(j,j,k) + eps_thknv;
    end

    % \hat{Y}_k = Y_k + B_kkt^T\hat{Q}_{k+1}A_kkt
    Y_hat(:,s_col,k) = Y(:,s_col,k) + G_ric*A(:,:,k+1);
    Y_hat(:,su_col,k) = Y(:,su_col,k);

    % \hat{r}_{pk} = \hat{r}_{sk+1} + \hat{Q}_{k+1}r_{pk}
    rp_hat(:,k) = rs_hat(:,k+1) + ...
                      Q_hat(:,:,k+1)*rp(:,k);

    % \hat{r}_{uk} = r_{uk} + B_kkt^T\hat{r}_{pk}
    ru_hat(:,k) = ru(:,k) + B(:,:,k+1)'*rp_hat(s_col,k) + ...
                      rp_hat(su_col,k);

    % Solve [K_k, d_k] in one factorization of the regularized block.
    solve_rhs_ric(:,se_col) = Y_hat(:,:,k);
    solve_rhs_ric(:,nse+1) = ru_hat(:,k);

    solve_result_ric(:,:) = linsolve(R_hat(:,:,k),...
                                     solve_rhs_ric,linsolve_opts);

    K_ric(:,:,k) = solve_result_ric(:,se_col);
    d_ric(:,k) = solve_result_ric(:,nse+1);

    % \hat{Q}_k = Q_k + A_kkt^T\hat{Q}_{k+1}A_kkt -
    %              \hat{Y}_k^TK_k
    QA_ric(:,:) = Q_hat(s_col,s_col,k+1)*A(:,:,k+1);

    Q_hat(:,:,k) = Q(:,:,k);
    Q_hat(s_col,s_col,k) = Q_hat(s_col,s_col,k) + ...
                               A(:,:,k+1)'*QA_ric;
    Q_hat(:,:,k) = Q_hat(:,:,k) - Y_hat(:,:,k)'*K_ric(:,:,k);

    % \hat{r}_{sk} = r_{sk} + A_kkt^T\hat{r}_{pk} - \hat{Y}_k^Td_k
    rs_hat(:,k) = rs(:,k);
    rs_hat(s_col,k) = rs_hat(s_col,k) + ...
                          A(:,:,k+1)'*rp_hat(s_col,k);
    rs_hat(:,k) = rs_hat(:,k) - Y_hat(:,:,k)'*d_ric(:,k);

end

G_ric(:,:) = B(:,:,1)'*Q_hat(s_col,s_col,1) + ...
             Q_hat(su_col,s_col,1);

R_hat_0(:,:) = R_0 + G_ric*B(:,:,1) + ...
                   B(:,:,1)'*Q_hat(s_col,su_col,1) + ...
                   Q_hat(su_col,su_col,1);

for j = 1:nu
    R_hat_0(j,j) = R_hat_0(j,j) + eps_thknv;
end

rp_hat_0(:) = rs_hat(:,1) + ...
                  Q_hat(:,:,1)*rp_0;

ru_hat_0(:) = ru_0 + B(:,:,1)'*rp_hat_0(s_col) + ...
                  rp_hat_0(su_col);

d_ric_0(:) = linsolve(R_hat_0,ru_hat_0,...
                      linsolve_opts);

% u_0 = -d_0
delta_u(:,1) = -d_ric_0;

% se_1 = [B;I]u_0 + r_{p0}
delta_se(s_col,1) = B(:,:,1)*delta_u(:,1) + rp_0(s_col);
delta_se(su_col,1) = delta_u(:,1) + rp_0(su_col);

% \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
%mu(:,1) = rs_hat_ric(:,1) + Q_hat_ric(:,:,1)*delta_se(:,1);

for k = 1:N-1

    ku = k+1; % u mu vectors shifted by 1 stage

    % u_k = -d_k - K_k*se_k
    delta_u(:,ku) = -d_ric(:,k) - ...
                    K_ric(:,:,k)*delta_se(:,k);

    % se_{k+1} = [A 0;0 0]se_k + [B;I]u_k + r_{pk}
    delta_se(s_col,k+1) = A(:,:,k+1)*delta_se(s_col,k) + ...
                              B(:,:,k+1)*delta_u(:,ku) + rp(s_col,k);
    delta_se(su_col,k+1) = delta_u(:,ku) + rp(su_col,k);

    % \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
    %mu(:,ku) = rs_hat_ric(:,k+1) + ...
    %           Q_hat_ric(:,:,k+1)*delta_se(:,k+1);

end


end
