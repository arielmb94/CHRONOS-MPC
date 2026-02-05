function mpc = expand_gradients_hessians(mpc)

    % Expand Equality constraints
    mpc = expand_equalities(mpc.Nx+mpc.Nu,mpc.Nv,mpc);

    % Expand the data for each constraint
    if ~isempty(mpc.s_cnstr)
        mpc.s_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.s_cnstr);
    end
    if ~isempty(mpc.s_ter_cnstr)
        mpc.s_ter_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.s_ter_cnstr);
    end
    if ~isempty(mpc.u_cnstr)
        mpc.u_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.u_cnstr);
    end
    if ~isempty(mpc.du_cnstr)
        mpc.du_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.du_cnstr);
    end
    if ~isempty(mpc.y_cnstr)
        mpc.y_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.y_cnstr);
    end
    if ~isempty(mpc.h_cnstr)
        mpc.h_cnstr = expand_constraint(mpc.Nx+mpc.Nu,mpc.Nv,mpc.h_cnstr);
    end

    % Expand Tracking Cost
    if ~isempty(mpc.Qe)
        mpc = expand_tracking_cost(mpc.Nx+mpc.Nu,mpc.Nv,mpc);
    end
    % Expand Control Cost
    if ~isempty(mpc.Ru) ||  ~isempty(mpc.ru)
        mpc = expand_control_cost(mpc.Nx+mpc.Nu,mpc.Nv,mpc);
    end
    % Expand Differential Control Cost
    if ~isempty(mpc.Rdu)
        mpc = expand_diff_control_cost(mpc.Nx+mpc.Nu,mpc.Nv,mpc);
    end
    % Expand Performance Cost
    if ~isempty(mpc.Qz) ||  ~isempty(mpc.qz)
        mpc = expand_performance_cost(mpc.Nx+mpc.Nu,mpc.Nv,mpc);
    end

    % Recompute the hessian of the cost function
    if ~isempty(mpc.hessCost)

        n = size(mpc.hessCost,1);
        length_diff = mpc.Nx+mpc.Nu+mpc.Nv - n;
        mpc.hessCost = [mpc.hessCost zeros(n,length_diff);
                            zeros(length_diff,n) zeros(length_diff)];

        mpc = update_mpc_f0_hess(mpc);
    end

end

function mpc = expand_equalities(N,Nv,mpc)

    length_diff = N+Nv - size(mpc.Aeq,2);

    if length_diff

        mpc.Aeq = [mpc.Aeq zeros(size(mpc.Aeq,1),length_diff)];
    end

end

function cnstr = expand_constraint(N,Nv,cnstr)

% gradients and Hessian min/max
if size(cnstr.grad_min,1) < N+Nv
    
    length_diff = N+Nv - size(cnstr.grad_min,1);
    cnstr.grad_min = [cnstr.grad_min;
                      zeros(length_diff,size(cnstr.grad_min,2))];


    cnstr.hess_min = genHessIneq(cnstr.grad_min);

end
if size(cnstr.grad_max,1) < N+Nv
    
    length_diff = N+Nv - size(cnstr.grad_max,1);
    cnstr.grad_max = [cnstr.grad_max;
                      zeros(length_diff,size(cnstr.grad_max,2))];

    cnstr.hess_max = genHessIneq(cnstr.grad_max);
end

% adapt slack variable constraints if existent
if cnstr.min_slack_nv
    
    length_diff = N+Nv - size(cnstr.min_slack_index,2);
    % adapt positivity constraint
    if length_diff

        cnstr.min_slack_index = [cnstr.min_slack_index, zeros(size(cnstr.min_slack_index,1),length_diff)];

        cnstr.min_slack_positivity_grad = [cnstr.min_slack_positivity_grad;
                                           zeros(length_diff,size(cnstr.min_slack_positivity_grad,2))];

        cnstr.min_slack_positivity_hess = genHessIneq(cnstr.min_slack_positivity_grad);
    end
    % adapt hard limit constraint
    if cnstr.min_slack_hard_limit && length_diff

        cnstr.min_slack_hard_limit_grad = [cnstr.min_slack_hard_limit_grad;
                                           zeros(length_diff,size(cnstr.min_slack_hard_limit_grad,2))];

        cnstr.min_slack_hard_limit_hess = genHessIneq(cnstr.min_slack_hard_limit_grad);
    end
end
if cnstr.max_slack_nv
    
    length_diff = N+Nv - size(cnstr.max_slack_index,2);
    % adapt positivity constraint
    if length_diff

        cnstr.max_slack_index = [cnstr.max_slack_index, zeros(size(cnstr.max_slack_index,1),length_diff)];

        cnstr.max_slack_positivity_grad = [cnstr.max_slack_positivity_grad;
                                           zeros(length_diff,size(cnstr.max_slack_positivity_grad,2))];

        cnstr.max_slack_positivity_hess = genHessIneq(cnstr.max_slack_positivity_grad);
    end
    % adapt hard limit constraint
    if cnstr.max_slack_hard_limit && length_diff

        cnstr.max_slack_hard_limit_grad = [cnstr.max_slack_hard_limit_grad;
                                           zeros(length_diff,size(cnstr.max_slack_hard_limit_grad,2))];

        cnstr.max_slack_hard_limit_hess = genHessIneq(cnstr.max_slack_hard_limit_grad);
    end
end

end

function mpc = expand_tracking_cost(N,Nv,mpc)
    length_diff = N+Nv - size(mpc.gradErrQe,1);

    if length_diff

        n = size(mpc.gradErrQe,1);
        m = size(mpc.gradErrQe,2);

        mpc.gradErrQe = [mpc.gradErrQe;
                         zeros(length_diff,m)];

        mpc.hessErrTerm = [mpc.hessErrTerm zeros(n,length_diff);
                           zeros(length_diff,n) zeros(length_diff)];
    end

end

function mpc = expand_control_cost(N,Nv,mpc)

if ~isempty(mpc.Ru)
    length_diff = N+Nv - size(mpc.gradCtlrRu,1);

    if length_diff

        n = size(mpc.gradCtlrRu,1);
        m = size(mpc.gradCtlrRu,2);

        mpc.gradCtlrRu = [mpc.gradCtlrRu;
            zeros(length_diff,m)];

        mpc.hessCtrlTerm = [mpc.hessCtrlTerm zeros(n,length_diff);
            zeros(length_diff,n) zeros(length_diff)];
    end
end

if ~isempty(mpc.ru)
    length_diff = N+Nv - size(mpc.gradCtlrru,1);

    if length_diff

        m = size(mpc.gradCtlrru,2);

        mpc.gradCtlrru = [mpc.gradCtlrru;
            zeros(length_diff,m)];
    end
end
    
end

function mpc = expand_diff_control_cost(N,Nv,mpc)
    length_diff = N+Nv - size(mpc.gradDiffCtlrR,1);

    if length_diff

        n = size(mpc.gradDiffCtlrR,1);
        m = size(mpc.gradDiffCtlrR,2);

        mpc.gradDiffCtlrR = [mpc.gradDiffCtlrR;
                             zeros(length_diff,m)];

        mpc.gradDiffCtlr = [mpc.gradDiffCtlr;
                            zeros(length_diff,m)];

        mpc.hessDiffCtrlTerm = [mpc.hessDiffCtrlTerm zeros(n,length_diff);
                                zeros(length_diff,n) zeros(length_diff)];
    end

end

function mpc = expand_performance_cost(N,Nv,mpc)

if ~isempty(mpc.Qz)
    length_diff = N+Nv - size(mpc.gradPerfQz,1);

    if length_diff

        n = size(mpc.gradPerfQz,1);
        m = size(mpc.gradPerfQz,2);

        mpc.gradPerfQz = [mpc.gradPerfQz;
                          zeros(length_diff,m)];

        mpc.hessPerfTerm = [mpc.hessPerfTerm zeros(n,length_diff);
                            zeros(length_diff,n) zeros(length_diff)];
    end
end

if ~isempty(mpc.qz)
    length_diff = N+Nv - size(mpc.gradPerfqz,1);

    if length_diff

        m = size(mpc.gradPerfqz,2);

        mpc.gradPerfqz = [mpc.gradPerfqz;
                          zeros(length_diff,m)];
    end
end
    
end
