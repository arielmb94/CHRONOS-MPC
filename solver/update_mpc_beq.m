function mpc = update_mpc_beq(mpc,x_prev,d)

    if isempty(d) || isempty(mpc.Bd)
        mpc.beq(1:mpc.nx) = -mpc.A*x_prev;
    else

        % copy d N times
        if length(d) < mpc.Nd
            d = repmat(d,mpc.N,1);
        end

        for k = 0:mpc.N-1        
            switch k 
                case 0
                    mpc.beq(1:mpc.nx) = -mpc.A*x_prev-mpc.Bd*d(1:mpc.nd);            
                otherwise        
                    mpc.beq(k*mpc.nx+1:(k+1)*mpc.nx) = ...
                                        -mpc.Bd*d(k*mpc.nd+1:(k+1)*mpc.nd);
            end
        end
    end

end