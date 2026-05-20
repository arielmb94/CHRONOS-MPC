function mpc = genEqualities(mpc,A,B,N,nx,nu)

du = mpc.has_du;

for k = 0:N-1

    if k == 0
        rows = 1:(nx+du*nu);
        columns = 1:nu+nx+du*nu;
        mpc.Aeq(rows,columns) = [B -eye(nx) zeros(nx,du*nu);
                                eye(du*nu) zeros(nu,du*nx) -eye(du*nu)];
    
    elseif k < N-1
        rows = k*(nx+du*nu)+1:(k+1)*(nx+du*nu);
        columns = (nu+(nx+du*nu))*k+1-(nx+du*nu):(nu+(nx+du*nu))*(k+1);
        mpc.Aeq(rows,columns) =...
            [A zeros(nx,du*nu) B -eye(nx) zeros(nx,du*nu);
             zeros(nu,du*(nx+nu)) eye(du*nu) zeros(nu,du*(nx)) -eye(du*nu)];

    else
        rows = k*(nx+du*nu)+1:(k+1)*(nx+du*nu)-du*nu;
        columns = (nu+(nx+du*nu))*k+1-(nx+du*nu):(nu+(nx+du*nu))*(k+1)-du*nu;
        mpc.Aeq(rows,columns) = [A zeros(nx,du*nu) B -eye(nx)];
    end
end

mpc.beq = zeros(size(mpc.Aeq,1),1);

end