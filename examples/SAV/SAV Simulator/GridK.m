function u = GridK(x,vx,rho,K)
            
    if abs(vx)<0.05
        vx = 0.05;
    end
    
    Ki = K{1} + K{2}*rho + K{3}*rho*vx + K{4}*rho*vx^2 + K{5}*rho*1/vx;
     
    u = Ki*x;

end

