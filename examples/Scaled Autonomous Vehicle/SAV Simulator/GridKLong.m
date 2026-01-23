function u = GridKLong(x,Vx,K)
            
    if abs(Vx)<0.05
        Vx = 0.05;
    end
    
    Ki = K{1}+Vx*K{2}+(1/Vx)*K{3};
        
    u = Ki*x;

end

