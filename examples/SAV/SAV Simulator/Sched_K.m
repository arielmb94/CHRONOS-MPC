function u = Sched_K(x,vx,K,Problem)

    if vx <= Problem.VaryingParam{1}(1)
        vx = Problem.VaryingParam{1}(1)+1e-4;
    end

    if vx >= Problem.VaryingParam{1}(end)
        vx = Problem.VaryingParam{1}(end)-1e-4;
    end

    if vx < 0.9
        
        VxMIN = 0.6;
        VxMAX = 0.9;
        
        O1 = (vx-VxMIN)/(VxMAX-VxMIN);

        Ki = K{1}*(1-O1) + K{2}*O1;

    elseif vx < 1.2
        
        VxMIN = 0.9;
        VxMAX = 1.2;

        O1 = (vx-VxMIN)/(VxMAX-VxMIN);

        Ki = K{2}*(1-O1) + K{3}*O1;

    elseif vx < 1.5
        
        VxMIN = 1.2;
        VxMAX = 1.5;

        O1 = (vx-VxMIN)/(VxMAX-VxMIN);

        Ki = K{3}*(1-O1) + K{4}*O1;

    else
        
        VxMIN = 1.5;
        VxMAX = 1.8;
        
        O1 = (vx-VxMIN)/(VxMAX-VxMIN);

        Ki = K{4}*(1-O1) + K{5}*O1;

    end

    u = Ki*x;
    
end

