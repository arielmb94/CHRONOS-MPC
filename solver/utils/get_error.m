function error = get_error(r,y,D,N,Ny)

    if length(r)<Ny
        if D == 0 
            r = repmat(r,N-1,1);
        else
            r = repmat(r,N,1);
        end
    end    

    error = r-y;
        
end