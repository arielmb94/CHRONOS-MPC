function dt_x = TRMS(Wh,Omh,Thth,Wv,Omv,Thtv,uh,uv)

    k_g = 0.2;
    l_t= 0.282;
    l_m = 0.246;
    l_b = 0.29;
    l_cb = 0.276;
    r_ts = 0.1;
    r_ms = 0.155;
    m_cb = 0.068;
    m_b = 0.022;
    m_m = 0.014;
    m_mr = 0.236;
    m_ms = 0.219;
    L_ah = 0.86e-3;
    L_av = 0.86e-3;
    k_a = 0.0202;
    B_tr = 0.0086;
    B_mr = 0.0026;
    k_tvp = 0.0168;
    k_tvn = 0.0155;
    k_fhp = 0.0566;
    k_fhn = 0.0660;
    k_fvp = 0.3819;
    k_fvn = 0.2197;
    k_chp = 0.0158;
    k_chn = 0.0111;
    k_cvp = 0.0623;
    m_t = 0.015;
    m_tr = 0.221;
    m_ts = 0.119;
    J_tr = 0.0059;
    J_mr = 0.0254;
    R_a = 8;
    k_t = 0;
    k_m = 0.0017;
    k_1 = 6.5;
    k_2 = 8.5;
    k_oh = 0.0185; 
    k_ov = 0.1026;
    k_cvn = 0.0563;
    g = 9.81;
    
    K_A = (0.5*m_t+m_tr+m_ts)*l_t;
    K_B = (0.5*m_m+m_mr+m_ms)*l_m;
    K_C = 0.5*m_b*l_b+m_cb*l_cb;
    K_D = (m_m/3+m_mr+m_ms)*(l_m^2)+(m_t/3+m_tr+m_ts)*(l_t^2);
    K_E = m_b*(l_b^2)/3+m_cb*(l_cb^2);
    K_F = m_ms*r_ms^2+0.5*m_ts*r_ts^2;
    K_H = K_A*l_t+K_B*l_m+0.5*m_b*(l_b^2)+m_cb*l_cb;
    
    dt_Wh = k_a*k_1/(J_tr*R_a)*uh-(B_tr/J_tr+(k_a^2)/(J_tr*R_a))*Wh-f1(Wh)/J_tr;
    dt_Omh = (l_t*f2(Wh)*cos(Thtv)-k_oh*Omh-f3(Thth)+f6(Thtv))/(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F)+...
        k_m*Wv*sin(Thtv)*Omv*(K_D*(cos(Thtv)^2)-K_E*(sin(Thtv)^2)-K_F-2*K_E*(cos(Thtv)^2))/((K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F)^2)+...
        k_m*cos(Thtv)*((k_a*k_2/R_a)*uv-(B_mr+(k_a^2)/R_a)*Wv-f4(Wv))/(J_mr*(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F));
    dt_Thth = Omh;
    
    Jv = m_mr*l_m^2+m_m*(l_m^2)/3+m_cb*l_cb^2+m_b*(l_b^2)/3+m_tr*l_t^2+m_t*(l_t^2)/3+...
        0.5*m_ms*r_ms^2+m_ms*l_m^2+0.5*m_ts*r_ts^2+m_ts*l_t^2;
    
    dt_Wv = k_a*k_2/(J_mr*R_a)*uv-(B_mr/J_mr+(k_a^2)/(J_mr*R_a))*Wv-f4(Wv)/J_mr;
    dt_Omv = (l_m*f5(Wv)+k_g*Omh*f5(Wv)*cos(Thtv)-k_ov*Omv)/Jv+...
        (g*((K_A-K_B)*cos(Thtv)-K_C*sin(Thtv))-(Omh^2)*K_H*sin(Thtv)*cos(Thtv))/Jv;
    dt_Thtv = Omv;     

    dt_x = [dt_Wh,dt_Omh,dt_Thth,dt_Wv,dt_Omv,dt_Thtv]';
end


function y = f1(Wh)

    k_thp = 0.0027;
    k_thn = 0.0028;

    if Wh>=0
        y = k_thp*Wh^2;
    else
        y = -k_thn*Wh^2;
    end
end

function y = f2(Wh)

    k_fhp = 0.0566;
    k_fhn = 0.0660;

    if Wh>=0
        y = k_fhp*Wh^2;
    else
        y = -k_fhn*Wh^2;
    end
end

function y = f3(Thth)

    k_chp = 0.0158;
    k_chn = 0.0111;

    if Thth>=0
        y = k_chp*Thth;
    else
        y = k_chn*Thth;
    end
end

function y = f4(Wv)

    k_tvp = 0.0168;
    k_tvn = 0.0155;

    if Wv>=0
        y = k_tvp*Wv^2;
    else
        y = -k_tvn*Wv^2;
    end
end

function y = f5(Wv)

    k_fvp = 0.3819;
    k_fvn = 0.2197;

    if Wv>=0
        y = k_fvp*Wv^2;
    else
        y = -k_fvn*Wv^2;
    end
end

function y = f6(Thtv)

    k_cvp = 0.0623;
    k_cvn = 0.0563;
    
    Thtv0 = -0.619426178368110;

    if Thtv>=Thtv0
        y = k_cvp*(Thtv-Thtv0)^2;
    else
        y = k_cvn*(Thtv-Thtv0)^2;
    end
end