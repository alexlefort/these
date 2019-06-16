function  f = model_capd()

    Vs = 1.5;
	model_sm   = build_model_lin(Vs);
    model_ctrl = build_controller(Vs);
    ZCo        = sym('ZCo','real');
    xCo        = sym([ZCo ; 0 ; 0 ; 0]);
    Z          = sym('Z','real');
    Theta      = sym('Theta','real');
    Pi         = sym('Pi','real');
    Q          = sym('Q','real');
    IZ         = sym('IZ','real');
    Beta       = sym('Beta','real');
    X          = sym([Z;Theta;Pi;Q;IZ;Beta]);
    
    A = model_sm.sys.a;
    B = model_sm.sys.b;
    
    Ak = model_ctrl.symss.a;
    Bk = model_ctrl.symss.b;
    Ck = model_ctrl.symss.c;
    Dk = model_ctrl.symss.d;
    
    M1  = sym([A-B*Dk B*Ck ; -Bk Ak]);
    M2  = sym([B*Dk*xCo ; Bk*xCo]);
    f   = M1*X + M2;
    
    res = ['"par:CmQ,CzW,CmB1,CzB1,ZCo,kpi,kq,ktheta,kz,iz;var:Z,Theta,Pi,Q,IZ,Beta;f:'];
    for ii=1:5
        res = [res char(f(ii)) ','];
    end
    res = [res char(f(6)) '"' ];
    file = fopen('res.txt','w');
    fprintf(file,res);
    fclose(file);
    
end


