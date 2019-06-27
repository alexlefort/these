function transferts = build_transferts(model_sys, model_ctrl)

    syms s;
    
    %% Construire la fonction de transfert en immersion
    
    %% X = [dZ DTheta DPi DQ]
    %% Y = [dZ DTheta DPi DQ]
    %% U = [Beta1 Beta2]

	K  = model_ctrl.symtf; %% controller transfert function
    G  = model_sys.G;      %% system transfert function
    
    [n,~] = size(G*K);
    %% Sensitivity
    H1 = inv(eye(n) + G*K);
    H2 = K*H1;
    
    transferts.zz  = simplify(H1(1,1));
    transferts.zt  = simplify(H1(1,2));
    transferts.tz  = simplify(H1(2,1));
    transferts.tt  = simplify(H1(2,2));
    transferts.zb1 = simplify(H2(1,1));
    transferts.zb2 = simplify(H2(2,1));
    transferts.tb1 = simplify(H2(1,2));
    transferts.tb2 = simplify(H2(2,2));
    
    [transferts.poly_stab,transferts.mat_stab] = symtbx_closed_loop_poly(model_sys.sys,model_ctrl.symss); %% Intern stability