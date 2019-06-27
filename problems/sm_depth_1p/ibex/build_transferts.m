function transferts = build_transferts(model_sys, model_ctrl)

    syms s;
    
    %% Construire la fonction de transfert en immersion
    
    %% X = [Z Theta Pi Q]
    %% Y = [Z Theta Pi Q]
    %% U = [Beta1]

	K  = model_ctrl.symtf; %% controller transfert function
    G  = model_sys.G;      %% system transfert function
    [n,~] = size(G*K);
    
    %% Sensitivity
    H1 = inv(eye(n) + G*K);
    H2 = K*H1;
    
    zz = H1(1,1);
    zt = H1(2,1);
    zb = H2(1,1);

    transferts.zz = simplify(zz);
    transferts.zt = simplify(zt);
    transferts.zb = simplify(zb);
    
    [transferts.poly_stab,transferts.mat_stab] = symtbx_closed_loop_poly(model_sys.sys,model_ctrl.symss); %% Intern stability