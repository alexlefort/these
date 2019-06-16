function transferts = build_transferts(model_sys, model_ctrl)

    syms s;
    
    %% Construire la fonction de transfert en immersion
    
    %% X = [Psi V R]
    %% Y = [Psi V R]
    %% U = [Alpha]

	K  = model_ctrl.symtf; %% controller transfert function
    G  = model_sys.G;      %% system transfert function

    %% Sensitivity
    H1 = inv(eye(3) + G*K);
    H2 = K*H1;
    
    psi  = H1(1,1);
    psia = H2(1,1);

    transferts.psi  = simplify(psi); %% Closed loop transfert function for the states
    transferts.psia = simplify(psia); %% Closed loop transfert function for the commands
    [transferts.poly_stab,transferts.mat_stab] = symtbx_closed_loop_poly(model_sys.sys,model_ctrl.symss); %% Intern stability