function p_ctrl = ctrl_imm_calculer_reglages(p)

    r = ctrl_imm_charger_param;

    p_ctrl.sat_Beta    = r.sat_Beta    ;
    p_ctrl.sat_Beta_d  = r.sat_Beta_d  ;
    p_ctrl.sat_Pi_min  = r.sat_Pi_min  ;
    p_ctrl.sat_Pi_max  = r.sat_Pi_max  ;
    p_ctrl.sat_Z_max   = r.sat_Z_max   ;
    p_ctrl.sat_Z_min   = r.sat_Z_min   ;
    p_ctrl.gain_antiwp = r.gain_antiwp ;
    p_ctrl.puls_Beta   = pi            ;

    p_ctrl.gains_KZ     = zeros(r.nVs,1);
    p_ctrl.gains_IZ     = zeros(r.nVs,1);
    p_ctrl.gains_KPi    = zeros(r.nVs,1);
    p_ctrl.gains_KTheta = zeros(r.nVs,1);
    p_ctrl.gains_KQ     = zeros(r.nVs,1);
    
    for ii = 1:r.nVs
    
        Vs   = r.Vs_Vect(ii);        
        Trep = 5*r.Trep_Vect(ii);
        
        [A,B,C,D] = modele_lineaire(Vs, p, 'IMM', r);
        
        sys_pb = ss(A,B,C,D);
        
        % Augmentation vers modele [IZ Z Theta Pente Q]
        A = [0 1 0 0 0 ; 0 0 0 -Vs 0 ; zeros(3,2) A];
        B = [0 ; 0 ; B];
        C = eye(5);
        D = zeros(5,1);

        ssm = ss(A,B,C,D);
        
        % Reglage du stabilisateur par retour LQ
    
        matQ  = diag([5 1 1 100 1]);
        matR  = Trep*ssm.b'*grampar(ssm,Trep,0.1)*ssm.b;
        Kstab = lqr(A,B,matQ,matR);
  
        p_ctrl.gains_IZ(ii)      = Kstab(1,1)/(Kstab(1,4)+Kstab(1,3));
        p_ctrl.gains_KZ(ii)      = Kstab(1,2)/(Kstab(1,4)+Kstab(1,3));
        p_ctrl.gains_KTheta(ii)  = Kstab(1,3);
        p_ctrl.gains_KPi(ii)     = Kstab(1,4);
        p_ctrl.gains_KQ(ii)      = Kstab(1,5);

        zthpiq = ss([0 tf(-Vs,[1 0]) 0]) ;
        
        ctrl.KQ     = p_ctrl.gains_KQ(ii)     ;
        ctrl.KTheta = p_ctrl.gains_KTheta(ii) ;
        ctrl.KPi    = p_ctrl.gains_KPi(ii)    ;
        ctrl.Kz     = p_ctrl.gains_KZ(ii)     ;
        ctrl.Iz     = p_ctrl.gains_IZ(ii)     ;
        
        ctrl.sat_Beta    = r.sat_Beta    ;
        ctrl.sat_Beta_d  = r.sat_Beta_d  ;
        ctrl.sat_Pi_min  = r.sat_Pi_min  ;
        ctrl.sat_Pi_max  = r.sat_Pi_max  ;
        ctrl.gain_antiwp = r.gain_antiwp ;
        ctrl.omega       = pi            ;
        ctrl.dt          = 0.1           ;
        ctrl.ZOrder      = 5.0           ;
    
        sys.zthpiq  = zthpiq;
        sys.sm      = sys_pb;

        simulation(sys,ctrl,200);
        
    end
