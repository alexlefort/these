function p = reglage_ctrl_imm(nav,bruit)

    warning('off');
    deg = pi/180;
    
    UVect  = [0.3;0.5;1.0;1.5;2.0;2.5;3.0;4.0;5.0];
    n      = length(UVect);
    deg = pi/180;

    p.sat_Beta    = 20*deg        ;
    p.sat_Beta_d  = 12*deg        ;
    p.sat_Pi_min  = -15*deg       ;
    p.sat_Pi_max  = 15*deg        ;
    p.sat_Z_max   = 100           ;
    p.sat_Z_min   = 1             ;
    p.gain_antiwp = -0.5          ;
    p.puls_Beta   = pi;

    p.gains_KZ     = zeros(n,1);
    p.gains_IZ     = zeros(n,1);
    p.gains_KPi    = zeros(n,1);
    p.gains_KTheta = zeros(n,1);
    p.gains_KQ     = zeros(n,1);

    opt = systuneOptions('RandomStart', 10  , ...
                         'UseParallel', true, ...
                         'Display', 'iter');
    
    bruit_barre_obj = 0.5*deg;
    coeff_att_z     = 1/(bruit_barre_obj/(bruit.Z*3));
    coeff_att_theta = 1/(bruit_barre_obj/(bruit.Theta*3));
    coeff_att_z
    coeff_att_theta

    for ii=4
    
        Vs = UVect(ii);
        disp(Vs);
        
	    %% Build linear model 
	    sm = build_model_imm(nav, Vs); 
        
        %% Sensors noise model

        %% Build controller model
	    kz      = realp('kz'     ,0);
        iz      = realp('iz'     ,0);
        kpi     = realp('kpi'    ,0);
        ktheta  = realp('ktheta' ,0);
        kq      = realp('kq'     ,0);
        
        ctrl            = [(tf(ktheta,1)+tf(kpi,1))*(tf(kz,1)+tf(iz,[1 0])) tf(ktheta,1) tf(kpi,1) tf(kq,1)];
        ctrl.InputName  = {'eZ','eTheta','ePi','eQ'};
        ctrl.OutputName = {'Beta'};
        
        ctrlf            = tf(1, [1/p.puls_Beta 1]);
        ctrlf.InputName  = 'Beta';
        ctrlf.OutputName = 'Betaf';
        
        %% Build weightings
        sm0 = sm.NominalValue;
        fc = 0.05*find_cutoff_frequency(sm0(1,1));
              
        wz = mk_pondw0(-60,fc/(2*pi),4,4,'wz');
        wz.InputName  = 'eZ';
        wz.OutputName = 'wZ';

        wbz = 1/coeff_att_z*tf(1,1);
        wbz.InputName  = 'Beta' ;
        wbz.OutputName = 'wBetaz'; 

        wbt = 1/coeff_att_theta*tf(1,1);
        wbt.InputName  = 'Beta' ;
        wbt.OutputName = 'wBetat';
        
        %% Build transfert functions
        refZ     = sumblk('eZ     = rZ  - Z - Z_s'    );
        refTheta = sumblk('eTheta =     - Theta - Theta_s');
        refPi    = sumblk('ePi    =     - Pi'   );
        refQ     = sumblk('eQ     =     - Q'    );
      
        %%figure; hold on; bode(sm(1,1)) ; bode(wz); grid on;
        sys = connect(sm, ctrl, ctrlf, refZ, refTheta, refPi, refQ, ...
            {'rZ','Z_s','Theta_s'}, {'eZ','Beta'});
    
        perf = TuningGoal.WeightedGain('rZ'      , 'eZ'   , [], 1/wz );   perf.focus = [1e-5 10];
        satz = TuningGoal.WeightedGain('Z_s'     , 'Beta' , [], 1/wbz);   satz.focus = [1e-5 10]; 
        satt = TuningGoal.WeightedGain('Theta_s' , 'Beta' , [], 1/wbt);   satt.focus = [1e-5 10]; 
        
        %% Appel Systune 
        [CL, ctr] = systune(sys, [perf, satz, satt], [satz, satt], opt);
        res = getBlockValue(CL);
        
        sysZ  = connect(CL , 'rZ'      , 'eZ'  );
        sysBZ = connect(CL , 'Z_s'     , 'Beta');
        sysBT = connect(CL , 'Theta_s' , 'Beta'); 
        
        p.gains_KZ    (ii) = res.kz     ;
        p.gains_IZ    (ii) = res.iz     ;
        p.gains_KPi   (ii) = res.kpi    ;
        p.gains_KTheta(ii) = res.ktheta ;
        p.gains_KQ    (ii) = res.kq     ;
        ctrVect(ii) = max(ctr);
        
        figure;  subplot(1, 3, 1);  stepplot(sysZ);
                 subplot(1, 3, 2);  impulseplot(sysBZ);
                 subplot(1, 3, 2);  impulseplot(sysBT);
                 
        figure;  bode(sysZ) ; hold on; bode(wz);
        figure;  bode(sysBZ); hold on; bode(wbz);
        figure;  bode(sysBT); hold on; bode(wbt);

        zthpiq = ss([0 tf(-Vs,[1 0]) 0]) ;
        
        p_ctrl.KQ     = p.gains_KQ(ii)     ;
        p_ctrl.KTheta = p.gains_KTheta(ii) ;
        p_ctrl.KPi    = p.gains_KPi(ii)    ;
        p_ctrl.Kz     = p.gains_KZ(ii)     ;
        p_ctrl.Iz     = p.gains_IZ(ii)     ;
        
        p_ctrl.sat_Beta    = p.sat_Beta    ;
        p_ctrl.sat_Beta_d  = p.sat_Beta_d  ;
        p_ctrl.sat_Pi_min  = p.sat_Pi_min  ;
        p_ctrl.sat_Pi_max  = p.sat_Pi_max  ;
        p_ctrl.gain_antiwp = p.gain_antiwp ;
        p_ctrl.omega       = pi            ;
        p_ctrl.dt          = 0.1           ;
        p_ctrl.ZOrder      = 5.0           ;
 
        r = ctrl_imm_charger_param;
        [A,B,C,D] = modele_lineaire(Vs, nav, 'IMM', r);
        sys_pb = ss(A,B,C,D);
        
        p_sys.zthpiq  = zthpiq;
        p_sys.sm      = sys_pb;

        simulation(p_sys,p_ctrl,200);
        
    end

    p.gains_vitesses = UVect ;

    disp(ctrVect);

