function p = reglage_ctrl_gir_rco()

    nav   = charger_modele_sm('iver2.nav');
    deg   = pi/180;
    UVect = [0.3;0.5;1.0;1.5;2.0;2.5;3.0;4.0;5.0];
    n     = length(UVect);

    p.sat_Alpha_min = -20*deg;
    p.sat_Alpha_max =  20*deg;
    p.sat_Alpha_d   =  15*deg;
    p.sat_RCo_max   =   3*deg;
    
    pulse        = pi;
    p.puls_Psi   = pulse;
    p.puls_Alpha = pulse;
    
    KPsiVect  = zeros(n,1);
    KRVect    = zeros(n,1);
    IRVect    = zeros(n,1);

    opt = systuneOptions('RandomStart', 10  , ...
                         'UseParallel', true, ...
                         'Display', 'off');

    for ii=1:n
    
	    Vs = UVect(ii);
        disp(Vs);
	    %% Build linear model 
	    sm = build_model_gir(nav, Vs);   
	
        %% Build controller model
        g.kpsi   = realp('kpsi'  ,0);
        g.kr     = realp('kr'    ,0);
        g.ir     = realp('ir'    ,0);
        
        ctrl_psi            = tf(g.kpsi,1);
        ctrl_psi.InputName  = {'ePsif'};
        ctrl_psi.OutputName = {'rco'};
        
        ctrl_r            = tf(g.kr,1)+tf(g.ir,[1 0]);
        ctrl_r.InputName  = {'er'};
        ctrl_r.OutputName = {'Alpha'};
        
        ctrlf            = tf(1, [1/p.puls_Alpha 1]);
        ctrlf.InputName  = 'Alpha';
        ctrlf.OutputName = 'Alphaf';
	
        %% Build filter model
        ppsi = tf(1    , [1/(p.puls_Psi^2) sqrt(2)/p.puls_Psi 1]);
        dpsi = tf([1 0], [1/(p.puls_Psi^2) sqrt(2)/p.puls_Psi 1]);
        
        filter            = [ppsi ; dpsi];
        filter.InputName  = 'ePsi';
        filter.OutputName = {'ePsif','dePsif'};
        
        %% Build Weightings
        sm0 = sm.NominalValue;
        fc = 0.1*find_cutoff_frequency(sm0(1,1));
              
        wpsi = mk_pondw0(-40,fc/(2*pi),2,4,'wPsi');

        wpsi.InputName  = 'ePsi';
        wpsi.OutputName = 'wPsi';
        
        wa = tf(1.6,1);
        wa.InputName  = 'Alphaf';
        wa.OutputName = 'wAlpha';

        %% Build transfert functions
        refPsi   = sumblk('ePsi   = rPsi  - Psi'     );
        refPsif  = sumblk('er     = rco   + dePsif'  );

        %%figure; hold on; bode(sm(1,1)) ; bode(wpsi); grid on;
        sys = connect(sm,filter,ctrl_psi, ctrl_r,ctrlf,refPsi, refPsif,{'rPsi'}, {'ePsi', 'Alphaf'});
    
        perf = TuningGoal.WeightedGain({'rPsi'}, {'ePsi'}, [], 1/wpsi) ;  perf.focus = [1e-5 10];
        sat  = TuningGoal.WeightedGain({'rPsi'}, {'Alphaf'}, [], 1/wa) ;  sat.focus  = [1e-5 10];
    
        [CL,~] = systune(sys, [perf,sat], opt);
        
        %%sysPsi   = connect(CL,{'rPsi'},{'ePsi'});
        %%sysAlpha = connect(CL,{'rPsi'},{'Alpha'});
        
        res = getBlockValue(CL);
        
        KPsiVect(ii)  = res.kpsi  ;
        KRVect(ii)    = res.kr    ;
        IRVect(ii)    = res.ir     ;
    
        %%figure;  subplot(1,2,1);  stepplot(sysPsi)   ;
        %%         subplot(1,2,2);  stepplot(sysAlpha) ;
        %%figure;  hold on ; bode(sysPsi)   ; bode(wpsi) ;
        %%figure;  hold on ; bode(sysAlpha) ; bode(wa)   ;
        %%showTunable(CL);

    end

    p.gains_vitesses    = UVect;
    p.gains_IR_vect     = IRVect;
    p.gains_KPPSI_vect  = KPsiVect;
    p.gains_KR_vect     = KRVect;
