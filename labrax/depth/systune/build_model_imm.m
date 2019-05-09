function model = build_model_imm(nav,Vs)
 
    p.Vs    = Vs ;
	
    p.L        = nav.L        ;
    p.mu       = ureal('Mu', 0.99, 'Percentage', 0.01);
    p.g        = nav.g        ;
    p.CRemp    = nav.CRemp    ;
    p.Xg       = ureal('Xg', 0, 'Range', [-0.01 0.01]);
    p.Zg       = nav.Zg       ;
    p.X_2      = nav.X_2      ;
    p.Mu_3     = nav.Mu_3     ;
    p.Nu_35    = nav.Nu_35    ;
    p.Lambda_2 = nav.Lambda_2 ;
	
    p.CzW  = ureal('CzW', nav.CzW, 'Percentage', 20);
	p.CmQ  = ureal('CmQ', nav.CmQ, 'Percentage', 20);
	
    p.CzQ      = nav.CzQ  ;
    p.CzB1     = ureal('CzB1', nav.CzB1, 'Percentage', 20);
    p.CmW      = nav.CmW  ;
    p.CmB1     = ureal('CmB1', nav.CmB1, 'Percentage', 20);
    
    %%[Z Theta W Q beta1 0]

    Cz = [0 0 p.CzW/p.Vs p.CzQ*p.L/p.Vs p.CzB1 0];
    Cm = [0 0 p.CmW/p.Vs p.CmQ*p.L/p.Vs p.CmB1 0];

    Fz_navire = p.Vs^2/(2*p.CRemp*p.L)*Cz;
    My_navire = p.Vs^2/(2*p.CRemp)*Cm;

    Fz_couplage = [0 0 0       p.mu*p.Vs 0 0];
    My_couplage = [0 0 0 -p.mu*p.Xg*p.Vs 0 0];

    Fz_hydrostat = [0              0 0 0 0   (p.mu-1)*p.g];
    My_hydrostat = [0 -p.mu*p.Zg*p.g 0 0 0 -p.mu*p.Xg*p.g];

    Fz = Fz_hydrostat + Fz_couplage + Fz_navire;
    My = My_hydrostat + My_couplage + My_navire;

    Mpropre  = [      p.mu  -p.mu*p.Xg  ;...
                -p.mu*p.Xg p.L^2*p.X_2] ;
                
    Majoutee = [       p.Mu_3      -p.L*p.Nu_35  ;...
                 -p.L*p.Nu_35  p.L^2*p.Lambda_2] ;
                 
    Masse = Mpropre + Majoutee;

    Aux   = inv(Masse)*[Fz ; My];

    Wpt   = Aux(1,:);
    Qpt   = Aux(2,:);
    
    thetapt = [0     0 0 1 0 0];
    Zpt     = [0 -p.Vs 1 0 0 0];

    F = [Zpt ; thetapt ; Wpt ; Qpt];
    A = F(:,1:4);
    B = F(:,5);

    ZThetaPiQ = [1 0 0 0 ; 0 1 0 0 ; 0 1 -1/p.Vs 0 ; 0 0 0 1]; %% Transformee en X = [Theta, Pente, Q];

    A = ZThetaPiQ*A*inv(ZThetaPiQ);
    B = ZThetaPiQ*B;

	sys.a = simplify(A);
	sys.b = simplify(B);
    sys.c = eye(4);
    sys.d = zeros(4,1);

	%% Save in struct
    model            = ss(sys.a, sys.b, sys.c, sys.d);
    model.InputName  = {'Betaf'};
    model.OutputName = {'Z','Theta','Pi','Q'};
