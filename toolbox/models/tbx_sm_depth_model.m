function model = tbx_sm_depth_model_one_plane(p)

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

    Mpropre  = [      p.mu   -p.mu*p.Xg ;...
                -p.mu*p.Xg p.L^2*p.X_2] ;
                
    Majoutee = [      p.Mu_3       -p.L*p.Nu_35 ;...
                -p.L*p.Nu_35  p.L^2*p.Lambda_2] ;
                 
    Masse = Mpropre + Majoutee;

    Aux = Masse\[Fz ; My];
    
    Wpt = Aux(1,:);
    Qpt = Aux(2,:);
    
    thetapt = [0     0 0 1 0 0];
    Zpt     = [0 -p.Vs 1 0 0 0];

    F = [Zpt ; thetapt ; Wpt ; Qpt];
    A = F(:,1:4);
    B = F(:,5);

    %% Conversion (Z Theta W Q) -> (Z Theta Pi Q)

    ZThetaPiQ = sym([1 0 0 0 ; 0 1 0 0 ; 0 1 -1/p.Vs 0 ; 0 0 0 1]); %% Transformee en X = [Theta, Pente, Q];

    A = ZThetaPiQ*A*inv(ZThetaPiQ);
    B = ZThetaPiQ*B;
    
    sys.a = simplify(A);
    sys.b = simplify(B);
    sys.c = sym(eye(4));
    sys.d = sym(zeros(4,1));
