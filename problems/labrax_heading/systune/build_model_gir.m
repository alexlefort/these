function model = build_model_gir(nav,Vs)
 
    p.Vs        = Vs           ;
    p.L         = nav.L        ;
    p.Mu        = ureal('Mu', 1, 'Percentage', 0.01);
    p.g         = nav.g        ;
    p.CRemp     = nav.CRemp    ;
    p.Xg        = ureal('Xg', 0, 'Range', [-0.01 0.01]);
    p.Zg        = nav.Zg       ;
    p.Mu_2      = nav.Mu_2     ;
    p.Nu_24     = nav.Nu_24    ;
    p.Nu_26     = nav.Nu_26    ;
    p.X_1       = nav.X_1      ;
    p.X_3       = nav.X_3      ;
    p.X_13      = nav.X_13     ;
    p.Lambda_1  = nav.Lambda_1 ;
    p.Lambda_3  = nav.Lambda_3 ;
    p.Lambda_13 = nav.Lambda_13;
	
    p.CyV      = ureal('CyV', nav.CyV, 'Percentage', 20);
    p.CyP      = nav.CyP;
    p.CyR      = nav.CyR;
    p.ClV      = nav.ClV;
    p.ClP      = ureal('ClP' , nav.ClP , 'Percentage', 20);
    p.ClR      = nav.ClR;
    p.CnV      = nav.CnV;
    p.CnP      = nav.CnP;
    p.CnR      = ureal('CnR' , nav.CnR , 'Percentage', 20);
    p.CyAL     = ureal('CyAL', nav.CyAL, 'Percentage', 20);
    p.ClAL     = nav.ClAL;
    p.CnAL     = ureal('CnAL', nav.CnAL, 'Percentage', 20);
    
    M = zeros(3,3);
    M(1,1) = p.Mu+p.Mu_2;
    M(1,2) = -p.Mu*p.Zg-p.L*p.Nu_24;
    M(1,3) = p.Mu*p.Xg-p.L*p.Nu_26;
    M(2,1) = -p.Mu*p.Zg-p.L*p.Nu_24;
    M(2,2) = p.L*p.L*(p.X_1+p.Lambda_1);
    M(2,3) = -p.L*p.L*(p.X_13+p.Lambda_13);
    M(3,1) = p.Mu*p.Xg-p.L*p.Nu_26;
    M(3,2) = -p.L*p.L*(p.X_13+p.Lambda_13);
    M(3,3) = p.L*p.L*(p.X_3+p.Lambda_3);
   
    Minvlat2 = [eye(2) zeros(2,3) ; zeros(3,2) inv(M)] ;
 
    Ainert2 = [  0 0 0 1 0 ; 0 0 0 0 1 ; 0 0 0 0 -p.Mu*p.Vs ; 0 0 0 0 p.Mu*p.Zg*p.Vs ; 0 0 0 0 -p.Mu*p.Xg*p.Vs  ] ;
    Ahydrostat2 = [ zeros(2,5) ; (p.Mu-1)*p.g 0 0 0 0 ; -p.Mu*p.Zg*p.g 0 0 0 0 ; p.Mu*p.Xg*p.g 0 0 0 0] ;
    Aetat2 = Vs/(2*p.CRemp)*[ zeros(2,5) ; 0 0 p.CyV/p.L p.CyP p.CyR  ;  0 0 p.ClV p.ClP*p.L p.ClR*p.L ; 0 0 p.CnV p.CnP*p.L p.CnR*p.L ] ;
    A = Minvlat2*(Ainert2+Aetat2+Ahydrostat2) ;
    B = Minvlat2*p.Vs*p.Vs/(2*p.CRemp)*[0 ; 0 ; p.CyAL/p.L ; p.ClAL ; p.CnAL] ;
    C = [0 1 0 0 0];
    D = 0;

	%% Save in struct
    model            = ss(A, B, C, D);
    model.InputName  = {'Alphaf'};
    model.OutputName = {'Psi'};
