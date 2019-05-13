%% State space representation for Aucher Model [Psi V R]

function model = tbx_sm_heading_model(p)
    
    M  = [1 0         0             ; ...
          0 p.Mu      p.Mu*p.Xg     ; ...
          0 p.Mu*p.Xg p.L*p.L*p.X_3];
    
    Ma = [0 0            0                  ; ...
          0 p.Mu_2       -p.L*p.Nu_26       ; ...
          0 -p.L*p.Nu_26 p.L*p.L*p.Lambda_3];

    Masse = M + Ma;

    Ainert = [0 0 1                ; ...
              0 0 -p.Mu*p.Vs       ; ...
              0 0 -p.Mu*p.Xg*p.Vs] ;

    A      = [0 0         0        ; ...
              0 p.CyV/p.L p.CyR    ; ...
              0 p.CnV     p.CnR*p.L]*p.Vs/(2*p.CRemp);

    A = A + Ainert;

    B = [0 ; p.CyAL/p.L ; p.CnAL]*p.Vs*p.Vs/(2*p.CRemp) ;

    A = Masse\A;
    B = Masse\B;

    
    model.a = A;
    model.b = B;
    model.c = eye(3);
    model.d = zeros(3,1);
