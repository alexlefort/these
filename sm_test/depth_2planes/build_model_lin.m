function model = build_model_lin()
 
    p.Vs        = sym('Vs'      ,'real'); p.Vs       = sym( '20000/10000');
    p.L         = sym('L'       ,'real'); p.L        = sym( '69601/10000');
    p.mu        = sym('mu'      ,'real'); p.mu       = sym( '10000/10000');
    p.g         = sym('g'       ,'real'); p.g        = sym( '98100/10000');
    p.CRemp     = sym('CRemp'   ,'real'); p.CRemp    = sym( '09750/10000');
    p.Xg        = sym('Xg'      ,'real'); p.Xg       = sym( '00000/10000');
    p.Zg        = sym('Zg'      ,'real'); p.Zg       = sym( '00231/10000'); 
         
    p.CzW       = sym('CzW'     ,'real'); p0.CzW     = sym('-27000/10000');
    p.CzQ       = sym('CzQ'     ,'real'); p.CzQ      = sym('-06020/10000');
    p.CzB1      = sym('CzB1'    ,'real'); p.CzB1     = sym('-07030/10000');
    p.CzB2      = sym('CzB2'    ,'real'); p.CzB2     = sym('-05040/10000');
    p.CmW       = sym('CmW'     ,'real'); p.CmW      = sym( '08300/10000');
    p.CmQ       = sym('CmQ'     ,'real'); p0.CmQ     = sym('-03800/10000');
    p.CmB1      = sym('CmB1'    ,'real'); p.CmB1     = sym('-03230/10000');
    p.CmB2      = sym('CmB2'    ,'real'); p.CmB2     = sym('-00960/10000');

    p.X_2       = sym('X_2'     ,'real'); p.X_2      = sym( '00665/10000'); 
    p.Mu_3      = sym('Mu_3'    ,'real'); p.Mu_3     = sym( '08780/10000');
    p.Nu_35     = sym('Nu_35'   ,'real'); p.Nu_35    = sym('-00114/10000');
    p.Lambda_2  = sym('Lambda_2','real'); p.Lambda_2 = sym( '00667/10000');

    %%[Z Theta W Q beta1 0]

    Cz = sym([0 0 p.CzW/p.Vs p.CzQ*p.L/p.Vs p.CzB1 p.CzB2 0]);
    Cm = sym([0 0 p.CmW/p.Vs p.CmQ*p.L/p.Vs p.CmB1 p.CmB2 0]);

    Fz_navire = sym(p.Vs^2/(2*p.CRemp*p.L)*Cz);
    My_navire = sym(p.Vs^2/(2*p.CRemp)*Cm);

    Fz_couplage = sym([0 0 0       p.mu*p.Vs 0 0 0]);
    My_couplage = sym([0 0 0 -p.mu*p.Xg*p.Vs 0 0 0]);

    Fz_hydrostat = sym([0              0 0 0 0 0   (p.mu-1)*p.g]);
    My_hydrostat = sym([0 -p.mu*p.Zg*p.g 0 0 0 0 -p.mu*p.Xg*p.g]);

    Fz = sym(Fz_hydrostat + Fz_couplage + Fz_navire);
    My = sym(My_hydrostat + My_couplage + My_navire);

    Mpropre  = sym([      p.mu   -p.mu*p.Xg  ;...
                    -p.mu*p.Xg p.L^2*p.X_2]) ;
                
    Majoutee = sym([       p.Mu_3       -p.L*p.Nu_35  ;...
                     -p.L*p.Nu_35  p.L^2*p.Lambda_2]) ;
                 
    Masse = sym(Mpropre + Majoutee);

    Aux = sym(inv(Masse)*[Fz ; My]);
    
    Wpt = sym(Aux(1,:));
    Qpt = sym(Aux(2,:));
    
    thetapt = sym([0     0 0 1 0 0 0]);
    Zpt     = sym([0 -p.Vs 1 0 0 0 0]);

    F = sym([Zpt ; thetapt ; Wpt ; Qpt]);
    A = sym(F(:,1:4));
    B = sym(F(:,5:6));
    

	sys.a = simplify(A);
	sys.b = simplify(B);
	sys.c = sym(eye(4));
	sys.d = sym(zeros(4,2));

	%% Build central model

	sys0.a = eval(symtbx_sym_subs_from_struct(sys.a, p0));
	sys0.b = eval(symtbx_sym_subs_from_struct(sys.b, p0));
	sys0.c = eval(symtbx_sym_subs_from_struct(sys.c, p0));
	sys0.d = eval(symtbx_sym_subs_from_struct(sys.d, p0));
	
	syms s;
	
	G  = simplify(symtbx_symss_to_symtf(sys));
	G0 = simplify(symtbx_symss_to_symtf(sys0));

	%% Save in struct
	model.sys   = sys  ;
	model.p     = p    ;
	model.G     = G    ;
    model.sys0  = sys0 ;
    model.p0    = p0   ;
	model.G0    = G0   ;
