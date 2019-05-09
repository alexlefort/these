function test_stab_nyquist(poly,param,itermax)

    poly = symtbx_fix_real_param(poly,param);
    rho  = symtbx_get_gershgorhin_radius(poly,param);
    
    diff_poly = diff(poly,'s');
    integ     = diff_poly/poly;
    w = sym('w','real');
    g1 = subs(g, 's', 1i*w);
    g2 = subs(g, 's',1e5*exp(1i*w));
    g_alt = subs(g_alt,'s',1i*w);


%% Take only real part :
g1 = real(g1);
g2 = g2*1e5*exp(1i*w);
    
    