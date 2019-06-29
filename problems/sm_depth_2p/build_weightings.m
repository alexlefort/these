function g = build_weightings(model_sys)

	delta_w = 1.0 ; %% Pulsation margin
    
    deg = sym(10/572);

	syms s;
	logw = -3:0.01:1;
	w    = 10.^(logw); %% Pulsation vector
    
    %% Todo make plot function
    b1z = symtbx_symtf_to_tf(model_sys.G0(1,1));
    b1t = symtbx_symtf_to_tf(model_sys.G0(1,2));
    b2z = symtbx_symtf_to_tf(model_sys.G0(2,1));
    b2t = symtbx_symtf_to_tf(model_sys.G0(2,2));
    
    %% state space representation
    [b1z_mag, ~, ~] = bode(b1z,w);  
    [b1t_mag, ~, ~] = bode(b1t,w);  
    [b2z_mag, ~, ~] = bode(b2z,w);  
    [b2t_mag, ~, ~] = bode(b2t,w);  

    b1z0  = squeeze(b1z_mag);
    b1t0  = squeeze(b1t_mag); 
    b2z0  = squeeze(b2z_mag); 
    b2t0  = squeeze(b2t_mag);     
    
	%% Build weightings
    g.b   = sym(20)*deg;
    g.bdt = sym(4)*deg/s;
    g.t   = deg;
    g.z   = sym(11/10);

    wc = symtbx_find_cutoff_pulsation(20*log10(b1z0),w);
    disp(wc);
    wd = 10^(log10(wc)-1);
    disp(1/wd);
	g.zz1 = s*sym(562/10); 
    

    
    figure; loglog(w,b1z0,'r',w,b1t0,'b',w,b2z0,'k',w,b2t0,'p');
    grid on;
