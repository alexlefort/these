function gabarits = build_weightings(model_sys)

	delta_w = 1.0 ; %% Pulsation margin

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

	%% Build weightings

	gabarits.zz1 = s*sym(10); 
	gabarits.zz2 = sym(16/10);
	gabarits.zb1 = sym(16/10);
    gabarits.zb2 = sym(16/10);

	tf_wzz1 = symtbx_symtf_to_tf(gabarits.zz1);
    tf_wzz2 = symtbx_symtf_to_tf(gabarits.zz2);
    tf_wzb1 = symtbx_symtf_to_tf(gabarits.zb1);
    tf_wzb2 = symtbx_symtf_to_tf(gabarits.zb2);

    [wzz1_mag, ~, ~] = bode(tf_wzz1, w);  
    [wzz2_mag, ~, ~] = bode(tf_wzz2, w);  
    [wzb1_mag, ~, ~] = bode(tf_wzb1, w); 
    [wzb2_mag, ~, ~] = bode(tf_wzb2, w); 
    
    b1z0  = squeeze(b1z_mag);
    b1t0  = squeeze(b1t_mag); 
    b2z0  = squeeze(b2z_mag); 
    b2t0  = squeeze(b2t_mag); 
    
    wzz1 = squeeze(wzz1_mag);
    wzz2 = squeeze(wzz2_mag);
    wzb1 = squeeze(wzb1_mag);
    wzb2 = squeeze(wzb2_mag);

    %%figure; loglog(w, zz0,'r',w , wzz1,'b', w, wzz2, 'b');
    %%grid on;
    
    figure; loglog(w,b1z0,'r',w,b1t0,'b',w,b2z0,'k',w,b2t0,'p');
    grid on;
