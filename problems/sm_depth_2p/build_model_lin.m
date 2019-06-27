function sym_model = build_model_lin(model, Vs)
    
    prec = 5;
    p = symtbx_convert_struct_to_sym(model,prec);
    p.Vs = symtbx_convert_double_to_sym(Vs,prec);
    
    ptmp = p;
    p.CzB1 = sym('CzB1');
    p.CmB1 = sym('CmB1');
    p.CzB2 = sym('CzB2');
    p.CmB2 = sym('CmB2');
    p.CzW  = sym('CzW' );
    p.CmQ  = sym('CmQ' );
    
    sys = tbx_sm_model_depth_2p(p);
    
	%% Build central model
    
    var_po = [symvar(sys.a) symvar(sys.b)];
    p0 = [];
    for ii=1:length(var_po)
        p0.(char(var_po(ii))) = ptmp.(char(var_po(ii)));
    end
    
	sys0.a = eval(symtbx_sym_subs_from_struct(sys.a, p0));
	sys0.b = eval(symtbx_sym_subs_from_struct(sys.b, p0));
	sys0.c = eval(symtbx_sym_subs_from_struct(sys.c, p0));
	sys0.d = eval(symtbx_sym_subs_from_struct(sys.d, p0));
    
	syms s;
	
	G  = simplify(symtbx_symss_to_symtf(sys));
	G0 = simplify(symtbx_symss_to_symtf(sys0));

	%% Save in struct
	sym_model.sys   = sys  ;
	sym_model.G     = G    ;
    sym_model.sys0  = sys0 ;
    sym_model.p0    = p0   ;
	sym_model.G0    = G0   ;