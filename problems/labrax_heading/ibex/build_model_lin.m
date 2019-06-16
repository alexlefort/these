function sym_model = build_model_lin(Vs)
    
    model = tbx_sm_load_model('iver2.nav');
    p = symtbx_convert_struct_to_sym(model,5);
    p.Vs = symtbx_convert_double_to_sym(Vs,5);
    
    ptmp = p;
    
    p.CyV   = sym('CyV');
    p.CnR   = sym('CnR');
    p.CyAL  = sym('CyAL');
    p.CnAL  = sym('CnAL');
    
    sys = tbx_sm_heading_model(p);
    
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
    sym_model.p     = p    ;
    sym_model.G     = G    ;
    sym_model.sys0  = sys0 ;
    sym_model.p0    = p0   ;
    sym_model.G0    = G0   ;
