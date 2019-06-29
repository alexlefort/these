function criterias = build_criterias(t, g)

    %% Build criterias separated for each variable
    nvpa = 5;
    criterias.zzc  = simplify(vpa(build_criteria_hinf(t.zz  / g.z  ),nvpa));
    criterias.zzs  = simplify(vpa(build_criteria_hinf(t.zz  / g.zz1),nvpa));
    criterias.ttc  = simplify(vpa(build_criteria_hinf(t.tt  / g.z  ),nvpa)); 
    
    criterias.zb1c = simplify(vpa(build_criteria_hinf(t.zb1 / g.b),nvpa));
    criterias.zb2c = simplify(vpa(build_criteria_hinf(t.zb2 / g.b),nvpa));
    criterias.zb1s = simplify(vpa(build_criteria_hinf(t.zb1 / g.bdt),nvpa));
    criterias.zb2s = simplify(vpa(build_criteria_hinf(t.zb2 / g.bdt),nvpa));
    
    criterias.tb1c  = simplify(vpa(build_criteria_hinf(t.zb1 / g.b*g.t),nvpa));
    criterias.tb2c  = simplify(vpa(build_criteria_hinf(t.zb2 / g.b*g.t),nvpa));
    criterias.tb1s  = simplify(vpa(build_criteria_hinf(t.zb1 / g.bdt*g.t),nvpa));
    criterias.tb2s  = simplify(vpa(build_criteria_hinf(t.zb2 / g.bdt*g.t),nvpa));
    
    %% Build stability criterion

    criterias.stab_coefs = build_criteria_stab_coeffs(t.poly_stab);
    criterias.stab_lc    = build_criteria_stab_Lienard_Chipart(t.poly_stab);
    
    for ii=1:length(criterias.stab_coefs)
        criterias.stab_coefs{ii} = simplify(vpa(criterias.stab_coefs{ii},nvpa));
    end