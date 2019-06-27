function criterias = build_criterias(transferts, gabarits)

    %% Build criterias separated for each variable

    criterias.zz1    = build_criteria_hinf(transferts.zz  / gabarits.zz1);
    criterias.zz2    = build_criteria_hinf(transferts.zz  / gabarits.zz2);
    criterias.zb     = build_criteria_hinf(transferts.zb  / gabarits.zb);
    
    %% Build stability criterion
    criterias.stab_coefs = build_criteria_stab_coeffs(transferts.poly_stab);