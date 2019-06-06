function criterias = build_criterias(transferts, gabarits)

    %% Build criterias separated for each variable
    nvpa = 5;
    criterias.zz1    = build_criteria_hinf(transferts.zz  / gabarits.zz1);
    disp('c1');
    criterias.zz2    = build_criteria_hinf(transferts.zz  / gabarits.zz2);
    criterias.zb1    = build_criteria_hinf(transferts.zb1 / gabarits.zb1);
    criterias.zb2    = build_criteria_hinf(transferts.zb2 / gabarits.zb2);
    criterias.zz_t   = build_criteria_hinf(transferts.zz);
    criterias.zb1_t  = build_criteria_hinf(transferts.zb1);
    criterias.zb2_t  = build_criteria_hinf(transferts.zb1);
    disp('c2');
    criterias.zz1_w  = build_criteria_hinf(gabarits.zz1 );
    criterias.zz2_w  = build_criteria_hinf(gabarits.zz2 );
    criterias.zb1_w  = build_criteria_hinf(gabarits.zb1 );
    criterias.zb2_w  = build_criteria_hinf(gabarits.zb2 );
    %% Build stability criterion

    criterias.stab_coefs = build_criteria_stab_coeffs(transferts.poly_stab);
    %%criterias.stab_lc    = build_criteria_stab_Lienard_Chipart(transferts.poly_stab);
    
    disp('c3');
    for ii=1:length(criterias.stab_coefs)
        criterias.stab_coefs{ii} = criterias.stab_coefs{ii};
    end