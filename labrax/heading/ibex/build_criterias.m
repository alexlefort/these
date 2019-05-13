function criterias = build_criterias(transferts, gabarits)

    %% Build criterias separated for each variable

    criterias.psi1    = build_criteria_hinf(transferts.psi  / gabarits.psi1);
    criterias.psi2    = build_criteria_hinf(transferts.psi  / gabarits.psi2);
    criterias.psia    = build_criteria_hinf(transferts.psia / gabarits.psia);
    criterias.psi_t   = build_criteria_hinf(transferts.psi);
    criterias.psia_t  = build_criteria_hinf(transferts.psia);
    criterias.psi1_w  = build_criteria_hinf(gabarits.psi1 );
    criterias.psi2_w  = build_criteria_hinf(gabarits.psi2 );
    criterias.psia_w  = build_criteria_hinf(gabarits.psia );
    %% Build stability criterion

    criterias.stab_coefs = build_criteria_stab_coeffs(transferts.poly_stab);
    criterias.stab_lc    = build_criteria_stab_Lienard_Chipart(transferts.poly_stab);