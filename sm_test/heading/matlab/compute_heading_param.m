function heading_param = compute_heading_param(param)

ifplot = false;
r = load_heading_user_param;

nVs = length(r.speeds) ;

% Vecteur des gains du controlleur
KP_optim_vect   = zeros(1,nVs);
KD_optim_vect   = zeros(1,nVs);
KPSI_optim_vect = zeros(1,nVs);

grid_regl = (1:r.NPTS)./r.NPTS;

for iiVs = 1:nVs

    %% Point de linearisation du modele
    Vs = r.speeds(iiVs); disp(Vs);
    MdB_R_min = r.MdB_R_min_vect(iiVs);

    [Xeq,Ueq] = calculer_equilibre(Vs, param, r);
    [A,B,~,~] = lineariser_modele(Xeq, Ueq, param, r.GIR);
    sysbo    = minreal(ss(A,B,[0 0 1 0 0],0)) ; % sorties = R
    sysboret = sysbo*r.ret                    ; % version du systeme avec retard pur

    %% REGLAGE PETITE BOUCLE

    KP_vect = (r.KP_max_vect(iiVs) - r.KP_min_vect(iiVs))*grid_regl + r.KP_min_vect(iiVs);
    KD_vect = (r.KD_max_vect(iiVs) - r.KD_min_vect(iiVs))*grid_regl + r.KD_min_vect(iiVs);

    OVERSH    =   1000*ones(r.NPTS) ;
    TEMPS_REP =   1000*ones(r.NPTS) ;
    MARGE     =  -1000*ones(r.NPTS) ;

    if (ifplot) figure; end
    for ii = 1:r.NPTS
        for jj = 1:r.NPTS

            KP = KP_vect(ii);
            KD = KD_vect(jj);

            controleur_PB = tf([KD KP],[1 0]);

            sysbf = feedback(sysboret*controleur_PB, 1) ;

            if (max(real(eig(sysbf))) < -1e-4) % si le systeme est stable en boucle fermee
                sysboret_corr = controleur_PB*sysboret;
                Mg = 20*log10(margin(sysboret_corr));

                if (Mg > MdB_R_min) % si la marge du systeme en boucle ouverte est suffisante

                    [Y,T] = step(sysbf,r.T_vect);
                    resultat_fin = abs(Y(end)-1);

                    if (resultat_fin <= 1e-2) % Si la cible est atteinte

                        % Calcul de l'oversh
                        oversh = max((max(Y)-1),0);

                        % Calcul du temps de reponse
                        z_max = size(Y,1);
                        z_rep = z_max;
                        while (abs(Y(z_rep,1)-1) <= 1e-2)
                            z_rep = z_rep - 1;
                        end

                        MARGE(ii,jj)     = Mg;
                        OVERSH(ii,jj)    = oversh;
                        TEMPS_REP(ii,jj) = T(z_rep);

                        if (ifplot) plot(T,Y,'g'); end
                        if (ifplot) hold on; end
                    end
                end
            end
        end
    end

    % Recherche du temps de reponse minimal
    val = min(min(OVERSH));

    for ii = 1:r.NPTS
        for jj = 1:r.NPTS
            if (OVERSH(ii,jj) > r.dep_val)
                TEMPS_REP(ii,jj) = 1000;
                OVERSH(ii,jj) = 1000;
            end
        end
    end

    [aux_1, ii_optim_tab] = min(TEMPS_REP);
    [val, jj_optim] = min(aux_1);
    ii_optim = ii_optim_tab(jj_optim);

    % Enregistrement des gains optimaux
    if (val ~= 1000)
        KP     = KP_vect  (ii_optim);
        KD     = KD_vect  (jj_optim);
        marge  = MARGE    (ii_optim,jj_optim);
        oversh = OVERSH   (ii_optim,jj_optim);
        trep   = TEMPS_REP(ii_optim,jj_optim);
        disp(['Vs = ',num2str(Vs),'; KP = ',num2str(KP),'; KD = ', num2str(KD)]);
        disp(['marge = ',num2str(marge),'; dep = ',num2str(oversh),'; trep = ', num2str(trep)]);
        KP_optim_vect(iiVs) = KP;
        KD_optim_vect(iiVs) = KD;
        ctrl_PB_optim = tf([KD KP],[1 0]);
        sysbf    = feedback(sysbo   *ctrl_PB_optim,1) ;
        sysbfret = feedback(sysboret*ctrl_PB_optim,1) ;
    else
        error('criteres pb trop restreints');
    end

    if (ifplot) hold off; end

    %% REGLAGE GRANDE BOUCLE

    KPSI_vect = r.KPSI_max_vect(iiVs)*grid_regl ;

    OVERSH    =   1000*ones(r.NPTS,1) ;
    TEMPS_REP =   1000*ones(r.NPTS,1) ;
    MARGE     =  -1000*ones(r.NPTS,1) ;

    sysgbbo = tf(1,[1 0])*sysbf ;
    sysgbboret = tf(1,[1 0])*sysbfret;

    % Marge de gain objectif

    MdB_Psi_min = r.MdB_Psi_min_vect(iiVs) ;

    if (ifplot) figure; end
    
    for ii = 1:r.NPTS

        KPSI = KPSI_vect(ii);
        sysbf = feedback(sysgbbo*KPSI,1) ;

        if (max(real(eig(sysbf))) < -1e-4) % si le systeme est stable en boucle fermee

            sysboret_corr = KPSI*sysgbboret;
            Mg = 20*log10(margin(sysboret_corr));

            if (Mg > MdB_Psi_min) % si la marge respecte l'objectif

                [Y,T] = step(sysbf,r.T_vect);
                resultat_fin = abs(Y(end)-1);

                if (resultat_fin <= 1e-2) % Si l'objectif est atteint

                    % Calcul de l'oversh
                    oversh = max((max(Y)-1),0);

                    % Calcul du temps de reponse
                    z_max = size(Y,1);
                    z_rep = z_max;
                    while (abs(Y(z_rep,1)-1) <= 1e-2)
                        z_rep = z_rep - 1;
                    end

                    MARGE(ii)     = Mg;
                    OVERSH(ii)    = oversh;
                    TEMPS_REP(ii) = T(z_rep);

                    if (ifplot) plot(T,Y); end
                    if (ifplot) hold on; end

                end
            end
        end
    end

    % Recherche du temps de reponse minimal
    for ii = 1:r.NPTS
        if (OVERSH(ii) > (r.dep_gb))
            TEMPS_REP(ii) = 1000;
            OVERSH(ii) = 1000;
        end
    end

    [val, ii_optim] = min(TEMPS_REP);

    % Enregistrement des gains optimaux
    if (val ~= 1000)
        KPSI   = KPSI_vect(ii_optim);
        marge  = MARGE    (ii_optim);
        oversh = OVERSH   (ii_optim);
        trep   = TEMPS_REP(ii_optim);
        disp(['Vs = ',num2str(Vs),'; KPSI = ',num2str(KPSI)]);
        disp(['marge = ',num2str(marge),'; dep = ',num2str(oversh),'; trep = ', num2str(trep)]);
        KPSI_optim_vect(iiVs) = KPSI;
    else
        error('criteres gb trop restreints')
    end
    if (ifplot) hold off; end
end

heading_param.gains_KR     = KD_optim_vect'   ;
heading_param.gains_KPSI   = KPSI_optim_vect' ;
heading_param.gains_IR     = KP_optim_vect'   ;
heading_param.gains_speeds = r.speeds'        ;
heading_param.sat_Alpha    = r.sat_Alpha      ;
heading_param.sat_Alpha_d  = r.sat_Alpha_d    ;
heading_param.sat_R_max    = r.sat_R_max      ;
heading_param.puls_Alpha   = 2*pi*5           ;
heading_param.omega        = 0.0              ;
heading_param.ksi          = 0.0              ;
heading_param.sat_dRho_min = 0.0              ;
heading_param.sat_dRho_max = 0.0              ;
heading_param.avance_phase = 0.0              ;

figure;
subplot(3,1,1);
plot(heading_param.gains_speeds, heading_param.gains_KR);
subplot(3,1,2);
plot(heading_param.gains_speeds, heading_param.gains_IR);
subplot(3,1,3);
plot(heading_param.gains_speeds, heading_param.gains_KPSI);