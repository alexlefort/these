function r = load_heading_user_param

deg = pi/180;

r.speeds           = [  0.3   0.5   1.0   1.5   2.0    2.5   3.0    4.0   5.0 ] ;
r.KP_max_vect      = [ 10.0  10.0   8.0   6.0   5.0    4.0   3.0    2.0   1.0 ] ;
r.KP_min_vect      = [  0.0   0.0   0.0   0.0   0.0    0.0   0.0    0.0   0.0 ] ;
r.KD_max_vect      = [ 40.0  20.0  12.0   8.0   4.5    3.5   2.0    1.0   1.0 ] ;
r.KD_min_vect      = [  0.0   0.0   0.0   0.0   0.0    0.0   0.0    0.0   0.0 ] ;
r.KPSI_max_vect    = [  1.0   1.0   1.0   1.0   1.0    1.0   1.0    0.5   0.5 ] ;
r.MdB_R_min_vect   = [ 12.0  11.0  10.0  10.0  10.0   10.0  10.0   10.0  10.0 ] ;
r.MdB_Psi_min_vect = [  5.0   5.0   5.0   5.0   5.0    5.0   5.0    5.0   5.0 ] ;
r.Trep_vect        = [ 10.0  10.0  10.0  10.0  10.0   10.0  10.0   10.0  10.0 ] ;

r.GIR.states = {'V';'P';'R';'Phi';'Psi'};
r.GIR.inputs = {'Alpha'};
r.GIR.mes    = {'V';'P';'R';'Phi';'Psi'};

r.NPTS    = 30            ; % nb de points de reglage evalues
r.T_MAX   = 30            ; % Temps de reponse maximal tolere
r.T_vect  = 0:0.1:r.T_MAX ;
r.dep_val = 0.10          ; % 5% du Minimum
r.dep_gb  = 0.05          ;

[num_pade,den_pade] = pade(0.5,4);
r.ret = tf(num_pade,den_pade);

r.sat_Alpha      = 20*deg ;
r.sat_Alpha_d    = 12*deg ;
r.sat_R_max      = 3*deg  ;