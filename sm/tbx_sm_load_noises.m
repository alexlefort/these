function noise = tbx_sm_load_noises()

deg = pi/180;

%% Modele de bruit utilise pour les reglages des estimateurs et des simulateurs

noise.Phi    = 0.0138;
noise.Theta  = 0.0125;
noise.Psi    = 0.0248;
noise.P      = 0.0015; 
noise.Q      = 0.0010;
noise.R      = 0.0010;
noise.U      = 0.2265;
noise.V      = 0.2268;
noise.W      = 0.2268;
noise.Z      = 0.06  ;
