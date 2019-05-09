function p = reglage_ctrl_imm()


nav = charger_modele_sm('iver2.nav');

    s = tf('s');

deg = pi/180;
    
p.UVect   = [0.3 0.5 1.0 1.5 2.0 2.5 3.0 4.0 5.0];
n         = length(p.UVect);

p.beta_max   = 20*deg;
p.dbeta_max  = 15*deg;

p.IZVect  = zeros(n,1);
p.KZVect  = zeros(n,1);
p.KDZVect = zeros(n,1);
p.KTVect  = zeros(n,1);
p.KDTVect = zeros(n,1);

p.DZ_max  = 10.0;
p.omega_z = 0.2*pi;
p.omega_t = 0.2*pi;

for ii=1:length(p.UVect)
    
	Vs = p.UVect(ii);
    
	%% Build linear model 
	
	    sm = build_model_imm(nav, Vs);   
	
    %% Build controller model
	
	    g.kpz  = realp('kpz' ,0);
        g.kdz  = realp('kdz' ,0);
        g.kpt  = realp('kpt' ,0);
        g.kdt  = realp('kdt' ,0);
        g.kiz  = realp('kiz' ,0);
        g.kiz.Maximum =  0.5;
        g.kiz.Minimum = -0.5;
        g.kpz.Maximum =  0.5;
        g.kpz.Minimum = -0.5;          
        ctrl            = [(g.kpz+g.kiz/s) g.kpt g.kdz g.kdt];
        ctrl.InputName  = {'eZf','eThetaf','eWf','eQf'};
        ctrl.OutputName = {'Beta'};
	
    %% Build filter model    

        oz = p.omega_z;
        ot = p.omega_t;
        denumz = 1 + sqrt(2)/oz*s + s^2/(oz^2);
        denumt = 1 + sqrt(2)/ot*s + s^2/(ot^2);
        
        filter = [1/denumz 0        ; ...
                  s/denumz 0        ; ...
                  0        1/denumt ; ...              
                  0        s/denumt];
        filter.InputName  = {'Z','Theta'};
        filter.OutputName = {'Zf','Thetaf','Wf','Qf'};
        sm0 = sm.NominalValue;
        fc = 0.01*find_cutoff_frequency(sm0(1,1));
              
        wz = mk_pondw0(-40,fc/(2*pi),2,4,'wz');

        wz.InputName  = 'eZ';
        wz.OutputName = 'wZ';
        
        wb = tf(5*pi/180,1);
        wb.InputName  = 'Beta' ;
        wb.OutputName = 'wBeta';

    %% Build transfert functions
    	
        ref       = sumblk('eZ      = rZ  - Z'     );
        refZf     = sumblk('eZf     = rZ  - Zf'    );
        refThetaf = sumblk('eThetaf =     - Thetaf');
        refWf     = sumblk('eWf     =     - Wf'    );
        refQf     = sumblk('eQf     =     - Qf'    );
        
        reff = connect(ref, refZf, refThetaf, refWf, refQf,  {'rZ', 'Z', 'Zf', 'Thetaf', 'Wf', 'Qf'}, {'eZ', 'eZf', 'eThetaf', 'eWf', 'eQf'});
        figure; hold on; bode(sm(1,1)) ; bode(wz); grid on;
        sys = connect(sm,filter,ctrl,reff,{'rZ'}, {'eZ', 'Beta'});
    
        perf = TuningGoal.WeightedGain({'rZ'}, {'eZ'}, [], 1/wz);
        perf.focus = [1e-5 10];
    
        sat  = TuningGoal.WeightedGain({'rZ'}, {'Beta'}, [], 1/wb);
        sat.focus  = [1e-5 10];
    
        opt = systuneOptions('RandomStart', 10);
        [CL,fSoft] = systune(sys, [perf,sat], opt);
        
        sysZ = connect(CL,{'rZ'},{'eZ'});
        sysB = connect(CL,{'rZ'},{'Beta'});
        
        res = getBlockValue(CL);
        
        p.IZVect(ii)  = res.kiz;
        p.KZVect(ii)  = res.kpz;
        p.KDZVect(ii) = res.kdz;
        p.KTVect(ii)  = res.kpt;
        p.KDTVect(ii) = res.kdt;
    
        figure;  subplot(1,2,1);  stepplot(sysZ);
                 subplot(1,2,2);  stepplot(sysB);
        figure;  hold on; bode(sysZ); bode(wz);
        figure;  hold on; bode(sysB); bode(wb);
        showTunable(CL);

end
