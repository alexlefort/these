clear all
clc
close all


Vs = 1.5;

%% Build linear model
	disp('--> build linear model');    
	model_sm = build_model_lin(Vs);

%% Build controller model
	disp('--> build controller model');    
	model_ctrl = build_controller(Vs);

%% Build transfert functions
	disp('--> build transferts'); 
	transferts = build_transferts(model_sm, model_ctrl);

%% Build criterias
	disp('--> build weightings'); 
	gabarits = build_weightings(model_sm);

%% Build criterias
	disp('--> build criterias'); 	
	criterias = build_criterias(transferts, gabarits);

    save('criterias.mat', 'criterias');
    save('transferts.mat', 'transferts');
    
%% Save in files

symtbx_save_criterion(criterias.psi1       , model_ctrl.gains , model_sm.p0 , 'Tpsi1.txt'      );
symtbx_save_criterion(criterias.psi2       , model_ctrl.gains , model_sm.p0 , 'Tpsi2.txt'      );
symtbx_save_criterion(criterias.psia       , model_ctrl.gains , model_sm.p0 , 'Tpsia.txt'      );
symtbx_save_criterion(criterias.stab_coefs , model_ctrl.gains , model_sm.p0 , 'Tstab_coefs.txt');
symtbx_save_criterion(criterias.stab_lc    , model_ctrl.gains , model_sm.p0 , 'Tstab_lc.txt'   );

symtbx_save_criterion(criterias.stab_coefs, model_ctrl.gains , model_sm.p0 , 'tstab.h');
save_lienard_chipart(length(criterias.stab_coefs), 'stab' , 'stab.h', sym(0));