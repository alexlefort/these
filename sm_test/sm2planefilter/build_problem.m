clear all
clc
close all

%% Build linear model
	disp('--> build linear model');    
	model_sm = build_model_lin();

%% Build controller model
	disp('--> build controller model');    
	model_ctrl = build_controller();

%% Build filter model
	disp('--> build filter model');    
	model_filter = build_filter();

%% Build transfert functions
	disp('--> build transferts'); 
	transferts = build_transferts(model_sm, model_ctrl, model_filter);

%% Merge gains filter and gain controllers
        gains = symtbx_merge_structs(model_ctrl.gains, model_filter.gains);

%% Build criterias
	disp('--> build weightings'); 
	gabarits = build_weightings(model_sm);

%% Build criterias
	disp('--> build criterias'); 	
	criterias = build_criterias(transferts, gabarits);

	save('criterias.mat', 'criterias');

%% Save in files

symtbx_save_criterion(criterias.zz1        , gains , model_sm.p0 , 'functions/Tzz1.txt'       );
symtbx_save_criterion(criterias.zz2        , gains , model_sm.p0 , 'functions/Tzz2.txt'       );
symtbx_save_criterion(criterias.zb1        , gains , model_sm.p0 , 'functions/Tzb1.txt'       );
symtbx_save_criterion(criterias.zb2        , gains , model_sm.p0 , 'functions/Tzb2.txt'       );
symtbx_save_criterion(criterias.stab_coefs , gains , model_sm.p0 , 'functions/Tstab_coefs.txt');
%%symtbx_save_criterion(criterias.stab_lc    , gains , model_sm.p0 , 'functions/Tstab_lc.txt'   );

symtbx_save_criterion(criterias.stab_coefs, model_ctrl.gains , model_sm.p0 , 'tstab.h');
save_lienard_chipart(length(criterias.stab_coefs), 'stab' , 'stab.h', sym(0));
