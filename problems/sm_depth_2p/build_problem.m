clear all
clc
close all

Vs =  2.5;

model = tbx_sm_load_model('../models/sm_test.nav');

%% Build linear model
	disp('--> build linear model');    
	model_sm = build_model_lin(model,Vs);

%% Build controller model
	disp('--> build controller model');    
	model_ctrl = build_controller();

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

%% Save in files

symtbx_save_criterion(criterias.zzc  , model_ctrl.gains , model_sm.p0 , 'functions/Tzzc.txt');
symtbx_save_criterion(criterias.zzs  , model_ctrl.gains , model_sm.p0 , 'functions/Tzzs.txt');
symtbx_save_criterion(criterias.ttc  , model_ctrl.gains , model_sm.p0 , 'functions/Tttc.txt');
    
symtbx_save_criterion(criterias.zb1c , model_ctrl.gains , model_sm.p0 , 'functions/Tzb1c.txt');
symtbx_save_criterion(criterias.zb2c , model_ctrl.gains , model_sm.p0 , 'functions/Tzb2c.txt');
%symtbx_save_criterion(criterias.zb1s , model_ctrl.gains , model_sm.p0 , 'functions/Tzb1s.txt');
%symtbx_save_criterion(criterias.zb2s , model_ctrl.gains , model_sm.p0 , 'functions/Tzb2s.txt');
    
symtbx_save_criterion(criterias.tb1c , model_ctrl.gains , model_sm.p0 , 'functions/Ttb1c.txt');
symtbx_save_criterion(criterias.tb2c , model_ctrl.gains , model_sm.p0 , 'functions/Ttb2c.txt');
%symtbx_save_criterion(criterias.tb1s , model_ctrl.gains , model_sm.p0 , 'functions/Ttb1s.txt');
%symtbx_save_criterion(criterias.tb2s , model_ctrl.gains , model_sm.p0 , 'functions/Ttb2s.txt');

symtbx_save_criterion(criterias.stab_coefs , model_ctrl.gains , model_sm.p0 , 'functions/Tstab_coefs.txt');
symtbx_save_criterion(criterias.stab_lc    , model_ctrl.gains , model_sm.p0 , 'functions/Tstab_lc.txt'   );

%symtbx_save_criterion(criterias.stab_coefs, model_ctrl.gains , model_sm.p0 , 'tstab.h');
save_lienard_chipart(length(criterias.stab_coefs), 'stab' , 'stab.h', sym(0));