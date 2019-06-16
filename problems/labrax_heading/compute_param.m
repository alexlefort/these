model = tbx_sm_load_model('../sm/iver2.nav');

param_depth_controller   = ctrl_imm_calculer_reglages(model);
param_heading_controller = compute_heading_param(model);

save('param_depth_controller.mat'  , 'param_depth_controller'  );
save('param_heading_controller.mat', 'param_heading_controller');

