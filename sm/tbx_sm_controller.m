function [command, intern] = tbx_sm_controller(measure, mission, pctrl, intern)

    [c_heading,intern.heading] = tbx_sm_heading_controller(measure, mission, intern.heading, pctrl.heading, pctrl.dt);
    [c_depth,intern.depth]     = tbx_sm_depth_controller  (measure, mission, intern.depth  , pctrl.depth  , pctrl.dt);

    command.Alpha = c_heading.Alpha ;
    command.Beta1 = c_depth.Beta1   ;
    command.Beta2 = c_depth.Beta2   ;
    command.Prop  = mission.Prop    ;