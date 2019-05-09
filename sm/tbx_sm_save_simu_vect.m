function res = tbx_sm_save_simu_vect(res, state, measure, command, mission, ii)
    
    f = fields(res.state);
    for nn=1:length(f)
        nom = f{ii};
        res.state.(nom)(ii) = state.(nom);
    end

    f = fields(res.measure);
    for nn=1:length(f)
        nom = f{ii};
        res.measure.(nom)(ii) = measure.(nom);
    end
    
    f = fields(res.command);
    for nn=1:length(f)
        nom = f{ii};
        res.command.(nom)(ii) = command.(nom);
    end
    
    f = fields(res.mission);
    for nn=1:length(f)
        nom = f{ii};
        res.mission.(nom)(ii) = mission.(nom);
    end
    