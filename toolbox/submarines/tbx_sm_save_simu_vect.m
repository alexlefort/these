function res = tbx_sm_save_simu_vect(res, state, measure, command, order, ii)
    
    f = fields(res.state);

    for nn=1:length(f)
        nom = f{nn};       
        res.state.(nom)(ii) = state.(nom);
    end

    f = fields(res.measure);
    for nn=1:length(f)
        nom = f{nn};
        res.measure.(nom)(ii) = measure.(nom);
    end
    
    f = fields(res.command);
    for nn=1:length(f)
        nom = f{nn};
        res.command.(nom)(ii) = command.(nom);
    end
    
    f = fields(res.order);
    for nn=1:length(f)
        nom = f{nn};
        res.order.(nom)(ii) = order.(nom);
    end
    