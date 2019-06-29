%% Find cutoff pulsation for a low pass symbolic tf 

function cut_off = symtbx_find_cutoff_pulsation(bode_t,w)

    cut_off  = interp1(bode_t, w,-3);