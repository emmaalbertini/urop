% master RHS
function dydt = master_RHS(t, y, p, bioreactor)

    y_act = convert_state_vec(y, bioreactor);
    
    Ecell_idx = getStateIdx('Ecell', bioreactor.Ecell);
    Mrcell_idx = getStateIdx('Mrcell', bioreactor.Mrcell);
    Mpcell_idx = getStateIdx('Mpcell', bioreactor.Mpcell);
    Mcell_idx = getStateIdx('Mcell', bioreactor.Mcell);
    
    buffer = p.N - (y_act.Ecell(Ecell_idx) + y_act.Mpcell(Mpcell_idx)+ y_act.Mrcell(Mrcell_idx)+ y_act.Mcell(Mcell_idx));

    dydt_Ecell = E_Cell_RHS(t, y_act, p, bioreactor, buffer);
    dydt_Mrcell = Mr_Cell_RHS(t, y_act, p, bioreactor, buffer);
    dydt_Mpcell = Mp_Cell_RHS(t, y_act, p, bioreactor, buffer);
    dydt_Mcell = M_Cell_RHS(t, y_act, p, bioreactor, buffer);    
    
    a = 555;
    
    dydt = [dydt_Ecell; dydt_Mrcell; dydt_Mpcell; dydt_Mcell];
    
end






