function [Rspr, rot_stiff, harden_fold, he_pl] = ElastoPlasticityHinge(he, h0, kpi, L0, harden_fold, he_pl, yield_fold, pl_mod_fold,del, icrm)
he_tr = he-he_pl+h0;
Rspr_tr = L0.*kpi*(he_tr - h0);

trial_yield_f = abs(Rspr_tr) - L0*(yield_fold + pl_mod_fold*harden_fold);

if trial_yield_f <=0
    Rspr = Rspr_tr;
    rot_stiff = L0.*kpi;
    harden_fold = 0;
elseif trial_yield_f > 0
    incr_pl_mult = 1/L0*trial_yield_f/(kpi + pl_mod_fold);
    he_el = he_tr - incr_pl_mult*Rspr_tr/abs(Rspr_tr);
    harden_fold = harden_fold + incr_pl_mult;
    Rspr = L0.*kpi*(he_el - h0);
    rot_stiff = L0.*kpi*pl_mod_fold/(kpi + pl_mod_fold);
    he_pl = he - he_el + h0;
else
    error('/// Wrong analysis type value.')
end
