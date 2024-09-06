function [Sx, Ct,harden,pstr_p,Wb] = ElastoPlasticityBar(Ex,harden,pstr_p,el_mod,yield_bar,pl_mod)
pstr = real(sqrt(2*Ex+1));
pstr_e=pstr/pstr_p;
strain_tr = log(pstr_e);
stress_trial = el_mod * strain_tr;

trial_yield_f = abs(stress_trial) - ( yield_bar + pl_mod * harden );

if trial_yield_f <= 0
    tau = stress_trial;
    mat_stiff = el_mod;
    harden = 0;
elseif trial_yield_f > 0
    %% Inelastic part of return-mapping algorithm
    % Incremental plastic multiplier:
    incr_pl_mult = trial_yield_f / (el_mod + pl_mod);
    % fprintf("incr_pl_mult = %6.4f, pstr = %6.4f, harden = %6.4f\n", incr_pl_mult, pstr, harden);
    Dstrain_pl = incr_pl_mult * stress_trial / (abs(stress_trial));
    %stress = stress_trial - young * Dstrain_pl;
    harden = harden + incr_pl_mult;
    Dstrain_el=strain_tr-Dstrain_pl;
    tau=el_mod*Dstrain_el;
    pstr_p=pstr/(exp(Dstrain_el));
    % Algorithmic tangent modulus:
    mat_stiff = (el_mod * pl_mod) / (el_mod + pl_mod);
else
    error('/// Wrong analysis type value.')
end

Sx =tau/(pstr^2);
Ct=(mat_stiff-2*tau)/(pstr^4);
if nargout>2
    Wb = Sx*Ex/2;
end
end
