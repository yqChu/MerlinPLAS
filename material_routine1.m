function [stress,mat_stiff,pstr_p,harden] = ...
    material_routine1(strain_tr,pstr,harden,pstr_p,iter)
%% MATERIAL_ROUTINE Returns stress and material stiffness.
% Inputs:
% params - structure containing material parameters.
% strain - total strain.
% analysis_type - string ('elastic' or 'plastic') to decide if plasticity is allowed.
% el - element number.

% Output:
% stress - Kirchhoff stress.
% mat_stiffness - derivative of stress w.r.t. strain.
% params - updated input structure (Matlab sends a copy).

% Young's modulus, [N/mm^2] (see p.89).
% if crease_flag == 0
%     young = 11001;
% else
%     young = 11000;
% end
young = 100;
yield = 8000000; % Yield stress, [N/mm^2].
plast_mod = 2000; % Plastic modulus H, [N/mm^2] (smth small compared to Young's modulus).

stress_trial = young * strain_tr;
%stress_trial = young * (strain - state.plast_strain);
trial_yield_f = abs(stress_trial) - ( yield + plast_mod * harden );

% trial_yield_f
if trial_yield_f <= 0
    stress = stress_trial;
    mat_stiff = young;
    harden = 0;
elseif trial_yield_f > 0
    % fprintf("Now there is plasticity!n");
    % aa = 1;
    % if iter<10
    %     aa=1;
    % else
    %     aa=aa/2;
    % end
    %% Inelastic part of return-mapping algorithm
    % Incremental plastic multiplier:
    incr_pl_mult = trial_yield_f / (young + plast_mod);
    % fprintf("incr_pl_mult = %6.4f, pstr = %6.4f, harden = %6.4f\n", incr_pl_mult, pstr, harden);
    Dstrain_pl = incr_pl_mult * stress_trial / (abs(stress_trial));
    %stress = stress_trial - young * Dstrain_pl;
 %   state.plast_strain = state.plast_strain + Dstrain_pl;
    harden = harden + incr_pl_mult;
    Dstrain_el=strain_tr-Dstrain_pl;
    stress=young*Dstrain_el;
    % stress = stress_trial - young*Dstrain_pl;
    pstr_p=pstr/(exp(Dstrain_el));
    % Algorithmic tangent modulus:
    mat_stiff = (young * plast_mod) / (young + plast_mod);
else
    error('/// Wrong analysis type value.')
end

end