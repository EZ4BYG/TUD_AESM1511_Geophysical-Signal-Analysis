% Generate the integrand 
function [Integrand] = Integrand_Func(kappa, gamma, rho, h)
    Gamma = sqrt( kappa.^2 + gamma^2);
    J0 = besselj(0, rho.*kappa);
    Integrand = exp(-Gamma.*h) .* J0 ./ Gamma;
end