function fld = GBkf(k,grg,zdist,rho)
%GBkf: this function computes a product of Bessel functions I0K0 as an
%inverse spatial Fourier-Bessel integral
%grg=1i*omega*sigma*mu is gamma squared in equation 1 with
%sigma - the electric conductivity, omega=2*pi*f  - the angula frequency, and
% mu - the magnetic permeability.
%rho is the horizontal distance.
%zdist is the vertical distance between transmitter and receiver.

%Computing the integrand
fld=besselj(0,k.*rho).*exp(-zdist*sqrt(k.^2+grg))./sqrt(k.^2+grg);
end

