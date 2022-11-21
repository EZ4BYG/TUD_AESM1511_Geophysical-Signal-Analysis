# 3.1 Frequency-Wavenumber Filtering

In the previous assignment, we tried to filter surface waves in the frequency domain, i.e., you applied a one-dimensional filter. You also performed two-dimensional Fourier transformation from the time-space domain to the frequency-wavenumber (f-k) domain. Now, we are going to apply a two-dimensional filter in the f-k domain, i.e., an f-k filter.

# 3.2 Interpolation and Integration of Diffusive Fields

In exploration geophysics, wave and diffusive-field signals are computed in difficult geometries. That means that the computational costs might be very high, and we would like to minimize our computational efforts. Time-domain signals are transformed to the frequency domain and time-domain results are obtained again. Controlled-source electromagnetic fields are diffusive fields and behave as smooth functions in the space-time and in the space-frequency domains. For modelling, usually logarithmic frequency spacing is sufficient to capture the field. To transform the frequency-domain data back to the time domain, though, linear frequency spacing is necessary.

Suppose the function we need is known in a transformed domain and the spacefrequency domain function must be computed from the following integral:

$$
\hat{F}=i \omega \int_{\kappa=0}^{\infty} \frac{\exp (-\Gamma h)}{\Gamma} J_0(\rho \kappa) d \kappa
$$

where $\Gamma=\sqrt{\kappa^2+\gamma^2}$, $\gamma^2=i \omega \sigma \mu$, with $\omega=2 \pi f$ denoting angular frequency, $f$ - frequency, $\sigma$ - conductivity, and $mu$ - magnetic permeability; $h$ is the vertical distance between transmitter and receiver, and $\kappa$ is the integration variable. The function $J_0(\rho \kappa)$, with $\rho$ being the horizontal distance, is the Bessel function of order zero. We can compute this function using the command *besselj*.
