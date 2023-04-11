    % Runge-Kutta scheme

    k1 = delta_t * lorenz63(x, sigma, rho, beta);
    k2 = delta_t * lorenz63(x + k1/2, sigma, rho, beta);
    k3 = delta_t * lorenz63(x + k2/2, sigma, rho, beta);
    k4 = delta_t * lorenz63(x + k3, sigma, rho, beta);
    x = x + (k1 + 2*k2 + 2*k3 + k4)/6;   