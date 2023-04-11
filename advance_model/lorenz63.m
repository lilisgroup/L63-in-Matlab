% Lorenz '63 systems
function lorenz = lorenz63(x, sigma, rho, beta)


x_dot  = -sigma*x(1) + sigma*x(2);
y_dot  = -x(1)*x(3) + rho*x(1) - x(2);
z_dot  =  x(1)*x(2) - beta*x(3);


lorenz = [x_dot y_dot z_dot]';