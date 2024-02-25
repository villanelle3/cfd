clear; clc;

%---------------------------- Condição Inicial ----------------------------
[u, v, nu_til, nu_t] = Canal.AplicarCondicaoInicial();
% -------------------------------------------------------------------------

% ----------------- Loop para avançar no tempo usando RK4 -----------------
dt = AvancoTemporal.dt;