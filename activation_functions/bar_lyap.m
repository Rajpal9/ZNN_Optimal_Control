% barrier lyapanov functions
 function [phi_err,k2] = bar_lyap(e,t,e0)
% e = error
% i =  adaptivity counter
%e0 = initial error
n = size(e,1);
s = 3*t+1;
%s = exp(t);
% Parameters of the function
k2 = 4*e0/s;
phi_err = zeros(n,1);
for i = 1:n
    phi_err(i,1) = (k2*e(i,1))./(k2^2 - e(i,1)^(2));
end