% hyperbolic sine function
function phi = hyp_sin(e)
n = size(e,1);
% Parameters of the function
zeta = 2;
phi = zeros(n,1);
for i = 1:n
    
        phi(i,1) = (exp(zeta*e(i,1))-exp(-zeta*e(i,1)))/2;
  
end