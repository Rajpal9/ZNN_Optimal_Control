% power sigmoid functions
function phi = pow_sig(e)
n = size(e,1);
% Parameters of the function
p = 3;
zeta = 4;
phi = zeros(n,1);
for i = 1:n
    if abs(e(i,1)) >= 1
        phi(i,1) = e(i,1)^(p);
    else
        phi(i,1) = ((1+exp(-zeta))/(1-exp(-zeta)))*((1-exp(-zeta*e(i,1)))/(1+exp(-zeta*e(i,1))));
    end
end