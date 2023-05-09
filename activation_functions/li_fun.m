% li function
function phi = li_fun(e)
n = size(e,1);
% Parameters of the function
r = 0.5;
phi = zeros(n,1);
for i = 1:n
    
    phi(i,1) = ((sign(e(i,1))*((abs(e(i,1)))^r)) + (sign(e(i,1))*((abs(e(i,1)))^(1/r))))/2;
end