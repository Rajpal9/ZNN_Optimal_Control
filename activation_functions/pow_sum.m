% power sum function
function phi = pow_sum(e)
n = size(e,1);
% Parameters of the function
N = 2;
phi = zeros(n,1);
for i = 1:n
    asum = 0;
    for k = 1:N
        asum = asum + e(i,1)^(2*k-1);
    end
    phi(i,1) = asum;
end