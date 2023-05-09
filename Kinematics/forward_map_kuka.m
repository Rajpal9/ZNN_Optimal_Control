function X = forward_map_kuka(a8,d1,d3,d5,theta)
    theta1 = theta(1);
    theta2 = theta(2);
    theta3 = theta(3);
    theta4 = theta(4);
    theta5 = theta(5);
    theta6 = theta(6);
    theta7 = theta(7);
    
    c1 = cos(theta1);
    c2 = cos(theta2);
    c3 = cos(theta3);
    c4 = cos(theta4);
    c5 = cos(theta5);
    c6 = cos(theta6);
    c7 = cos(theta7);
    
    s1 = sin(theta1);
    s2 = sin(theta2);
    s3 = sin(theta3);
    s4 = sin(theta4);
    s5 = sin(theta5);
    s6 = sin(theta6);
    s7 = sin(theta7);
    
    sigma2x = s1*s3-c1*c2*c3;
    sigma1x = s4*sigma2x+c1*c4*s2;
    x = a8*(c6*sigma1x-s6*(c5*(c4*sigma2x-c1*s2*s4)+s5*(c3*s1+c1*c2*s3)))+d5*sigma1x+d3*c1*s2;
    
    sigma2y = c1*s3+c2*c3*s1;
    sigma1y = s4*sigma2y-c4*s1*s2;
    y = d3*s1*s2-d5*sigma1y-a8*(c6*sigma1y-s6*(c5*(c4*sigma2y+s1*s2*s4)+s5*(c1*c3-c2*s1*s3)));
    
    sigma1z = c2*c4+c3*s2*s4;
    z = d1+a8*(s6*(c5*(c2*s4-c3*c4*s2)+s2*s3*s5)+c6*sigma1z)+d5*sigma1z+d3*c2;
    
    X = [x;y;z];
end