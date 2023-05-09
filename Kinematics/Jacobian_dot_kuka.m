function J_dot = Jacobian_dot_kuka(a8,d1,d3,d5,theta,theta_dot)
    theta1 = theta(1);
    theta2 = theta(2);
    theta3 = theta(3);
    theta4 = theta(4);
    theta5 = theta(5);
    theta6 = theta(6);
    theta7 = theta(7);
    
    theta1_dot = theta_dot(1);
    theta2_dot = theta_dot(2);
    theta3_dot = theta_dot(3);
    theta4_dot = theta_dot(4);
    theta5_dot = theta_dot(5);
    theta6_dot = theta_dot(6);
    theta7_dot = theta_dot(7);

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
    
    %time derivatives
    c1_dot = -sin(theta1)*theta1_dot;
    c2_dot = -sin(theta2)*theta2_dot;
    c3_dot = -sin(theta3)*theta3_dot;
    c4_dot = -sin(theta4)*theta4_dot;
    c5_dot = -sin(theta5)*theta5_dot;
    c6_dot = -sin(theta6)*theta6_dot;
    c7_dot = -sin(theta7)*theta7_dot;
    
    s1_dot = cos(theta1)*theta1_dot;
    s2_dot = cos(theta2)*theta2_dot;
    s3_dot = cos(theta3)*theta3_dot;
    s4_dot = cos(theta4)*theta4_dot;
    s5_dot = cos(theta5)*theta5_dot;
    s6_dot = cos(theta6)*theta6_dot;
    s7_dot = cos(theta7)*theta7_dot;
    
    sigma15 = c1*s3+c2*c3*s1;
    sigma14 = s1*s3-c1*c2*c3;
    sigma13 = c1*c3-c2*s1*s3;
    sigma12 = c4*sigma15+s1*s2*s4;
    sigma11 = c3*s1+c1*c2*s3;
    sigma10 = c4*sigma14-c1*s2*s4;
    sigma9 = c5*sigma12+s5*sigma13;
    sigma8 = c5*sigma10+s5*sigma11;
    sigma7 = s4*sigma14+c1*c4*s2;
    sigma6 = s4*sigma15-c4*s1*s2;
    sigma5 = c2*c4+c3*s2*s4;
    sigma4 = c4*s2-c2*c3*s4;
    sigma3 = c1*c2*c4+c1*c3*s2*s4;
    sigma2 = c2*c4*s1+c3*s1*s2*s4;
    sigma1 = c2*s4-c3*c4*s2;
    
    %derivatives
    sigma15_dot = c1_dot*s3+c1*s3_dot+c2_dot*c3*s1+c2*c3_dot*s1+c2*c3*s1_dot;
    sigma14_dot = s1_dot*s3+s1*s3_dot-(c1_dot*c2*c3+c1*c2_dot*c3+c1*c2*c3_dot);
    sigma13_dot = c1_dot*c3+c1*c3_dot-(c2_dot*s1*s3+c2*s1_dot*s3+c2*s1*s3_dot);
    sigma12_dot = c4_dot*sigma15+c4*sigma15_dot+s1_dot*s2*s4+s1*s2_dot*s4+s1*s2*s4_dot;
    sigma11_dot = c3_dot*s1+c3*s1_dot+(c1_dot*c2*s3+c1*c2_dot*s3+c1*c2*s3_dot);
    sigma10_dot = c4_dot*sigma14+c4*sigma14_dot-(c1_dot*s2*s4+c1*s2_dot*s4+c1*s2*s4_dot);
    sigma9_dot = c5_dot*sigma12+c5*sigma12_dot+s5_dot*sigma13+s5*sigma13_dot;
    sigma8_dot = c5_dot*sigma10+c5*sigma10_dot+s5_dot*sigma11+s5*sigma11_dot;
    sigma7_dot = s4_dot*sigma14+s4*sigma14_dot+c1_dot*c4*s2+c1*c4_dot*s2+c1*c4*s2_dot;
    sigma6_dot = s4_dot*sigma15+s4*sigma15_dot-(c4_dot*s1*s2+c4*s1_dot*s2+c4*s1*s2_dot);
    sigma5_dot = c2_dot*c4+c2*c4_dot+c3_dot*s2*s4+c3*s2_dot*s4+c3*s2*s4_dot;
    sigma4_dot = c4_dot*s2+c4*s2_dot-(c2_dot*c3*s4+c2*c3_dot*s4+c2*c3*s4_dot);
    sigma3_dot = c1_dot*c2*c4+c1*c2_dot*c4+c1*c2*c4_dot+c1_dot*c3*s2*s4+c1*c3_dot*s2*s4+c1*c3*s2_dot*s4+c1*c3*s2*s4_dot;
    sigma2_dot = c2_dot*c4*s1+c2*c4_dot*s1+c2*c4*s1_dot+c3_dot*s1*s2*s4+c3*s1_dot*s2*s4+c3*s1*s2_dot*s4+c3*s1*s2*s4_dot;
    sigma1_dot =c2_dot*s4+c2*s4_dot-(c3_dot*c4*s2+c3*c4_dot*s2+c3*c4*s2_dot);
  
    j11_dot = a8*(c6_dot*sigma6+c6*sigma6_dot-(s6_dot*sigma9+s6*sigma9_dot))+d5*sigma6_dot-d3*(s1_dot*s2+s1*s2_dot);
    j12_dot = a8*(c6_dot*sigma3+c6*sigma3_dot+(s6_dot*c5*c1*c2*s4+s6*c5_dot*c1*c2*s4+s6*c5*c1_dot*c2*s4+s6*c5*c1*c2_dot*s4+s6*c5*c1*c2*s4_dot)-(s6_dot*c5*c1*c3*c4*s2+s6*c5_dot*c1*c3*c4*s2+s6*c5*c1_dot*c3*c4*s2+s6*c5*c1*c3_dot*c4*s2+s6*c5*c1*c3*c4_dot*s2+s6*c5*c1*c3*c4*s2_dot)+(s6_dot*c1*s2*s3*s5+s6*c1_dot*s2*s3*s5+s6*c1*s2_dot*s3*s5+s6*c1*s2*s3_dot*s5+s6*c1*s2*s3*s5_dot))+d5*sigma3_dot+d3*(c1_dot*c2+c1*c2_dot);
    j13_dot = a8*((s6_dot*s5*sigma14+s6*s5_dot*sigma14+s6*s5*sigma14_dot)-(s6_dot*c4*c5*sigma11+s6*c4_dot*c5*sigma11+s6*c4*c5_dot*sigma11+s6*c4*c5*sigma11_dot)+(c6_dot*s4*sigma11+c6*s4_dot*sigma11+c6*s4*sigma11_dot))+d5*s4_dot*sigma11+d5*s4*sigma11_dot;
    j14_dot = a8*(c6_dot*sigma10+c6*sigma10_dot+(c5_dot*s6*sigma7+c5*s6_dot*sigma7+c5*s6*sigma7_dot))+d5*sigma10_dot;
    j15_dot = a8*((s6_dot*s5*sigma10+s6*s5_dot*sigma10+s6*s5*sigma10_dot)-(s6_dot*c5*sigma11+s6*c5_dot*sigma11+s6*c5*sigma11_dot));
    j16_dot = -a8*(s6_dot*sigma7+s6*sigma7_dot+c6_dot*sigma8+c6*sigma8_dot);
    j17_dot = 0;
    
    j21_dot = a8*(c6_dot*sigma7+c6*sigma7_dot-(s6_dot*sigma8+s6*sigma8_dot))+d5*sigma7_dot+d3*c1_dot*s2+d3*c1*s2_dot;
    j22_dot = d5*sigma2_dot+a8*(s6_dot*c5*c2*s1*s4+s6*c5_dot*c2*s1*s4+s6*c5*c2_dot*s1*s4+s6*c5*c2*s1_dot*s4+s6*c5*c2*s1*s4_dot)-a8*(s6_dot*c5*c3*c4*s1*s2+s6*c5_dot*c3*c4*s1*s2+s6*c5*c3_dot*c4*s1*s2+s6*c5*c3*c4_dot*s1*s2+s6*c5*c3*c4*s1_dot*s2+s6*c5*c3*c4*s1*s2_dot)+a8*(s6_dot*s1*s2*s3*s5+s6*s1_dot*s2*s3*s5+s6*s1*s2_dot*s3*s5+s6*s1*s2*s3_dot*s5+s6*s1*s2*s3*s5_dot)+a8*(c6_dot*sigma2+c6*sigma2_dot)+d3*(c2_dot*s1+c2*s1_dot);
    j23_dot = -a8*((s6_dot*s5*sigma15+s6*s5_dot*sigma15+s6*s5*sigma15_dot)-(s6_dot*c4*c5*sigma13+s6*c4_dot*c5*sigma13+s6*c4*c5_dot*sigma13+s6*c4*c5*sigma13_dot)+(c6_dot*s4*sigma13+c6*s4_dot*sigma13+c6*s4*sigma13_dot))-d5*(s4_dot*sigma13+s4*sigma13_dot);
    j24_dot = -a8*(c6_dot*sigma12+c6*sigma12_dot+c5_dot*s6*sigma6+c5*s6_dot*sigma6+c5*s6*sigma6_dot)-d5*sigma12_dot;
    j25_dot = -a8*((s6_dot*s5*sigma12+s6*s5_dot*sigma12+s6*s5*sigma12_dot)-(s6_dot*c5*sigma13+s6*c5_dot*sigma13+s6*c5*sigma13_dot));
    j26_dot = a8*(s6_dot*sigma6+s6*sigma6_dot+c6_dot*sigma9+c6*sigma9_dot);
    j27_dot = 0;
    
    j31_dot = 0;
    j32_dot = -a8*((s6_dot*c5*s2*s4+s6*c5_dot*s2*s4+s6*c5*s2_dot*s4+s6*c5*s2*s4_dot)+(s6_dot*c5*c2*c3*c4+s6*c5_dot*c2*c3*c4+s6*c5*c2_dot*c3*c4+s6*c5*c2*c3_dot*c4+s6*c5*c2*c3*c4_dot)-(s6_dot*c2*s3*s5+s6*c2_dot*s3*s5+s6*c2*s3_dot*s5+s6*c2*s3*s5_dot)+(s6_dot*c6*sigma4+s6*c6_dot*sigma4+s6*c6*sigma4_dot))-d5*sigma4_dot-d3*s2_dot;
    j33_dot = a8*((s6_dot*c3*s2*s5+s6*c3_dot*s2*s5+s6*c3*s2_dot*s5+s6*c3*s2*s5_dot)+(s6_dot*c4*c5*s2*s3+s6*c4_dot*c5*s2*s3+s6*c4*c5_dot*s2*s3+s6*c4*c5*s2_dot*s3+s6*c4*c5*s2*s3_dot)-(c6_dot*s2*s3*s4+c6*s2_dot*s3*s4+c6*s2*s3_dot*s4+c6*s2*s3*s4_dot)-(s2_dot*s3*s4+s2*s3_dot*s4+s2*s3*s4_dot));
    j34_dot = -d5*sigma1_dot-a8*((c6_dot*sigma1+c6*sigma1_dot)-(c5_dot*s6*sigma5+c5*s6_dot*sigma5+c5*s6*sigma5_dot));
    j35_dot = -a8*((s6_dot*s5*sigma1+s6*s5_dot*sigma1+s6*s5*sigma1_dot)-(s6_dot*c5*s2*s3+s6*c5_dot*s2*s3+s6*c5*s2_dot*s3+s6*c5*s2*s3_dot));
    j36_dot = a8*((c6_dot*c5*sigma1+c6*c5_dot*sigma1+c6*c5*sigma1_dot)+(c6_dot*s2*s3*s5+c6*s2_dot*s3*s5+c6*s2*s3_dot*s5+c6*s2*s3*s5_dot)-(s6_dot*sigma5+s6*sigma5_dot));
    j37_dot = 0;
    
    J_dot = [j11_dot j12_dot j13_dot j14_dot j15_dot j16_dot j17_dot;
             j21_dot j22_dot j23_dot j24_dot j25_dot j26_dot j27_dot;
             j31_dot j32_dot j33_dot j34_dot j35_dot j36_dot j37_dot];
end