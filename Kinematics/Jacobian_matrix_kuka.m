function J = Jacobian_matrix_kuka(a8,d1,d3,d5,theta)
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

%defining elements of the jacobian matrix
    j11 = a8*(c6*sigma6-s6*sigma9)+d5*sigma6-d3*s1*s2;
    j12 = a8*(c6*sigma3+s6*(c5*(c1*c2*s4-c1*c3*c4*s2)+c1*s2*s3*s5))+d5*sigma3+d3*c1*c2;
    j13 = a8*(s6*(s5*sigma14-c4*c5*sigma11)+c6*s4*sigma11)+d5*s4*sigma11;
    j14 = a8*(c6*sigma10+c5*s6*sigma7)+d5*sigma10;
    j15 = a8*s6*(s5*sigma10-c5*sigma11);
    j16 = -a8*(s6*sigma7+c6*sigma8);
    j17=0;
    
    j21 = a8*(c6*sigma7-s6*sigma8)+d5*sigma7+d3*c1*s2;
    j22 = d5*sigma2+a8*(s6*(c5*(c2*s1*s4-c3*c4*s1*s2)+s1*s2*s3*s5)+c6*sigma2)+d3*c2*s1;
    j23 = -a8*(s6*(s5*sigma15-c4*c5*sigma13)+c6*s4*sigma13)-d5*s4*sigma13;
    j24 = -a8*(c6*sigma12+c5*s6*sigma6)-d5*sigma12;
    j25 = -a8*s6*(s5*sigma12-c5*sigma13);
    j26 = a8*(s6*sigma6+c6*sigma9);
    j27 =0;
    
    j31 = 0;
    j32 = -a8*(s6*(c5*(s2*s4+c2*c3*c4)-c2*s3*s5)+c6*sigma4)-d5*sigma4-d3*s2;
    j33= a8*(s6*(c3*s2*s5+c4*c5*s2*s3)-c6*s2*s3*s4)-d5*s2*s3*s4;
    j34 = -d5*sigma1-a8*(c6*sigma1-c5*s6*sigma5);
    j35 = -a8*s6*(s5*sigma1-c5*s2*s3);
    j36 = a8*(c6*(c5*sigma1+s2*s3*s5)-s6*sigma5);
    j37 = 0;
    
    J = [j11 j12 j13 j14 j15 j16 j17;
         j21 j22 j23 j24 j25 j26 j27;
         j31 j32 j33 j34 j35 j36 j37];
end