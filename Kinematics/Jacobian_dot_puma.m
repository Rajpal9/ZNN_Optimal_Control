function J_dot = Jacobian_dot_puma(a,d,theta,theta_dot)
    % robot params
    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    a4 = a(4);
    a5 = a(5);
    a6 = a(6);

    d1 = d(1);
    d2 = d(2);
    d3 = d(3);
    d4 = d(4);
    d5 = d(5);
    d6 = d(6);





    th1 = theta(1);
    th2 = theta(2);
    th3 = theta(3);
    th4 = theta(4);
    th5 = theta(5);
    th6 = theta(6);
    
    th1_dot = theta_dot(1);
   th2_dot = theta_dot(2);
   th3_dot = theta_dot(3);
   th4_dot = theta_dot(4);
   th5_dot = theta_dot(5);
   th6_dot = theta_dot(6);

    c1 = cos(th1);
    c2 = cos(th2);
    c3 = cos(th3);
    c4 = cos(th4);
    c5 = cos(th5);
    c6 = cos(th6);

    s1 = sin(th1);
    s2 = sin(th2);
    s3 = sin(th3);
    s4 = sin(th4);
    s5 = sin(th5);
    s6 = sin(th6);

    c23 = c2*c3-s2*s3;
    s23 = c2*s3+s2*c3;
    
   d_s1 = c1*th1_dot;
    d_s2 = c2*th2_dot;
    d_s3 = c3*th3_dot;
    d_s4 = c4*th4_dot;
    d_s5 = c5*th5_dot;
    d_s6 = c6*th6_dot;

    d_c1 = -s1*th1_dot;
    d_c2 = -s2*th2_dot;
    d_c3 = -s3*th3_dot;
    d_c4 = -s4*th4_dot;
    d_c5 = -s5*th5_dot;
    d_c6 = -s6*th6_dot;

    d_c23 = d_c2*c3+c2*d_c3-d_s2*s3-s2*d_s3;
    d_s23 = d_c2*s3+c2*d_s3+d_s2*c3+s2*d_c3;

    j11_dot = d6*(-(d_s1*c23*c4*c5*c6+s1*d_c23*c4*c5*c6+s1*c23*d_c4*c5*c6+s1*c23*c4*d_c5*c6+s1*c23*c4*c5*d_c6)+(d_s1*c23*s4*s6+s1*d_c23*s4*s6+s1*c23*d_s4*s6+s1*c23*s4*d_s6)+(d_s1*s23*s5*c6+s1*d_s23*s5*c6+s1*s23*d_s5*c6+s1*s23*s5*d_c6)+(d_c1*s4*c5*c6+c1*d_s4*c5*c6+c1*s4*d_c5*c6+c1*s4*c5*d_c6)+(d_c1*c4*s6+c1*d_c4*s6+c1*c4*d_s6))-(a2*d_s1*c2+a2*s1*d_c2)-(a3*d_s1*c23+a3*s1*d_c23)+(d4*d_s1*s23+d4*s1*d_s23)-(d3*d_c1);
    j12_dot = d6*(-(d_c1*s2*c3*c4*c5*c6+c1*d_s2*c3*c4*c5*c6+c1*s2*d_c3*c4*c5*c6+c1*s2*c3*d_c4*c5*c6+c1*s2*c3*c4*d_c5*c6+c1*s2*c3*c4*c5*d_c6)+(d_c1*s2*c3*s4*s6+c1*d_s2*c3*s4*s6+c1*s2*d_c3*s4*s6+c1*s2*c3*d_s4*s6+c1*s2*c3*s4*d_s6)-(d_c1*c2*s3*c4*c5*c6+c1*d_c2*s3*c4*c5*c6+c1*c2*d_s3*c4*c5*c6+c1*c2*s3*d_c4*c5*c6+c1*c2*s3*c4*d_c5*c6+c1*c2*s3*c4*c5*d_c6)+(d_c1*c2*s3*s4*s6+c1*d_c2*s3*s4*s6+c1*c2*d_s3*s4*s6+c1*c2*s3*d_s4*s6+c1*c2*s3*s4*d_s6)+(d_c1*s2*s3*s5*c6+c1*d_s2*s3*s5*c6+c1*s2*d_s3*s5*c6+c1*s2*s3*d_s5*c6+c1*s2*s3*s5*d_c6)-(d_c1*c2*c3*s5*c6+c1*d_c2*c3*s5*c6+c1*c2*d_c3*s5*c6+c1*c2*c3*d_s5*c6+c1*c2*c3*s5*d_c6))-(a2*d_c1*s2+a2*c1*d_s2)-(a3*d_c1*s2*s3+a3*c1*d_s2*s3+a3*c1*s2*d_s3)-(a3*d_c1*c2*s3+a3*c1*d_c2*s3+a3*c1*c2*d_s3)+(d4*d_c1*s2*s3+d4*c1*d_s2*s3+d4*c1*s2*d_s3)-(d4*d_c1*c2*c3+d4*c1*d_c2*c3+d4*c1*c2*d_c3);
    j13_dot = d6*(-(d_c1*c2*s3*c4*c5*c6+c1*d_c2*s3*c4*c5*c6+c1*c2*d_s3*c4*c5*c6+c1*c2*s3*d_c4*c5*c6+c1*c2*s3*c4*d_c5*c6+c1*c2*s3*c4*c5*d_c6)+(d_c1*c2*s3*s4*s6+c1*d_c2*s3*s4*s6+c1*c2*d_s3*s4*s6+c1*c2*s3*d_s4*s6+c1*c2*s3*s4*d_s6)-(d_c1*s2*c3*c4*c5*c6+c1*d_s2*c3*c4*c5*c6+c1*s2*d_c3*c4*c5*c6+c1*s2*c3*d_c4*c5*c6+c1*s2*c3*c4*d_c5*c6+c1*s2*c3*c4*c5*d_c6)+(d_c1*s2*c3*s4*s6+c1*d_s2*c3*s4*s6+c1*s2*d_c3*s4*s6+c1*s2*c3*d_s4*s6+c1*s2*c3*s4*d_s6)-(d_c1*c2*c3*s5*c6+c1*d_c2*c3*s5*c6+c1*c2*d_c3*s5*c6+c1*c2*c3*d_s5*c6+c1*c2*c3*s5*d_c6)+(d_c1*s2*s3*s5*c6+c1*d_s2*s3*s5*c6+c1*s2*d_s3*s5*c6+c1*s2*s3*d_s5*c6+c1*s2*s3*s5*d_c6))-(a3*d_c1*c2*s3+a3*c1*d_c2*s3+a3*c1*c2*d_s3)-(a3*d_c1*s2*c3+a3*c1*d_s2*c3+a3*c1*s2*d_c3)-(d4*d_c1*c2*c3+d4*c1*d_c2*c3+d4*c1*c2*d_c3)+(d4*d_c1*s2*s3+d4*c1*d_s2*s3+d4*c1*s2*d_s3);
    j14_dot = d6*(-(d_c1*c23*s4*c5*c6+c1*d_c23*s4*c5*c6+c1*c23*d_s4*c5*c6+c1*c23*s4*d_c5*c6+c1*c23*s4*c5*d_c6)-(d_c1*c23*c4*s6+c1*d_c23*c4*s6+c1*c23*d_c4*s6+c1*c23*c4*d_s6)+(d_s1*c4*c5*c6+s1*d_c4*c5*c6+s1*c4*c5*c6+s1*c4*d_c5*c6+s1*c4*c5*d_c6)-(d_s1*s4*s6+s1*d_s4*s6+s1*s4*d_s6));
    j15_dot = d6*(-(d_c1*c23*c4*s5*c6+c1*d_c23*c4*s5*c6+c1*c23*d_c4*s5*c6+c1*c23*c4*d_s5*c6+c1*c23*c4*s5*d_c6)+(d_c1*s23*c5*c6+c1*d_s23*c5*c6+c1*s23*d_c5*c6+c1*s23*c5*d_c6)-(d_s1*s4*s5*c6+s1*d_s4*s5*c6+s1*s4*d_s5*c6+s1*s4*s5*d_c6));
    j16_dot = d6*(-(d_c1*c23*c4*c5*c6+c1*d_c23*c4*c5*c6+c1*c23*d_c4*c5*c6+c1*c23*c4*d_c5*c6+c1*c23*c4*c5*d_c6)-(d_c1*c23*s4*c6+c1*d_c23*s4*c6+c1*c23*d_s4*c6+c1*c23*s4*d_c6)+(d_c1*s23*s5*s6+c1*d_s23*s5*s6+c1*s23*d_s5*s6+c1*s23*s5*d_s6)-(d_s1*s4*c5*s6+s1*d_s4*c5*s6+s1*s4*d_c5*s6+s1*s4*c5*d_s6)+(d_s1*c4*c6+s1*d_c4*c6+s1*c4*d_c6));

    j21_dot = d6*((d_c1*c23*c4*c5*c6+c1*d_c23*c4*c5*c6+c1*c23*d_c4*c5*c6+c1*c23*c4*d_c5*c6+c1*c23*c4*c5*d_c6)-(d_c1*c23*s4*s6+c1*d_c23*s4*s6+c1*c23*d_s4*s6+c1*c23*s4*d_s6)-(d_c1*s23*s5*c6+c1*d_s23*s5*c6+c1*s23*d_s5*c6+c1*s23*s5*d_c6)+(d_s1*s4*c5*c6+s1*d_s4*c5*c6+s1*s4*d_c5*c6+s1*s4*c5*d_c6)+(d_s1*c4*s6+s1*d_c4*s6+s1*c4*d_s6))+(a2*d_c1*c2+a2*c1*d_c2)+(a3*d_c1*c23+a3*c1*d_c23)-(d4*d_c1*s23+d4*c1*d_s23)-(d3*d_s1);
    j22_dot = d6*(-(d_s1*s2*c3*c4*c5*c6+s1*d_s2*c3*c4*c5*c6+s1*s2*d_c3*c4*c5*c6+s1*s2*c3*d_c4*c5*c6+s1*s2*c3*c4*d_c5*c6+s1*s2*c3*c4*c5*d_c6)+(d_s1*s2*c3*s4*s6+s1*d_s2*c3*s4*s6+s1*s2*d_c3*s4*s6+s1*s2*c3*d_s4*s6+s1*s2*c3*s4*d_s6)-(d_s1*c2*s3*c4*c5*c6+s1*d_c2*s3*c4*c5*c6+s1*c2*d_s3*c4*c5*c6+s1*c2*s3*d_c4*c5*c6+s1*c2*s3*c4*d_c5*c6+s1*c2*s3*c4*c5*d_c6)+(d_s1*c2*s3*s4*s6+s1*d_c2*s3*s4*s6+s1*c2*d_s3*s4*s6+s1*c2*s3*d_s4*s6+s1*c2*s3*s4*d_s6)+(d_s1*s2*s3*s5*c6+s1*d_s2*s3*s5*c6+s1*s2*d_s3*s5*c6+s1*s2*s3*d_s5*c6+s1*s2*s3*s5*d_c6)-(d_s1*c2*c3*s5*c6+s1*d_c2*c3*s5*c6+s1*c2*d_c3*s5*c6+s1*c2*c3*d_s5*c6+s1*c2*c3*s5*d_c6))-(a2*d_s1*s2+a2*s1*d_s2)-(a3*d_s1*s2*c3+a3*s1*d_s2*c3+a3*s1*s2*d_c3)-(a3*d_s1*c2*s3+a3*s1*d_c2*s3+a3*s1*c2*d_s3)+(d4*d_s1*s2*s3+d4*s1*d_s2*s3+d4*s1*s2*d_s3)-(d4*d_s1*c2*c3+d4*s1*d_c2*c3+d4*s1*c2*d_c3);
    j23_dot = d6*(-(d_s1*c2*s3*c4*c5*c6+s1*d_c2*s3*c4*c5*c6+s1*c2*d_s3*c4*c5*c6+s1*c2*s3*d_c4*c5*c6+s1*c2*s3*c4*d_c5*c6+s1*c2*s3*c4*c5*d_c6)+(d_s1*c2*s3*s4*s6+s1*d_c2*s3*s4*s6+s1*c2*d_s3*s4*s6+s1*c2*s3*d_s4*s6+s1*c2*s3*s4*d_s6)-(d_s1*s2*c3*c4*c5*c6+s1*d_s2*c3*c4*c5*c6+s1*s2*d_c3*c4*c5*c6+s1*s2*c3*d_c4*c5*c6+s1*s2*c3*c4*d_c5*c6+s1*s2*c3*c4*c5*d_c6)+(d_s1*s2*c3*s4*s6+s1*d_s2*c3*s4*s6+s1*s2*d_c3*s4*s6+s1*s2*c3*d_s4*s6+s1*s2*c3*s4*d_s6)+(d_s1*s2*s3*s5*c6+s1*d_s2*s3*s5*c6+s1*s2*d_s3*s5*c6+s1*s2*s3*d_s5*c6+s1*s2*s3*s5*d_c6)-(d_s1*c2*c3*s5*c6+s1*d_c2*c3*s5*c6+s1*c2*d_c3*s5*c6+s1*c2*c3*d_s5*c6+s1*c2*c3*s5*d_c6))-(a3*d_s1*c2*s3+a3*s1*d_c2*s3+a3*s1*c2*d_s3)-(a3*d_s1*s2*c3+a3*s1*d_s2*c3+a3*s1*s2*d_c3)-(d4*d_s1*c2*c3+d4*s1*d_c2*c3+d4*s1*c2*d_c3)+(d4*d_s1*s2*s3+d4*s1*d_s2*s3+d4*s1*s2*d_s3);
    j24_dot = d6*(-(d_s1*c23*s4*c5*c6+s1*d_c23*s4*c5*c6+s1*c23*d_s4*c5*c6+s1*c23*s4*d_c5*c6+s1*c23*s4*c5*d_c6)-(d_s1*c4*s6*c23+s1*d_c4*s6*c23+s1*c4*d_s6*c23+s1*c4*s6*d_c23)-(d_c1*c4*c5*c6+c1*d_c4*c5*c6+c1*c4*d_c5*c6+c1*c4*c5*d_c6)+(d_c1*s4*s6+c1*d_s4*s6+c1*s4*d_s6));
    j25_dot = d6*(-(d_s1*c23*c4*s5*c6+s1*d_c23*c4*s5*c6+s1*c23*d_c4*s5*c6+s1*c23*c4*d_s5*c6+s1*c23*c4*s5*d_c6)-(d_s1*s23*c5*c6+s1*d_s23*c5*c6+s1*s23*d_c5*c6+s1*s23*c5*d_c6)+(d_c1*s4*s5*c6+c1*d_s4*s5*c6+c1*s4*d_s5*c6+c1*s4*s5*d_c6));
    j26_dot = d6*(-(d_s1*c23*c4*c5*s6+s1*d_c23*c4*c5*s6+s1*c23*d_c4*c5*s6+s1*c23*c4*d_c5*s6+s1*c23*c4*c5*d_s6)-(d_s1*c23*s4*c6+s1*d_c23*s4*c6+s1*c23*d_s4*c6+s1*c23*s4*d_c6)+(d_s1*s23*s5*s6+s1*d_s23*s5*s6+s1*s23*d_s5*s6+s1*s23*s5*d_s6)+(d_c1*s4*c5*s6+c1*d_s4*c5*s6+c1*s4*d_c5*s6+c1*s4*c5*d_s6)-(d_c1*c4*c6+c1*d_c4*c6+c1*c4*d_c6));

    j31_dot = 0;
    j32_dot = d6*((s2*s3*c4*c5*c6+s2*s3*c4*c5*c6+s2*s3*c4*c5*c6+s2*s3*c4*c5*c6+s2*s3*c4*c5*c6)-(d_s2*s3*s4*s6+s2*d_s3*s4*s6+s2*s3*d_s4*s6+s2*s3*s4*d_s6)-(d_c2*c3*c4*c5*c6+c2*d_c3*c4*c5*c6+c2*c3*d_c4*c5*c6+c2*c3*c4*d_c5*c6+c2*c3*c4*c5*d_c6)+(d_c2*c3*s4*s6+c2*d_c3*s4*s6+c2*c3*d_s4*s6+c2*c3*s4*d_s6)+(d_s2*c3*s5*c6+s2*d_c3*s5*c6+s2*c3*d_s5*c6+s2*c3*s5*d_c6)+(d_c2*s3*s5*c6+c2*d_s3*s5*c6+c2*s3*d_s5*c6+c2*s3*s5*d_c6))+a3*(d_s2*s3+s2*d_s3)-a3*(d_c2*c3+c2*d_c3)-a2*d_c2-d4*(-(d_s2*c3+s2*d_c3)-(d_c2*s3+c2*d_s3));
    j33_dot = d6*(-(d_c2*c3*c4*c5*c6+c2*d_c3*c4*c5*c6+c2*c3*d_c4*c5*c6+c2*c3*c4*d_c5*c6+c2*c3*c4*c5*d_c6)+(d_c2*c3*s4*s6+c2*d_c3*s4*s6+c2*c3*d_s4*s6+c2*c3*s4*d_s6)+(d_s2*s3*c4*c5*c6+s2*d_s3*c4*c5*c6+s2*s3*d_c4*c5*c6+s2*s3*c4*d_c5*c6+s2*s3*c4*c5*d_c6)-(d_s2*s3*s4*s6+s2*d_s3*s4*s6+s2*s3*d_s4*s6+s2*s3*s4*d_s6)+(d_c2*s3*s5*c6+c2*d_s3*s5*c6+c2*s3*d_s5*c6+c2*s3*s5*d_c6)+(d_s2*c3*s5*c6+s2*d_c3*s5*c6+s2*c3*d_s5*c6+s2*c3*s5*d_c6))-a3*(d_c2*c3+c2*d_c3)+a3*(d_s2*s3+s2*d_s3)+d4*(d_s2*c3+s2*d_c3)+d4*(d_c2*s3+c2*d_s3);
    j34_dot = d6*((d_s23*s4*c5*c6+s23*d_s4*c5*c6+s23*s4*d_c5*c6+s23*s4*c5*d_c6)+(d_s23*c4*s6+s23*d_c4*s6+s23*c4*d_s6));
    j35_dot = d6*((d_s23*c4*s5*c6+s23*d_c4*s5*c6+s23*c4*d_s5*c6+s23*c4*s5*d_c6)-(d_c23*c5*c6+c23*d_c5*c6+c23*c5*d_c6));
    j36_dot = d6*((d_s3*c4*c5*s6+s3*d_c4*c5*s6+s3*c4*d_c5*s6+s3*c4*c5*d_s6)+(d_s23*s4*c6+s23*d_s4*c6+s23*s4*d_c6)+(d_c23*s5*s6+c23*d_s5*s6+c23*s5*d_s6));
% Derivative of jacobian
    J_dot = [j11_dot j12_dot j13_dot j14_dot j15_dot j16_dot;
             j21_dot j22_dot j23_dot j24_dot j25_dot j26_dot;
             j31_dot j32_dot j33_dot j34_dot j35_dot j36_dot];
end