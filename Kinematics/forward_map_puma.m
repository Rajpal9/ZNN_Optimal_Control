function X = forward_map_puma(a,d,theta)

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
    
   X =  [d6*(c1*(c23*(c4*c5*c6-s4*s6)-s23*s5*c6)+s1*(s4*c5*c6+c4*s6))+c1*(a2*c2+a3*c23-d4*s23)-d3*s1;
         d6*(s1*(c23*(c4*c5*c6-s4*s6)-s23*s5*c6)-c1*(s4*c5*c6+c4*s6))+s1*(a2*c2+a3*c23-d4*s23)+d3*c1;
         d6*(-s23*(c4*c5*c6-s4*s6)-c23*s5*c6)-a3*s23-a2*s2-d4*c23];
    
    
end