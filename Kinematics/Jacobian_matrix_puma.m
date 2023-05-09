function J = Jacobian_matrix_puma(a,d,theta)
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
    
   j11 = d6*(-s1*(c23*(c4*c5*c6-s4*s6)-s23*s5*c6)+c1*(s4*c5*c6+c4*s6))+(a2*c2+a3*c23-d4*s23)*(-s1)-d3*c1;
  j12 = d6*(c1*((c2*(-s3)-s2*c3)*(c4*c5*c6-s4*s6)-(c2*c3-s2*s3)*s5*c6))+c1*(a2*(-s2)+a3*(-s2*c3-c2*s3)-d4*(-s2*s3+c2*c3));
  j13 = d6*(c1*((c2*(-s3)-s2*c3)*(c4*c5*c6-s4*s6)-(c2*c3-s2*s3)*s5*c6))+c1*(a3*(c2*(-s3)-s2*c3)-d4*(c2*c3-s2*s3));
  j14 = d6*(c1*(c23*(-s4*c5*c6-c4*s6))+s1*(c4*c5*c6-s4*s6));
  j15 = d6*(c1*(c23*(c4*(-s5)*c6)-s23*c5*c6)+s1*s4*(-s5)*c6);
  j16 = d6*(c1*(c23*(-c4*c5*s6-s4*c6)-s23*s5*(-s6))+s1*(s4*c5*(-s6)+c4*c6));

  j21 = d6*(c1*(c23*(c4*c5*c6-s4*s6)-s23*s5*c6)+s1*(s4*c5*c6+c4*s6))+c1*(a2*c2+a3*c23-d4*s23)-d3*s1;
  j22 = d6*(s1*((-s2*c3-c2*s3)*(c4*c5*c6-s4*s6)-(-s2*s3+c2*c3)*s5*c6))+s1*(a2*(-s2)+a3*(-s2*c3-c2*s3)-d4*(-s2*s3+c2*c3));
  j23 = d6*(s1*((-c2*s3-s2*c3)*(c4*c5*c6-s4*s6)-(-s2*s3+c2*c3)*s5*c6))+s1*(a3*(c2*(-s3)-s2*c3)-d4*(c2*c3-s2*s3));
  j24 = d6*(s1*(c23*(-s4*c5*c6-c4*s6))-c1*(c4*c5*c6-s4*s6));
  j25 = d6*(s1*(c23*(c4*(-s5)*c6)-s23*c5*c6)-c1*(s4*(-s5)*c6));
  j26 = d6*(s1*(c23*(c4*c5*(-s6)-s4*c6)+s23*s5*s6)-c1*(s4*c5*(-s6)+c4*c6));

  j31 = 0;
  j32 = d6*((s2*s3-c2*c3)*(c4*c5*c6-s4*s6)-(-s2*c3-c2*s3)*s5*c6)-a3*(-s2*s3+c2*c3)-a2*c2-d4*(-s2*c3-c2*s3);
  j33 = d6*(-(c2*c3-s2*s3)*(c4*c5*c6-s4*s6)-(c2*(-s3)-s2*c3)*s5*c6)-a3*(c2*c3-s2*s3)-d4*(-c2*s3-s2*c3);
  j34 = d6*(-s23*(-s4*c5*c6-c4*s6));
  j35 = d6*(-s23*(-c4*s5*c6)-c23*c5*c6);
  j36 = d6*(-s23*(-c4*c5*s6-s4*c6)+c23*s5*s6);
  % jacobian
  J = [j11,j12,j13,j14,j15,j16;
       j21,j22,j23,j24,j25,j26;
       j31,j32,j33,j34,j35,j36];
    
end