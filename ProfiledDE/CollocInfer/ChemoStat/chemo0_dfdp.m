function dfdpval = chemo0_dfdp(times, y, logp)

p = exp(logp);

[nt, ny] = size(y);  np = length(p);

p3  = p(3);
p4  = p(4);
p5  = p(5);
p6  = p(6);
p9  = p(9);
p10 = p(10);
p11 = p(11);
p12 = p(12);
p13 = p(13);
p14 = p(14);
p15 = p(15);

y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
y4 = y(:,4);
y5 = y(:,5);

p10py1 = p10 + y1;
p11py1 = p11 + y1;

y4py5  = y4 + y5;

D1 = p15 + p3*y2 + p4*y3;

dfdpval = zeros(nt,ny,np);

dfdpval(:,2,3) = -(p13*y2*(p15 + p4*y3)*y4py5)/D1^2;
dfdpval(:,3,3) =  (p13*p4*y2*y3*y4py5)/D1^2;
dfdpval(:,4,3) =  (p13*p14*p15*y2*y4)/D1^2;

dfdpval(:,2,4) =  (p13*p3*y2*y3*y4py5)/D1^2;
dfdpval(:,3,4) = -(p13*y3*(p15 + p3*y2)*y4py5)/D1^2;
dfdpval(:,4,4) =  (p13*p14*p15*y3*y4)/D1^2;

dfdpval(:,1,5) = p6;

dfdpval(:,1,6) = p5 - y1;
dfdpval(:,2,6) = -y2;
dfdpval(:,3,6) = -y3;
dfdpval(:,4,6) = -y4;
dfdpval(:,5,6) = -y5;

dfdpval(:,4,7) = -y4;
dfdpval(:,5,7) = -y5;

dfdpval(:,4,8) = -y4;
dfdpval(:,5,8) =  y4;

dfdpval(:,2,9) = (p12*y1*y2)/p10py1;
dfdpval(:,3,9) = (p12*y1*y3)/p11py1;

dfdpval(:,1,10) = (p12*y1*y2)/p10py1^2;
dfdpval(:,2,10) = -(p12*p9*y1*y2)/p10py1^2;

dfdpval(:,1,11) = (p12*y1*y3)/p11py1^2;
dfdpval(:,3,11) = -(p12*p9*y1*y3)/p11py1^2;

dfdpval(:,1,12) = -(y1*y2)/p10py1 - (y1*y3)/p11py1;
dfdpval(:,2,12) = (p9*y1*y2)/p10py1;
dfdpval(:,3,12) = (p9*y1*y3)/p11py1;

dfdpval(:,2,13) = -(p3*y2*y4py5)/D1;
dfdpval(:,3,13) = -(p4*y3*y4py5)/D1;
dfdpval(:,4,13) =  (p14*y4*(p3*y2 + p4*y3))/D1;

dfdpval(:,4,14) = (p13*y4*(p3*y2 + p4*y3))/D1;

dfdpval(:,2,15) =  (p13*p3*y2*y4py5)/D1^2;
dfdpval(:,3,15) =  (p13*p4*y3*y4py5)/D1^2;
dfdpval(:,4,15) = -(p13*p14*y4*(p3*y2 + p4*y3))/D1^2;

end

