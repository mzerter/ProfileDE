function d2fdx2val = chemo0_d2fdx2(times, y, logp, more)

[nt, ny] = size(y);

p = exp(logp);

p3  = p(3);
p4  = p(4);
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

d2fdx2val = zeros(nt,ny,ny,ny);

% y1-y1
d2fdx2val(:,1,1,1) =  (2*p12*y2)/( p10py1)^2 + ...
    (2*p12*y3)/( p11py1)^2 - (2*p12*y1*y2)/( p10py1)^3 - ...
    (2*p12*y1*y3)/( p11py1)^3;
d2fdx2val(:,2,1,1) = -(2*p10*p12*p9*y2)/p10py1^3;
d2fdx2val(:,3,1,1) = -(2*p11*p12*p9*y3)/p11py1^3;

% y2-y1
d2fdx2val(:,1,2,1) = -(p10*p12)/p10py1^2;
d2fdx2val(:,2,2,1) =  (p10*p12*p9)/p10py1^2;

% y3-y1
d2fdx2val(:,1,3,1) = -(p11*p12)/( p11py1)^2;
d2fdx2val(:,3,3,1) =  (p11*p12*p9)/( p11py1)^2;

% y1-y2
d2fdx2val(:,1,1,2) = -(p10*p12)/p10py1^2;
d2fdx2val(:,2,1,2) =  (p10*p12*p9)/p10py1^2;

% y2-y2
d2fdx2val(:,2,2,2) = (2*p13*p3^2*(p15 + p4*y3)*y4py5)/D1^3;
d2fdx2val(:,3,2,2) = -(2*p13*p3^2*p4*y3*y4py5)/D1^3;
d2fdx2val(:,4,2,2) = -(2*p13*p14*p15*p3^2*y4)/D1^3;

%  y3-y2
d2fdx2val(:,2,3,2) = (p13*p3*p4*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdx2val(:,3,3,2) = (p13*p3*p4*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdx2val(:,4,3,2) = -(2*p13*p14*p15*p3*p4*y4)/D1^3;

%  y4-y2
d2fdx2val(:,2,4,2) = -(p13*p3*(p15 + p4*y3))/D1^2;
d2fdx2val(:,3,4,2) = (p13*p3*p4*y3)/D1^2;
d2fdx2val(:,4,4,2) = (p13*p14*p15*p3)/D1^2;

% y5-y2
d2fdx2val(:,2,5,2) = -(p13*p3*(p15 + p4*y3))/D1^2;
d2fdx2val(:,3,5,2) = (p13*p3*p4*y3)/D1^2;

% y1-y3
d2fdx2val(:,1,1,3) = -(p11*p12)/p11py1^2;
d2fdx2val(:,3,1,3) =  (p11*p12*p9)/p11py1^2;

% y2-y3
d2fdx2val(:,2,2,3) = d2fdx2val(:,2,3,2);
d2fdx2val(:,3,2,3) = d2fdx2val(:,3,3,2);
d2fdx2val(:,4,2,3) = d2fdx2val(:,4,3,2);

% y3-y3
d2fdx2val(:,2,3,3) = -(2*p13*p3*p4^2*y2*y4py5)/D1^3;
d2fdx2val(:,3,3,3) = (2*p13*p4^2*(p15 + p3*y2)*y4py5)/D1^3;
d2fdx2val(:,4,3,3) = -(2*p13*p14*p15*p4^2*y4)/D1^3;

% y4-y3
d2fdx2val(:,2,4,3) = (p13*p3*p4*y2)/D1^2;
d2fdx2val(:,3,4,3) = -(p13*p4*(p15 + p3*y2))/D1^2;
d2fdx2val(:,4,4,3) = (p13*p14*p15*p4)/D1^2;

% y5-y3
d2fdx2val(:,2,5,3) = (p13*p3*p4*y2)/D1^2;
d2fdx2val(:,3,5,3) = -(p13*p4*(p15 + p3*y2))/D1^2;

% y2-y4
d2fdx2val(:,2,2,4) = -(p13*p3*(p15 + p4*y3))/D1^2;
d2fdx2val(:,3,2,4) = (p13*p3*p4*y3)/D1^2;
d2fdx2val(:,4,2,4) = (p13*p14*p15*p3)/D1^2;

% y3-y4
d2fdx2val(:,2,3,4) = (p13*p3*p4*y2)/D1^2;
d2fdx2val(:,3,3,4) = -(p13*p4*(p15 + p3*y2))/D1^2;
d2fdx2val(:,4,3,4) = (p13*p14*p15*p4)/D1^2;

% y2-y5
d2fdx2val(:,2,2,5) = d2fdx2val(:,2,5,2);
d2fdx2val(:,3,2,5) = d2fdx2val(:,3,5,2);

% y3-y5
d2fdx2val(:,2,3,5) = d2fdx2val(:,2,5,3);
d2fdx2val(:,3,3,5) = d2fdx2val(:,3,5,3);

end

