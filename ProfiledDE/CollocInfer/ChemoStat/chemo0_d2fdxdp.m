function d2fdxdpval = chemo0_d2fdxdp(times, y, logp, more)

[nt, ny] = size(y);
np = length(logp);

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
p10my1 = p10 - y1;
p11my1 = p11 - y1;

y4py5  = y4 + y5;

D1 = p15 + p3*y2 + p4*y3;

d2fdxdpval = zeros(nt,ny,ny,np);

% y2-p3
d2fdxdpval(:,2,2,3) = -(p13*(p15 + p4*y3)*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdxdpval(:,3,2,3) = (p13*p4*y3*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdxdpval(:,4,2,3) = (p13*p14*p15*y4*(p15 - p3*y2 + p4*y3))/D1^3;

% y3-p3
d2fdxdpval(:,2,3,3) = (p13*p4*y2*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdxdpval(:,3,3,3) = (p13*p4*y2*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdxdpval(:,4,3,3) = -(2*p13*p14*p15*p4*y2*y4)/D1^3;

% y4-p3
d2fdxdpval(:,2,4,3) = -(p13*y2*(p15 + p4*y3))/D1^2;
d2fdxdpval(:,3,4,3) = (p13*p4*y2*y3)/D1^2;
d2fdxdpval(:,4,4,3) = (p13*p14*p15*y2)/D1^2;

% y5-p3
d2fdxdpval(:,2,5,3) = -(p13*y2*(p15 + p4*y3))/D1^2;
d2fdxdpval(:,3,5,3) = (p13*p4*y2*y3)/D1^2;

% y2-p4
d2fdxdpval(:,2,2,4) = (p13*p3*y3*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdxdpval(:,3,2,4) = (p13*p3*y3*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdxdpval(:,4,2,4) = -(2*p13*p14*p15*p3*y3*y4)/D1^3;

% y3-p4
d2fdxdpval(:,2,3,4) = (p13*p3*y2*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdxdpval(:,3,3,4) = -(p13*(p15 + p3*y2)*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdxdpval(:,4,3,4) = (p13*p14*p15*y4*(p15 + p3*y2 - p4*y3))/D1^3;

% y4-p4
d2fdxdpval(:,2,4,4) = (p13*p3*y2*y3)/D1^2;
d2fdxdpval(:,3,4,4) = -(p13*y3*(p15 + p3*y2))/D1^2;
d2fdxdpval(:,4,4,4) = (p13*p14*p15*y3)/D1^2;

% y5-p4
d2fdxdpval(:,2,5,4) = (p13*p3*y2*y3)/D1^2;
d2fdxdpval(:,3,5,4) = -(p13*y3*(p15 + p3*y2))/D1^2;

% y1-p6
d2fdxdpval(:,1,1,6) = -1;

% y2-p6
d2fdxdpval(:,2,2,6) = -1;

% y3-p6
d2fdxdpval(:,3,3,6) = -1;

% y4-p6
d2fdxdpval(:,4,4,6) = -1;

% y5-p6
d2fdxdpval(:,5,5,6) = -1;

% y4-p7
d2fdxdpval(:,4,4,7) = -1;

% y5-p7
d2fdxdpval(:,5,5,7) = -1;

% y4-p8
d2fdxdpval(:,4,4,8) = -1;
d2fdxdpval(:,5,4,8) = 1;

% y1-p9
d2fdxdpval(:,2,1,9) = (p10*p12*y2)/p10py1^2;
d2fdxdpval(:,3,1,9) = (p11*p12*y3)/p11py1^2;

% y2-p9
d2fdxdpval(:,2,2,9) = (p12*y1)/p10py1;

% y3-p9
d2fdxdpval(:,3,3,9) = (p12*y1)/p11py1;

% y1-p10
d2fdxdpval(:,1,1,10) = (p12*y2*p10my1)/p10py1^3;
d2fdxdpval(:,2,1,10) = -(p12*p9*y2*p10my1)/p10py1^3;

% y2-p10
d2fdxdpval(:,1,2,10) = (p12*y1)/p10py1^2;
d2fdxdpval(:,2,2,10) = -(p12*p9*y1)/p10py1^2;

% y1-p11
d2fdxdpval(:,1,1,11) = (p12*y3*p11my1)/p11py1^3;
d2fdxdpval(:,3,1,11) = -(p12*p9*y3*p11my1)/p11py1^3;

% y3-p11
d2fdxdpval(:,1,3,11) = (p12*y1)/p11py1^2;
d2fdxdpval(:,3,3,11) = -(p12*p9*y1)/p11py1^2;

% y1-p12
d2fdxdpval(:,1,1,12) = -(p11*y3)/p11py1^2 - ...
           (p10*(y2*p11^2 + 2*y2*p11*y1 + y2*y1^2))/(p10py1^2*p11py1^2);
d2fdxdpval(:,2,1,12) = (p10*p9*y2)/p10py1^2;
d2fdxdpval(:,3,1,12) = (p11*p9*y3)/p11py1^2;

% y2-p12
d2fdxdpval(:,1,2,12) = -y1/p10py1;
d2fdxdpval(:,2,2,12) = (p9*y1)/p10py1;

% y3-p12
d2fdxdpval(:,1,3,12) = -y1/p11py1;
d2fdxdpval(:,3,3,12) = (p9*y1)/p11py1;

% y2-p13
d2fdxdpval(:,2,2,13) = -(p3*(p15 + p4*y3)*y4py5)/D1^2;
d2fdxdpval(:,3,2,13) = (p3*p4*y3*y4py5)/D1^2;
d2fdxdpval(:,4,2,13) = (p14*p15*p3*y4)/D1^2;

% y3-p13
d2fdxdpval(:,2,3,13) = (p3*p4*y2*y4py5)/D1^2;
d2fdxdpval(:,3,3,13) = -(p4*(p15 + p3*y2)*y4py5)/D1^2;
d2fdxdpval(:,4,3,13) = (p14*p15*p4*y4)/D1^2;

% y4-p13
d2fdxdpval(:,2,4,13) = -(p3*y2)/D1;
d2fdxdpval(:,3,4,13) = -(p4*y3)/D1;
d2fdxdpval(:,4,4,13) = p14 - (p14*p15)/D1;

% y5-p13
d2fdxdpval(:,2,5,13) = -(p3*y2)/D1;
d2fdxdpval(:,3,5,13) = -(p4*y3)/D1;

% y2-p14
d2fdxdpval(:,4,2,14) = (p13*p15*p3*y4)/D1^2;

% y3-p14
d2fdxdpval(:,4,3,14) = (p13*p15*p4*y4)/D1^2;

% y4-p14
d2fdxdpval(:,4,4,14) = p13 - (p13*p15)/D1;

% y2-p15
d2fdxdpval(:,2,2,15) = (p13*p3*y4py5*(p15 - p3*y2 + p4*y3))/D1^3;
d2fdxdpval(:,3,2,15) = -(2*p13*p3*p4*y3*y4py5)/D1^3;
d2fdxdpval(:,4,2,15) = (p13*p14*p3*y4*(p3*y2 - p15 + p4*y3))/D1^3;

% y3-p15
d2fdxdpval(:,2,3,15) = -(2*p13*p3*p4*y2*y4py5)/D1^3;
d2fdxdpval(:,3,3,15) = (p13*p4*y4py5*(p15 + p3*y2 - p4*y3))/D1^3;
d2fdxdpval(:,4,3,15) = (p13*p14*p4*y4*(p3*y2 - p15 + p4*y3))/D1^3;

% y4-p15
d2fdxdpval(:,2,4,15) = (p13*p3*y2)/D1^2;
d2fdxdpval(:,3,4,15) = (p13*p4*y3)/D1^2;
d2fdxdpval(:,4,4,15) = -(p13*p14*(p3*y2 + p4*y3))/D1^2;

% y5-p15
d2fdxdpval(:,2,5,15) = (p13*p3*y2)/D1^2;
d2fdxdpval(:,3,5,15) = (p13*p4*y3)/D1^2;

end

