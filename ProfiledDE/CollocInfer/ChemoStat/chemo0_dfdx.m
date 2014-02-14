function dfdxval = chemo0_dfdx(times, y, logp)

p = exp(logp);

[nt, ny] = size(y);

p3  = p(3);
p4  = p(4);
p5  = p(5);
p6  = p(6);
p7  = p(7);
p8  = p(8);
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

dfdxval = zeros(nt,ny,ny);

dfdxval(:,1,1) = (p12*y1*y2)/p10py1^2 - (p12*y2)/p10py1 - ...
                 (p12*y3)/p11py1 - p6 + (p12*y1*y3)/p11py1^2;
dfdxval(:,1,2) = (p10.*p12.*p9.*y2)./p10py1.^2;
dfdxval(:,1,3) = (p11.*p12.*p9.*y3)./p11py1.^2;

dfdxval(:,2,1) = -(p12.*y1)./p10py1;
dfdxval(:,2,2) = (p12*p9*y1)/p10py1 - (p13*p3*y4py5)/D1 - p6 + ...
                 (p13*p3^2*y2*y4py5)/D1^2;
dfdxval(:,2,3) = (p13*p3*p4*y3*y4py5)/D1^2;
dfdxval(:,2,4) = (p13*p14*p15*p3*y4)/D1^2;

dfdxval(:,3,1) = -(p12.*y1)./p11py1;
dfdxval(:,3,2) = (p13*p3*p4*y2*y4py5)/D1^2;
dfdxval(:,3,3) = (p12*p9*y1)/p11py1 - (p13*p4*y4py5)/D1 - p6 + ...
dfdxval(:,3,4) = (p13*p14*p15*p4*y4)/D1^2; 

dfdxval(:,4,2) = -(p13*p3*y2)/D1;
dfdxval(:,4,3) = (p13*p14*p15*p4*y4)/D1^2; 
dfdxval(:,4,4) =  (p13*p14*(p3*y2 + p4*y3))/D1 - p7 - p8 - p6;
dfdxval(:,4,5) =  p8;

dfdxval(:,5,2) = -(p13*p3*y2)/D1;
dfdxval(:,5,5) = -p6 - p7;

% for i=1:nt
%     for j=1:ny
%         dfdxval(i,j,:) = dfdxval(i,j,:)/y(i,j);
%         dfdxval(i,j,j) = dfdxval(i,j,j) - fnval(i,j)/y(i,j)^2;
%     end
% end

end

