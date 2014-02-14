Q = p3*y2 + p4*y3;
Qs = Q*exp(100*(Q - p16))/(1 + exp(100*(Q - p16))) + ...
     p16/(1 + exp(100*(Q - p16)));

%  define derivatives of Q and Qs wrt five y's and 2 Q's

simplify(diff(Q,'y1'))
simplify(diff(Qs,'y1'))

simplify(diff(Q,'y2'))
simplify(diff(Qs,'y2'))

simplify(diff(Q,'y3'))
simplify(diff(Qs,'y3'))

simplify(diff(Q,'y4'))
simplify(diff(Qs,'y4'))

simplify(diff(Q,'y5'))
simplify(diff(Qs,'y5'))

%  define second derivatives of Q and Qs wrt five y's and 2 Q's

%  y1

simplify(diff(diff(Q,'y1'),'y1'))
simplify(diff(diff(Qs,'y1'),'y1'))

simplify(diff(diff(Q,'y2'),'y1'))
simplify(diff(diff(Qs,'y2'),'y1'))

simplify(diff(diff(Q,'y3'),'y1'))
simplify(diff(diff(Qs,'y3'),'y1'))

simplify(diff(diff(Q,'y4'),'y1'))
simplify(diff(diff(Qs,'y4'),'y1'))

simplify(diff(diff(Q,'y5'),'y1'))
simplify(diff(diff(Qs,'y5'),'y1'))

%  y2

simplify(diff(diff(Q,'y1'),'y2'))
simplify(diff(diff(Qs,'y1'),'y2'))

simplify(diff(diff(Q,'y2'),'y2'))
simplify(diff(diff(Qs,'y2'),'y2'))

simplify(diff(diff(Q,'y3'),'y2'))
simplify(diff(diff(Qs,'y3'),'y2'))

simplify(diff(diff(Q,'y4'),'y2'))
simplify(diff(diff(Qs,'y4'),'y2'))

simplify(diff(diff(Q,'y5'),'y2'))
simplify(diff(diff(Qs,'y5'),'y2'))

%  y3

simplify(diff(diff(Q,'y1'),'y3'))
simplify(diff(diff(Qs,'y1'),'y3'))

simplify(diff(diff(Q,'y2'),'y3'))
simplify(diff(diff(Qs,'y2'),'y3'))

simplify(diff(diff(Q,'y3'),'y3'))
simplify(diff(diff(Qs,'y3'),'y3'))

simplify(diff(diff(Q,'y4'),'y3'))
simplify(diff(diff(Qs,'y4'),'y3'))

simplify(diff(diff(Q,'y5'),'y3'))
simplify(diff(diff(Qs,'y5'),'y3'))

%  y4

simplify(diff(diff(Q,'y1'),'y4'))
simplify(diff(diff(Qs,'y1'),'y4'))

simplify(diff(diff(Q,'y2'),'y4'))
simplify(diff(diff(Qs,'y2'),'y4'))

simplify(diff(diff(Q,'y3'),'y4'))
simplify(diff(diff(Qs,'y3'),'y4'))

simplify(diff(diff(Q,'y4'),'y4'))
simplify(diff(diff(Qs,'y4'),'y4'))

simplify(diff(diff(Q,'y5'),'y4'))
simplify(diff(diff(Qs,'y5'),'y4'))

%  y5

simplify(diff(diff(Q,'y1'),'y5'))
simplify(diff(diff(Qs,'y1'),'y5'))

simplify(diff(diff(Q,'y2'),'y5'))
simplify(diff(diff(Qs,'y2'),'y5'))

simplify(diff(diff(Q,'y3'),'y5'))
simplify(diff(diff(Qs,'y3'),'y5'))

simplify(diff(diff(Q,'y4'),'y5'))
simplify(diff(diff(Qs,'y4'),'y5'))

simplify(diff(diff(Q,'y5'),'y5'))
simplify(diff(diff(Qs,'y5'),'y5'))

%  derivatives of 5 f's wrt last 14 p's

simplify(diff(Q,'p3'))
simplify(diff(Qs,'p3'))

simplify(diff(Q,'p4'))
simplify(diff(Qs,'p4'))

simplify(diff(Q,'p5'))
simplify(diff(Qs,'p5'))

simplify(diff(Q,'p6'))
simplify(diff(Qs,'p6'))

simplify(diff(Q,'p7'))
simplify(diff(Qs,'p7'))

simplify(diff(Q,'p8'))
simplify(diff(Qs,'p8'))

simplify(diff(Q,'p9'))
simplify(diff(Qs,'p9'))

simplify(diff(Q,'p10'))
simplify(diff(Qs,'p10'))

simplify(diff(Q,'p11'))
simplify(diff(Qs,'p11'))

simplify(diff(Q,'p12'))
simplify(diff(Qs,'p12'))

simplify(diff(Q,'p13'))
simplify(diff(Qs,'p13'))

simplify(diff(Q,'p14'))
simplify(diff(Qs,'p14'))

simplify(diff(Q,'p15'))
simplify(diff(Qs,'p15'))

simplify(diff(Q,'p16'))
simplify(diff(Qs,'p16'))

%     second derivatives of f's wrt y's and Q's and wrt last 14 p's

%  p3

simplify(diff(diff(Q,'y1'),'p3'))
simplify(diff(diff(Qs,'y1'),'p3'))

simplify(diff(diff(Q,'y2'),'p3'))
simplify(diff(diff(Qs,'y2'),'p3'))

simplify(diff(diff(Q,'y3'),'p3'))
simplify(diff(diff(Qs,'y3'),'p3'))

simplify(diff(diff(Q,'y4'),'p3'))
simplify(diff(diff(Qs,'y4'),'p3'))

simplify(diff(diff(Q,'y5'),'p3'))
simplify(diff(diff(Qs,'y5'),'p3'))

%  p4

simplify(diff(diff(Q,'y1'),'p4'))
simplify(diff(diff(Qs,'y1'),'p4'))

simplify(diff(diff(Q,'y2'),'p4'))
simplify(diff(diff(Qs,'y2'),'p4'))

simplify(diff(diff(Q,'y3'),'p4'))
simplify(diff(diff(Qs,'y3'),'p4'))

simplify(diff(diff(Q,'y4'),'p4'))
simplify(diff(diff(Qs,'y4'),'p4'))

simplify(diff(diff(Q,'y5'),'p4'))
simplify(diff(diff(Qs,'y5'),'p4'))

%  p5

simplify(diff(diff(Q,'y1'),'p5'))
simplify(diff(diff(Qs,'y1'),'p5'))

simplify(diff(diff(Q,'y2'),'p5'))
simplify(diff(diff(Qs,'y2'),'p5'))

simplify(diff(diff(Q,'y3'),'p5'))
simplify(diff(diff(Qs,'y3'),'p5'))

simplify(diff(diff(Q,'y4'),'p5'))
simplify(diff(diff(Qs,'y4'),'p5'))

simplify(diff(diff(Q,'y5'),'p5'))
simplify(diff(diff(Qs,'y5'),'p5'))

%  p6

simplify(diff(diff(Q,'y1'),'p6'))
simplify(diff(diff(Qs,'y1'),'p6'))

simplify(diff(diff(Q,'y2'),'p6'))
simplify(diff(diff(Qs,'y2'),'p6'))

simplify(diff(diff(Q,'y3'),'p6'))
simplify(diff(diff(Qs,'y3'),'p6'))

simplify(diff(diff(Q,'y4'),'p6'))
simplify(diff(diff(Qs,'y4'),'p6'))

simplify(diff(diff(Q,'y5'),'p6'))
simplify(diff(diff(Qs,'y5'),'p6'))

%  p7

simplify(diff(diff(Q,'y1'),'p7'))
simplify(diff(diff(Qs,'y1'),'p7'))

simplify(diff(diff(Q,'y2'),'p7'))
simplify(diff(diff(Qs,'y2'),'p7'))

simplify(diff(diff(Q,'y3'),'p7'))
simplify(diff(diff(Qs,'y3'),'p7'))

simplify(diff(diff(Q,'y4'),'p7'))
simplify(diff(diff(Qs,'y4'),'p7'))

simplify(diff(diff(Q,'y5'),'p7'))
simplify(diff(diff(Qs,'y5'),'p7'))

%  p8

simplify(diff(diff(Q,'y1'),'p8'))
simplify(diff(diff(Qs,'y1'),'p8'))

simplify(diff(diff(Q,'y2'),'p8'))
simplify(diff(diff(Qs,'y2'),'p8'))

simplify(diff(diff(Q,'y3'),'p8'))
simplify(diff(diff(Qs,'y3'),'p8'))

simplify(diff(diff(Q,'y4'),'p8'))
simplify(diff(diff(Qs,'y4'),'p8'))

simplify(diff(diff(Q,'y5'),'p8'))
simplify(diff(diff(Qs,'y5'),'p8'))

%  p9

simplify(diff(diff(Q,'y1'),'p9'))
simplify(diff(diff(Qs,'y1'),'p9'))

simplify(diff(diff(Q,'y2'),'p9'))
simplify(diff(diff(Qs,'y2'),'p9'))

simplify(diff(diff(Q,'y3'),'p9'))
simplify(diff(diff(Qs,'y3'),'p9'))

simplify(diff(diff(Q,'y4'),'p9'))
simplify(diff(diff(Qs,'y4'),'p9'))

simplify(diff(diff(Q,'y5'),'p9'))
simplify(diff(diff(Qs,'y5'),'p9'))

simplify(diff(diff(Q,'Q'),'p9'))
simplify(diff(diff(Qs,'Q'),'p9'))

simplify(diff(diff(Q,'Qs'),'p9'))
simplify(diff(diff(Qs,'Qs'),'p9'))

%  p10

simplify(diff(diff(Q,'y1'),'p10'))
simplify(diff(diff(Qs,'y1'),'p10'))

simplify(diff(diff(Q,'y2'),'p10'))
simplify(diff(diff(Qs,'y2'),'p10'))

simplify(diff(diff(Q,'y3'),'p10'))
simplify(diff(diff(Qs,'y3'),'p10'))

simplify(diff(diff(Q,'y4'),'p10'))
simplify(diff(diff(Qs,'y4'),'p10'))

simplify(diff(diff(Q,'y5'),'p10'))
simplify(diff(diff(Qs,'y5'),'p10'))

simplify(diff(diff(Q,'Q'),'p10'))
simplify(diff(diff(Qs,'Q'),'p10'))

simplify(diff(diff(Q,'Qs'),'p10'))
simplify(diff(diff(Qs,'Qs'),'p10'))

%  p11

simplify(diff(diff(Q,'y1'),'p11'))
simplify(diff(diff(Qs,'y1'),'p11'))

simplify(diff(diff(Q,'y2'),'p11'))
simplify(diff(diff(Qs,'y2'),'p11'))

simplify(diff(diff(Q,'y3'),'p11'))
simplify(diff(diff(Qs,'y3'),'p11'))

simplify(diff(diff(Q,'y4'),'p11'))
simplify(diff(diff(Qs,'y4'),'p11'))

simplify(diff(diff(Q,'y5'),'p11'))
simplify(diff(diff(Qs,'y5'),'p11'))

%  p12

simplify(diff(diff(Q,'y1'),'p12'))
simplify(diff(diff(Qs,'y1'),'p12'))

simplify(diff(diff(Q,'y2'),'p12'))
simplify(diff(diff(Qs,'y2'),'p12'))

simplify(diff(diff(Q,'y3'),'p12'))
simplify(diff(diff(Qs,'y3'),'p12'))

simplify(diff(diff(Q,'y4'),'p12'))
simplify(diff(diff(Qs,'y4'),'p12'))

simplify(diff(diff(Q,'y5'),'p12'))
simplify(diff(diff(Qs,'y5'),'p12'))

%  p13

simplify(diff(diff(Q,'y1'),'p13'))
simplify(diff(diff(Qs,'y1'),'p13'))

simplify(diff(diff(Q,'y2'),'p13'))
simplify(diff(diff(Qs,'y2'),'p13'))

simplify(diff(diff(Q,'y3'),'p13'))
simplify(diff(diff(Qs,'y3'),'p13'))

simplify(diff(diff(Q,'y4'),'p13'))
simplify(diff(diff(Qs,'y4'),'p13'))

simplify(diff(diff(Q,'y5'),'p13'))
simplify(diff(diff(Qs,'y5'),'p13'))

%  p14

simplify(diff(diff(Q,'y1'),'p14'))
simplify(diff(diff(Qs,'y1'),'p14'))

simplify(diff(diff(Q,'y2'),'p14'))
simplify(diff(diff(Qs,'y2'),'p14'))

simplify(diff(diff(Q,'y3'),'p14'))
simplify(diff(diff(Qs,'y3'),'p14'))

simplify(diff(diff(Q,'y4'),'p14'))
simplify(diff(diff(Qs,'y4'),'p14'))

simplify(diff(diff(Q,'y5'),'p14'))
simplify(diff(diff(Qs,'y5'),'p14'))

%  p15

simplify(diff(diff(Q,'y1'),'p15'))
simplify(diff(diff(Qs,'y1'),'p15'))

simplify(diff(diff(Q,'y2'),'p15'))
simplify(diff(diff(Qs,'y2'),'p15'))

simplify(diff(diff(Q,'y3'),'p15'))
simplify(diff(diff(Qs,'y3'),'p15'))

simplify(diff(diff(Q,'y4'),'p15'))
simplify(diff(diff(Qs,'y4'),'p15'))

simplify(diff(diff(Q,'y5'),'p15'))
simplify(diff(diff(Qs,'y5'),'p15'))

%  p16

simplify(diff(diff(Q,'y1'),'p16'))
simplify(diff(diff(Qs,'y1'),'p16'))

simplify(diff(diff(Q,'y2'),'p16'))
simplify(diff(diff(Qs,'y2'),'p16'))

simplify(diff(diff(Q,'y3'),'p16'))
simplify(diff(diff(Qs,'y3'),'p16'))

simplify(diff(diff(Q,'y4'),'p16'))
simplify(diff(diff(Qs,'y4'),'p16'))

simplify(diff(diff(Q,'y5'),'p16'))
simplify(diff(diff(Qs,'y5'),'p16'))

