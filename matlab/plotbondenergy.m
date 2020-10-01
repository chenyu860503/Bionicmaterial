r0=7.34754;
r=6:0.01:8;
e1=4144620*(r-r0).^2;


r1=3:0.01:5;
x1=272000000;
y1=3.17;
y2=y1./r1;
e2=y2.^12;
e3=y2.^6 ;
e4=4*x1*(e2-e3);




plot(r1,e4,r,e1)
grid on
title('strain rate 10^-2mm')
xlabel('micrometer')
ylabel('E(J)')