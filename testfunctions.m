sqr_ro=2;
sqr_w=500;
h=414.26;
l=721.01;
tet=0;
step=0.01;
end_xi=20;

temp1=Lambda_1(sqr_ro,sqr_w,l,h,tet,step,end_xi);
disp(temp1);


plot(temp1.xi,temp1.Lambda);
legend({'Lambda 1 Rod'});
grid;

plot(temp1.xi,temp1.r);
legend({'r'});
grid;

plot(temp1.xi,temp1.teta);
legend({'teta'});
grid;


temp3=Lambda_34(sqr_ro,sqr_w,l,h,tet,step,end_xi);
disp(temp3);

plot(temp3.xi,temp3.Lambda1);
legend({'Lambda 1 Rod'});
grid;

plot(temp3.xi,temp3.r1);
legend({'r'});
grid;

plot(temp3.xi,temp3.teta1);
legend({'teta'});
grid;

plot(temp3.xi,temp3.Lambda2);
legend({'Lambda 2 Rod'});
grid;

plot(temp3.xi,temp3.r2);
legend({'R'});
grid;

plot(temp3.xi,temp3.teta2);
legend({'teta'});
grid;
