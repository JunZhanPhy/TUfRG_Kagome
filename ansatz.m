%% SC order
fsc1nn=zeros(19,9);
fsc1nn(1,4)=1;fsc1nn(2,4)=-1;fsc1nn(1,2)=-1;fsc1nn(3,2)=1;
fsc1nn(1,8)=1;fsc1nn(6,8)=-1;fsc1nn(1,6)=-1;fsc1nn(7,6)=1;
fsc1nn(1,3)=1;fsc1nn(5,3)=-1;fsc1nn(1,7)=-1;fsc1nn(4,7)=1;
fsc1nn=fsc1nn/sqrt(12);

fsc3nn1=zeros(19,9);fsc3nn2=zeros(19,9);
fsc3nn1(2,1)=1; fsc3nn1(3,1)=-1; fsc3nn1(4,1)=-1; fsc3nn1(5,1)=1;fsc3nn2(6,1)=1; fsc3nn2(7,1)=-1;
fsc3nn1(2,5)=1; fsc3nn1(3,5)=-1; fsc3nn2(4,5)=-1; fsc3nn2(5,5)=1;fsc3nn1(6,5)=1; fsc3nn1(7,5)=-1;
fsc3nn2(2,9)=1; fsc3nn2(3,9)=-1; fsc3nn1(4,9)=-1; fsc3nn1(5,9)=1;fsc3nn1(6,9)=1; fsc3nn1(7,9)=-1;
fsc3nn1=fsc3nn1/sqrt(12);fsc3nn2=fsc3nn2/sqrt(6);

fsc4nn1=zeros(19,9);
fsc4nn1(5,4)=1;fsc4nn1(10,4)=-1;fsc4nn1(10,7)=1;fsc4nn1(3,7)=-1;
fsc4nn1(18,2)=1;fsc4nn1(7,2)=-1;fsc4nn1(2,8)=1;fsc4nn1(18,8)=-1;
fsc4nn1(6,3)=1;fsc4nn1(15,3)=-1;fsc4nn1(15,6)=1;fsc4nn1(4,6)=-1;
fsc4nn1=fsc4nn1/sqrt(12);
fsc4nn2=zeros(19,9);
fsc4nn2(19,4)=1;fsc4nn2(6,4)=-1;fsc4nn2(7,7)=1;fsc4nn2(14,7)=-1;
fsc4nn2(4,2)=1;fsc4nn2(11,2)=-1;fsc4nn2(14,8)=1;fsc4nn2(5,8)=-1;
fsc4nn2(11,3)=1;fsc4nn2(2,3)=-1;fsc4nn2(3,6)=1;fsc4nn2(19,6)=-1;
fsc4nn2=fsc4nn2/sqrt(12);

fsc5nn=zeros(19,9);
fsc5nn(3,4)=1;fsc5nn(8,4)=-1;fsc5nn(2,2)=-1;fsc5nn(9,2)=1;
fsc5nn(7,8)=1;fsc5nn(16,8)=-1;fsc5nn(6,6)=-1;fsc5nn(17,6)=1;
fsc5nn(4,3)=1;fsc5nn(13,3)=-1;fsc5nn(5,7)=-1;fsc5nn(12,7)=1;
fsc5nn=fsc5nn/sqrt(12);

fsc6nn=zeros(19,9);
fsc6nn(10,1)=0; fsc6nn(11,1)=0; fsc6nn(14,1)=1; fsc6nn(15,1)=-1;fsc6nn(18,1)=1; fsc6nn(19,1)=-1;
fsc6nn(10,5)=-1;fsc6nn(11,5)=1; fsc6nn(14,5)=-1;fsc6nn(15,5)=1; fsc6nn(18,5)=0; fsc6nn(19,5)=0;
fsc6nn(10,9)=1; fsc6nn(11,9)=-1;fsc6nn(14,9)=0; fsc6nn(15,9)=0; fsc6nn(18,9)=-1;fsc6nn(19,9)=1;
fsc6nn=fsc6nn/sqrt(12);

%% Charge order
cdw1=zeros(19,9);cdw2=zeros(19,9);
cdw1(1,1)=2;cdw1(1,5)=-1;cdw1(1,9)=-1;
cdw2(1,1)=0;cdw2(1,5)=1;cdw2(1,9)=-1;
cdw1=cdw1/sqrt(6);cdw2=cdw2/sqrt(2);

cbo1nna1=zeros(19,9);cbo1nna2=zeros(19,9);cbo1nna3=zeros(19,9);
cbo1nna1(1,4)=1;cbo1nna1(1,2)=1;cbo1nna1(2,4)=-1;cbo1nna1(3,2)=-1;
cbo1nna2(1,8)=1;cbo1nna2(1,6)=1;cbo1nna2(6,8)=-1;cbo1nna2(7,6)=-1;
cbo1nna3(1,3)=1;cbo1nna3(1,7)=1;cbo1nna3(5,3)=-1;cbo1nna3(4,7)=-1;
cbo1nna1=cbo1nna1/sqrt(4);cbo1nna2=cbo1nna2/sqrt(4);cbo1nna3=cbo1nna3/sqrt(4);

cbo2nna1=zeros(19,9);cbo2nna2=zeros(19,9);cbo2nna3=zeros(19,9);
cbo2nna1(4,4)=1;cbo2nna1(5,2)=-1;cbo2nna1(7,4)=-1;cbo2nna1(6,2)=1;
cbo2nna2(3,8)=1;cbo2nna2(2,6)=-1;cbo2nna2(4,8)=-1;cbo2nna2(5,6)=1;
cbo2nna3(3,3)=1;cbo2nna3(2,7)=-1;cbo2nna3(7,3)=-1;cbo2nna3(6,7)=1;
cbo2nna1=cbo2nna1/sqrt(4);cbo2nna2=cbo2nna2/sqrt(4);cbo2nna3=cbo2nna3/sqrt(4);

lco1nna1=zeros(19,9);lco1nna2=zeros(19,9);lco1nna3=zeros(19,9);
lco1nna1(1,4)=1;lco1nna1(1,2)=-1;lco1nna1(2,4)=-1;lco1nna1(3,2)=1;
lco1nna2(1,8)=1;lco1nna2(1,6)=-1;lco1nna2(6,8)=-1;lco1nna2(7,6)=1;
lco1nna3(1,3)=1;lco1nna3(1,7)=-1;lco1nna3(5,3)=-1;lco1nna3(4,7)=1;
lco1nna1=lco1nna1/sqrt(4);lco1nna2=lco1nna2/sqrt(4);lco1nna3=lco1nna3/sqrt(4);

lco2nna1=zeros(19,9);lco2nna2=zeros(19,9);lco2nna3=zeros(19,9);
lco2nna1(4,4)=1;lco2nna1(5,2)=1;lco2nna1(7,4)=-1;lco2nna1(6,2)=-1;
lco2nna2(3,8)=1;lco2nna2(2,6)=1;lco2nna2(4,8)=-1;lco2nna2(5,6)=-1;
lco2nna3(3,3)=1;lco2nna3(2,7)=1;lco2nna3(7,3)=-1;lco2nna3(6,7)=-1;
lco2nna1=lco2nna1/sqrt(4);lco2nna2=lco2nna2/sqrt(4);lco2nna3=lco2nna3/sqrt(4);

%%
ansatz_order=[cdw1;cdw2;lco2nna1];