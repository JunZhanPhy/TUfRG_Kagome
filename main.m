%% runing parameters
runtype=1; %% 1 for cluster run and 0 for local run 
parallelrun=1; %% 1 for parallel run and 0 for non-parallel run
Ncores=36; %% # of cores for paralell run on cluster
divergence_criterion=18; % divergence criterion of RG flow
ledingorderflow_track=1; %% track the 5 leading order flow or not 
ansatzorderflow_track=1; %% track the effective interaction for ansatz order flow or not 

%% model parameters
%% interaction parameters
U=0;V1=0.0;V2=0.2;V3=0;

%% q/k/formfactor mesh parameters 
qmeshtype=3; 
Lq=8; % qmeshtype=3;
% Lq=18; % qmeshtype=1;
ng_qmesh=8;wmax_qmesh=4;wmin_qmesh=1e-3;

kmeshtype=0;
Lk=Lq*30;
ng_kmesh=10;wmax_kmesh=6;wmin_kemsh=1e-3;
% ng_kmesh=8;wmax_kmesh=4;wmin_kmesh=1e-3;

shell=2; %% numbers of formfactor shells
Nf=1+3*shell + 3*shell^2; %Total number of formfactors

%% tignht binding parameters
latticetype='Kagome';
mu=0; %% chemical potential
Norbit=1; %% number of orbit
Nsub=3; %% number of sublattice
Nb=Norbit*Nsub; %% number of band

%% mesh setting
%% q mesh
% % % qvecs=hexBZmesh(ng_qmesh,wmax_qmesh,wmin_qmesh,'qmesh',1);
% % % qvecs=setqmesh(ng_qmesh,wmax_qmesh,wmin_qmesh,0);
% qvecs=setqmesh(ng_qmesh,wmax_qmesh,wmin_qmesh,1);
% wq=qvecs(:,3);qvecs=qvecs(:,1:2);Nq=length(qvecs);
if (qmeshtype==0)
    %% nonuniform q mesh
    % qvecs=hexBZmesh(ng_qmesh,wmax_qmesh,wmin_qmesh,'qmesh',1);
    qvecs=setqmesh(ng_qmesh,wmax_qmesh,wmin_qmesh,0);
    % wq=qvecs(:,3);qvecs=qvecs(:,1:2);Nq=length(qvecs);
    
    epsilon=1e-10;
    pos= (qvecs(:,2)>=0-epsilon)  & (qvecs(:,1)/sqrt(3)-qvecs(:,2)>=-epsilon) & (abs(qvecs(:,2)+sqrt(3)*qvecs(:,1))<=4*pi/sqrt(3)+epsilon);
    qsam=qvecs(pos,:);
    Nsam=length(qsam);
    qvecs=[qvecs(pos,:);qvecs(~pos,:)];
    wq=qvecs(:,3);qvecs=qvecs(:,1:2);Nq=length(qvecs);
    samid=1:Nsam;
    
    figure
    hold on
    scatter(qvecs(:,1),qvecs(:,2),'b.')
    scatter(qsam(:,1),qsam(:,2),'r.')
    plot(4*pi/3*cos(0+linspace(0,2*pi,6+1)),4*pi/3*sin(0+linspace(0,2*pi,6+1)))
    axis equal
    box on
    % ylim([-3.9,3.9])
    axis([-4.5,4.5,-4,4])
    set(gca, 'xTick', -4:2:4);
    set(gca, 'yTick', -4:2:4);
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    
elseif (qmeshtype==1)
    %% uniform q mesh in PZ
    [X1,Y1]=ndgrid(0:1/(Lq):1-1/Lq);
    v1=[0,1]*4*pi/sqrt(3);v2=[sqrt(3)/2,-1/2]*4*pi/sqrt(3);
    qX=X1*v1(1)+Y1*v2(1);qY=X1*v1(2)+Y1*v2(2);
    qvecs=[qX(:),qY(:)];
    Nq=length(qvecs); % Nq=(Lq+1)^2;
    
    epsilon=1e-10;
    pos=qY>=0-epsilon & qX/sqrt(3)-qY>=-epsilon & abs(qY+sqrt(3)*qX)<=4*pi/sqrt(3)+epsilon;
    qsamX=qX(pos);qsamY=qY(pos);
    qsam=[qsamX(:),qsamY(:)];
    Nsam=length(qsam);
    qvecs=[qX(pos),qY(pos);qX(~pos),qY(~pos)];
    samid=1:Nsam;
    wq=ones(Nq,1)/Nq;
    
    figure
    hold on
%     scatter(kX(:),kY(:),'g.')
    scatter(qvecs(:,1),qvecs(:,2),'b.')
    scatter(qsamX,qsamY,'r.')
    plot(4*pi/3*cos(0+linspace(0,2*pi,6+1)),4*pi/3*sin(0+linspace(0,2*pi,6+1)))
    plot([0,0,2*pi,2*pi,0],[0,4*pi/sqrt(3),2*pi/sqrt(3),-2*pi/sqrt(3),0])
    axis equal
    box on
    ylim([-4,8])
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    
elseif (qmeshtype==2)
    %% uniform q mesh in BZ
    [X1,Y1]=meshgrid(-1:1/Lq:1);
    %     v1=[sqrt(3)/2,1/2]*4*pi/3;v2=[sqrt(3)/2,-1/2]*4*pi/3;
    v1=[1,0]*4*pi/3;v2=[1/2,sqrt(3)/2]*4*pi/3;
    X=X1*v1(1)+Y1*v2(1);Y=X1*v1(2)+Y1*v2(2);
    %     pos=(abs(X)<=2*pi/sqrt(3)+1e-14) &  (abs(X/sqrt(3))+abs(Y)<=4*pi/3+1e-14 );
    pos=(abs(Y)<=2*pi/sqrt(3)+1e-14) &  (abs(Y/sqrt(3))+abs(X)<=4*pi/3+1e-14 );
    qvecs=[X(pos),Y(pos)];
    Nq=length(qvecs);
    
    epsilon=1e-10;pos=qvecs(:,2)>=0-epsilon & qvecs(:,1)/sqrt(3)-qvecs(:,2)>=-epsilon;
    qsam=qvecs(pos,:);Nsam=length(qsam);
    qvecs=[qvecs(pos,:);qvecs(~pos,:)];
    samid=1:Nsam;
    
    %%%calculate weights of momentum mesh
    wq=zeros(Nq,1);
    [X1,Y1]=meshgrid(4*pi/3*linspace(-1,1,1001));
%     pos= abs(X1)<=2*pi/sqrt(3)+eps & abs(X1/sqrt(3))+abs(Y1)<=4*pi/3+eps;
    pos= abs(Y1)<=2*pi/sqrt(3)+eps & abs(Y1/sqrt(3))+abs(X1)<=4*pi/3+eps;
    X=X1(pos);Y=Y1(pos);
    umesh=[X,Y];
    for i=1:length(umesh)
        pos=fcp(qvecs,umesh(i,:));
        wq(pos)=wq(pos)+1;
    end
    wq=round(wq/length(umesh),5,'significant');
    wq=wq/sum(wq);
    
    figure
    hold on
    scatter(qvecs(:,1),qvecs(:,2),'b.')
    scatter(qsam(:,1),qsam(:,2),'r.')
%     scatter(umesh(:,1),umesh(:,2),'g.')
%     plot(4*pi/3*cos(pi/2+linspace(0,2*pi,6+1)),4*pi/3*sin(pi/2+linspace(0,2*pi,6+1)))
    plot(4*pi/3*cos(linspace(0,2*pi,6+1)),4*pi/3*sin(linspace(0,2*pi,6+1)))
    axis equal
    axis([-4.5,4.5,-4,4])
    box on
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');  
    
elseif (qmeshtype==3)
    %% nonuniform q mesh in BZ
    [X1,Y1]=meshgrid(linspace(0,1,Lq+1));
    X2=X1(Y1<=X1+eps);Y2=Y1(Y1<=X1+eps);
%     M=[2*pi/sqrt(3),0];K=[2*pi/sqrt(3),2*pi/3];
    M=[pi,pi/sqrt(3)];K=[4*pi/3,0];
%     delta=1e-12*sqrt(3);
    delta=0;v0=[2*delta,delta];
    v1=M-[delta,0];v2=K-M-[0,4*delta];
    X=X2*v1(1)+Y2*v2(1);Y=X2*v1(2)+Y2*v2(2);
    qsam=v0+[X,Y];
    
    qtem=[qsam;mirror(qsam,0)];
    qvecs=[qtem;rot(qtem,pi/3);rot(qtem,pi*2/3);rot(qtem,pi*3/3);rot(qtem,pi*4/3);rot(qtem,pi*5/3)];
    qvecs=norepeat(qvecs);
    Nq=length(qvecs);
    
    epsilon=1e-10;pos=qvecs(:,2)>=0-epsilon & qvecs(:,1)/sqrt(3)-qvecs(:,2)>=-epsilon;
    qsam=qvecs(pos,:);Nsam=length(qsam);
    qvecs=[qvecs(pos,:);qvecs(~pos,:)];samid=1:Nsam;
    
    %%%calculate weights of momentum mesh
    wq=zeros(Nq,1);
    [X1,Y1]=meshgrid(4*pi/3*linspace(-1,1,1001));
%     pos= abs(X1)<=2*pi/sqrt(3)+eps & abs(X1/sqrt(3))+abs(Y1)<=4*pi/3+eps;
    pos= abs(Y1)<=2*pi/sqrt(3)+eps & abs(Y1/sqrt(3))+abs(X1)<=4*pi/3+eps;
    X=X1(pos);Y=Y1(pos);
    umesh=[X,Y];
    for i=1:length(umesh)
        pos=fcp(qvecs,umesh(i,:));
        wq(pos)=wq(pos)+1;
    end
    wq=round(wq/length(umesh),5,'significant');
    wq=wq/sum(wq);

    
    figure
    hold on
    scatter(qvecs(:,1),qvecs(:,2),'b.')
    scatter(qsam(:,1),qsam(:,2),'r.')
%     scatter(umesh(:,1),umesh(:,2),'g.')
%     plot(4*pi/3*cos(pi/2+linspace(0,2*pi,6+1)),4*pi/3*sin(pi/2+linspace(0,2*pi,6+1)))
    plot(4*pi/3*cos(linspace(0,2*pi,6+1)),4*pi/3*sin(linspace(0,2*pi,6+1)))
    axis equal
%     axis([-4,4,-4.5,4.5])
    axis([-4.5,4.5,-4,4])
    box on
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');  
end

%% k mesh
if(kmeshtype==0)
    % kvecs=hexBZmesh(ng_kmesh,wmax_kmesh,wmin_kemsh,'kmesh',1);
    % kvecs=setkmesh(ng_kmesh,wmax_kmesh,wmin_kemsh,0);
    kvecs=setkmesh(ng_kmesh,wmax_kmesh,wmin_kemsh,1,'Kagome',mu);
    % wk=kvecs(:,3);kvecs=kvecs(:,1:2);Nk=length(kvecs);
elseif (kmeshtype==1)
    [X1,Y1]=ndgrid(0+1/Lk:1/Lk:1);
    v1=[0,1]*4*pi/sqrt(3);v2=[sqrt(3)/2,-1/2]*4*pi/sqrt(3);
    kX=X1*v1(1)+Y1*v2(1);kY=X1*v1(2)+Y1*v2(2);
    kvecs=[kX(:),kY(:),ones(Lk^2,1)/Lk^2];

    figure
    hold on
    scatter(kX(:),kY(:),'g.')
    % % scatter(qX(:),qY(:),'b.')
    % scatter(qvecs(:,1),qvecs(:,2),'b.')
    % scatter(qsamX,qsamY,'r.')
    plot(4*pi/3*cos(0+linspace(0,2*pi,6+1)),4*pi/3*sin(0+linspace(0,2*pi,6+1)))
    axis equal
    box on
    ylim([-4,8])
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
end

%% formfactor mesh
% ps=setrmesh(2*shell,0);
sites=formfactor_grid(2*shell);
vec1=[1,0];vec2=[1/2,sqrt(3)/2];
ps=sites*[vec1;vec2];% ps=sites(:,1)*vec1+sites(:,2)*vec2;

%% symmertry transformations on formfactors and kmesh points
Nsym=12;
ksymto=zeros(Nq,Nsym);
for i=1:Nq
    for j=1:6
        %% rotations
        kto=rot(qvecs(i,:),(j-1)*pi/3);
        if (qmeshtype==1)
            kto=backtoPZ(kto);
        end
        ksymto(i,j)=fcp(qvecs,kto);
        %% mirror relections
        % kto=mirror(qvecs(i,:),(j-1)*pi/6);
        kto=mirror(qvecs(i,:),j*pi/6);
        if (qmeshtype==1)
            kto=backtoPZ(kto);
        end
        ksymto(i,j+6)=fcp(qvecs,kto);
    end
end

rsymto=zeros(Nf,Nsym);
for i=1:Nf
    for j=1:6
        %% rotations
        pto=rot(ps(i,:),(j-1)*pi/3);
        rsymto(i,j)=fcp(ps,pto);
        %% mirror relections
        % pto=mirror(ps(i,:),(j-1)*pi/6);
        pto=mirror(ps(i,:),j*pi/6);
        rsymto(i,j+6)=fcp(ps,pto);
    end
end

% Nb=3;
osymto=zeros(Nb,Nsym);
for i=1:Nb
    for j=[1,4]
        osymto(i,j)=i;
    end
    for j=[2,5]
        osymto(i,j)=i-1+(i==1)*3;
    end
    for j=[3,6]
        osymto(i,j)=i+1-(i==3)*3;
    end
    for j=[7,10]
        if i==1
            osymto(i,j)=1;
        elseif i==2
            osymto(i,j)=3;
        else
            osymto(i,j)=2;
        end
    end
    for j=[8,11]
        if i==2
            osymto(i,j)=2;
        elseif i==1
            osymto(i,j)=3;
        else
            osymto(i,j)=1;
        end
    end
    for j=[9,12]
        if i==3
            osymto(i,j)=3;
        elseif i==1
            osymto(i,j)=2;
        else
            osymto(i,j)=1;
        end
    end
end

usymto=zeros(Nb,Nsym);
usymto(1,1)=1;usymto(2,1)=1;usymto(3,1)=1;
usymto(1,2)=5;usymto(2,2)=1;usymto(3,2)=3;
usymto(1,3)=1;usymto(2,3)=1;usymto(3,3)=1;
usymto(1,4)=1;usymto(2,4)=3;usymto(3,4)=5;
usymto(1,5)=1;usymto(2,5)=1;usymto(3,5)=1;
usymto(1,6)=3;usymto(2,6)=5;usymto(3,6)=1;

usymto(1,7)=1;usymto(2,7)=1;usymto(3,7)=1;
usymto(1,8)=5;usymto(2,8)=3;usymto(3,8)=1;
usymto(1,9)=1;usymto(2,9)=1;usymto(3,9)=1;
usymto(1,10)=1;usymto(2,10)=5;usymto(3,10)=3;
usymto(1,11)=1;usymto(2,11)=1;usymto(3,11)=1;
usymto(1,12)=3;usymto(2,12)=1;usymto(3,12)=5;

rsymtoshift=zeros(Nf,Nsym,Nb,Nb);
for i=1:Nf
    for j=1:Nsym
        for o1=1:Nb
            for o2=1:Nb
                pto= ps(rsymto(i,j),:) + ps(usymto(o1,j ),:)- ps(usymto(o2,j ),:) ;
                rsymtoshift(i,j,o1,o2)=fcp(ps,pto);
            end
        end
    end
end

phasefactormesh=zeros(Nq,Nsym,Nb,Nb);
for i=1:Nq
    for j=1:Nsym
        for o1=1:Nb
            for o3=1:Nb
                %% rotations
                phasefactormesh(i,j,o1,o3)=exp(-1j*dot( qvecs(ksymto(i,j),:) ,ps(usymto(o1,j ),:)-ps(usymto(o3,j ),:)));
            end
        end
    end
end

%% Set initial conditions
% fprintf('Set initial conditions\n')
n1sites=[2,3,4,5,6,7];
n2sites=[10,11,14,15,18,19];
n3sites=[8,9,12,13,16,17];

P0=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);C0=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);D0=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);
%% V1
P0(1,1,1,2,1,2,:)= P0(1,1,1,2,1,2,:)+V1;P0(1,1,1,2,2,1,:)= P0(1,1,1,2,2,1,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,2,1,1,2,:)= P0(1,1,2,1,1,2,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,2,1,2,1,:)= P0(1,1,2,1,2,1,:)+V1;
C0(1,1,1,2,1,2,:)= C0(1,1,1,2,1,2,:)+V1;C0(1,1,1,2,2,1,:)= C0(1,1,1,2,2,1,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,2,1,1,2,:)= C0(1,1,2,1,1,2,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,2,1,2,1,:)= C0(1,1,2,1,2,1,:)+V1;
D0(1,1,1,2,1,2,:)= D0(1,1,1,2,1,2,:)+V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);D0(1,1,1,2,2,1,:)= D0(1,1,1,2,2,1,:)-V1;D0(1,1,2,1,1,2,:)= D0(1,1,2,1,1,2,:)-V1;D0(1,1,2,1,2,1,:)= D0(1,1,2,1,2,1,:)+V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);
P0(2,2,1,2,1,2,:)= P0(2,2,1,2,1,2,:)+V1;P0(2,3,1,2,2,1,:)= P0(2,3,1,2,2,1,:)-V1*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);P0(3,2,2,1,1,2,:)= P0(3,2,2,1,1,2,:)-V1*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);P0(3,3,2,1,2,1,:)= P0(3,3,2,1,2,1,:)+V1;
C0(2,2,1,2,1,2,:)= C0(2,2,1,2,1,2,:)+V1;C0(1,1,1,2,2,1,:)= C0(1,1,1,2,2,1,:)-V1*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,2,1,1,2,:)= C0(1,1,2,1,1,2,:)-V1*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);C0(3,3,2,1,2,1,:)= C0(3,3,2,1,2,1,:)+V1;
D0(1,1,1,2,1,2,:)= D0(1,1,1,2,1,2,:)+V1*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);D0(2,2,1,2,2,1,:)= D0(2,2,1,2,2,1,:)-V1;D0(3,3,2,1,1,2,:)= D0(3,3,2,1,1,2,:)-V1;D0(1,1,2,1,2,1,:)= D0(1,1,2,1,2,1,:)+V1*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);

P0(1,1,1,3,1,3,:)= P0(1,1,1,3,1,3,:)+V1;P0(1,1,1,3,3,1,:)= P0(1,1,1,3,3,1,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,3,1,1,3,:)= P0(1,1,3,1,1,3,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,3,1,3,1,:)= P0(1,1,3,1,3,1,:)+V1;
C0(1,1,1,3,1,3,:)= C0(1,1,1,3,1,3,:)+V1;C0(1,1,1,3,3,1,:)= C0(1,1,1,3,3,1,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,1,1,3,:)= C0(1,1,3,1,1,3,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,1,3,1,:)= C0(1,1,3,1,3,1,:)+V1;
D0(1,1,1,3,1,3,:)= D0(1,1,1,3,1,3,:)+V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);D0(1,1,1,3,3,1,:)= D0(1,1,1,3,3,1,:)-V1;D0(1,1,3,1,1,3,:)= D0(1,1,3,1,1,3,:)-V1;D0(1,1,3,1,3,1,:)= D0(1,1,3,1,3,1,:)+V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);
P0(4,4,1,3,1,3,:)= P0(4,4,1,3,1,3,:)+V1;P0(4,5,1,3,3,1,:)= P0(4,5,1,3,3,1,:)-V1*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,4,3,1,1,3,:)= P0(5,4,3,1,1,3,:)-V1*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,5,3,1,3,1,:)= P0(5,5,3,1,3,1,:)+V1;
C0(4,4,1,3,1,3,:)= C0(4,4,1,3,1,3,:)+V1;C0(1,1,1,3,3,1,:)= C0(1,1,1,3,3,1,:)-V1*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,1,1,3,:)= C0(1,1,3,1,1,3,:)-V1*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(5,5,3,1,3,1,:)= C0(5,5,3,1,3,1,:)+V1;
D0(1,1,1,3,1,3,:)= D0(1,1,1,3,1,3,:)+V1*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);D0(4,4,1,3,3,1,:)= D0(4,4,1,3,3,1,:)-V1;D0(5,5,3,1,1,3,:)= D0(5,5,3,1,1,3,:)-V1;D0(1,1,3,1,3,1,:)= D0(1,1,3,1,3,1,:)+V1*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);

P0(1,1,2,3,2,3,:)= P0(1,1,2,3,2,3,:)+V1;P0(1,1,2,3,3,2,:)= P0(1,1,2,3,3,2,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,3,2,2,3,:)= P0(1,1,3,2,2,3,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);P0(1,1,3,2,3,2,:)= P0(1,1,3,2,3,2,:)+V1;
C0(1,1,2,3,2,3,:)= C0(1,1,2,3,2,3,:)+V1;C0(1,1,2,3,3,2,:)= C0(1,1,2,3,3,2,:)-V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,2,2,3,:)= C0(1,1,3,2,2,3,:)-V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,2,3,2,:)= C0(1,1,3,2,3,2,:)+V1;
D0(1,1,2,3,2,3,:)= D0(1,1,2,3,2,3,:)+V1*reshape(exp(1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);D0(1,1,2,3,3,2,:)= D0(1,1,2,3,3,2,:)-V1;D0(1,1,3,2,2,3,:)= D0(1,1,3,2,2,3,:)-V1;D0(1,1,3,2,3,2,:)= D0(1,1,3,2,3,2,:)+V1*reshape(exp(-1j*qvecs*ps(1,:).'),[1,1,1,1,1,1,Nq]);
P0(6,6,2,3,2,3,:)= P0(6,6,2,3,2,3,:)+V1;P0(6,7,2,3,3,2,:)= P0(6,7,2,3,3,2,:)-V1*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);P0(7,6,3,2,2,3,:)= P0(7,6,3,2,2,3,:)-V1*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);P0(7,7,3,2,3,2,:)= P0(7,7,3,2,3,2,:)+V1;
C0(6,6,2,3,2,3,:)= C0(6,6,2,3,2,3,:)+V1;C0(1,1,2,3,3,2,:)= C0(1,1,2,3,3,2,:)-V1*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,2,2,3,:)= C0(1,1,3,2,2,3,:)-V1*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);C0(7,7,3,2,3,2,:)= C0(7,7,3,2,3,2,:)+V1;
D0(1,1,2,3,2,3,:)= D0(1,1,2,3,2,3,:)+V1*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);D0(6,6,2,3,3,2,:)= D0(6,6,2,3,3,2,:)-V1;D0(7,7,3,2,2,3,:)= D0(7,7,3,2,2,3,:)-V1;D0(1,1,3,2,3,2,:)= D0(1,1,3,2,3,2,:)+V1*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);

%% V2
P0(7,7,1,2,1,2,:)= P0(7,7,1,2,1,2,:)+V2;P0(7,6,1,2,2,1,:)= P0(7,6,1,2,2,1,:)-V2*reshape(exp(1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);P0(6,7,2,1,1,2,:)= P0(6,7,2,1,1,2,:)-V2*reshape(exp(-1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);P0(6,6,2,1,2,1,:)= P0(6,6,2,1,2,1,:)+V2;
C0(7,7,1,2,1,2,:)= C0(7,7,1,2,1,2,:)+V2;C0(1,1,1,2,2,1,:)= C0(1,1,1,2,2,1,:)-V2*reshape(exp(1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,2,1,1,2,:)= C0(1,1,2,1,1,2,:)-V2*reshape(exp(-1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);C0(6,6,2,1,2,1,:)= C0(6,6,2,1,2,1,:)+V2;
D0(1,1,1,2,1,2,:)= D0(1,1,1,2,1,2,:)+V2*reshape(exp(1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);D0(7,7,1,2,2,1,:)= D0(7,7,1,2,2,1,:)-V2;D0(6,6,2,1,1,2,:)= D0(6,6,2,1,1,2,:)-V2;D0(1,1,2,1,2,1,:)= D0(1,1,2,1,2,1,:)+V2*reshape(exp(-1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);
P0(4,4,1,2,1,2,:)= P0(4,4,1,2,1,2,:)+V2;P0(4,5,1,2,2,1,:)= P0(4,5,1,2,2,1,:)-V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,4,2,1,1,2,:)= P0(5,4,2,1,1,2,:)-V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,5,2,1,2,1,:)= P0(5,5,2,1,2,1,:)+V2;
C0(4,4,1,2,1,2,:)= C0(4,4,1,2,1,2,:)+V2;C0(1,1,1,2,2,1,:)= C0(1,1,1,2,2,1,:)-V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,2,1,1,2,:)= C0(1,1,2,1,1,2,:)-V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(5,5,2,1,2,1,:)= C0(5,5,2,1,2,1,:)+V2;
D0(1,1,1,2,1,2,:)= D0(1,1,1,2,1,2,:)+V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);D0(4,4,1,2,2,1,:)= D0(4,4,1,2,2,1,:)-V2;D0(5,5,2,1,1,2,:)= D0(5,5,2,1,1,2,:)-V2;D0(1,1,2,1,2,1,:)= D0(1,1,2,1,2,1,:)+V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);
 
P0(6,6,1,3,1,3,:)= P0(6,6,1,3,1,3,:)+V2;P0(6,7,1,3,3,1,:)= P0(6,7,1,3,3,1,:)-V2*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);P0(7,6,3,1,1,3,:)= P0(7,6,3,1,1,3,:)-V2*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);P0(7,7,3,1,3,1,:)= P0(7,7,3,1,3,1,:)+V2;
C0(6,6,1,3,1,3,:)= C0(6,6,1,3,1,3,:)+V2;C0(1,1,1,3,3,1,:)= C0(1,1,1,3,3,1,:)-V2*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,1,1,3,:)= C0(1,1,3,1,1,3,:)-V2*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);C0(7,7,3,1,3,1,:)= C0(7,7,3,1,3,1,:)+V2;
D0(1,1,1,3,1,3,:)= D0(1,1,1,3,1,3,:)+V2*reshape(exp(1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);D0(6,6,1,3,3,1,:)= D0(6,6,1,3,3,1,:)-V2;D0(7,7,3,1,1,3,:)= D0(7,7,3,1,1,3,:)-V2;D0(1,1,3,1,3,1,:)= D0(1,1,3,1,3,1,:)+V2*reshape(exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);
P0(2,2,1,3,1,3,:)= P0(2,2,1,3,1,3,:)+V2;P0(2,3,1,3,3,1,:)= P0(2,3,1,3,3,1,:)-V2*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);P0(3,2,3,1,1,3,:)= P0(3,2,3,1,1,3,:)-V2*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);P0(3,3,3,1,3,1,:)= P0(3,3,3,1,3,1,:)+V2;
C0(2,2,1,3,1,3,:)= C0(2,2,1,3,1,3,:)+V2;C0(1,1,1,3,3,1,:)= C0(1,1,1,3,3,1,:)-V2*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,1,1,3,:)= C0(1,1,3,1,1,3,:)-V2*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);C0(3,3,3,1,3,1,:)= C0(3,3,3,1,3,1,:)+V2;
D0(1,1,1,3,1,3,:)= D0(1,1,1,3,1,3,:)+V2*reshape(exp(1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);D0(2,2,1,3,3,1,:)= D0(2,2,1,3,3,1,:)-V2;D0(3,3,3,1,1,3,:)= D0(3,3,3,1,1,3,:)-V2;D0(1,1,3,1,3,1,:)= D0(1,1,3,1,3,1,:)+V2*reshape(exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);
 
P0(3,3,2,3,2,3,:)= P0(3,3,2,3,2,3,:)+V2;P0(3,2,2,3,3,2,:)= P0(3,2,2,3,3,2,:)-V2*reshape(exp(1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);P0(2,3,3,2,2,3,:)= P0(2,3,3,2,2,3,:)-V2*reshape(exp(-1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);P0(2,2,3,2,3,2,:)= P0(2,2,3,2,3,2,:)+V2;
C0(3,3,2,3,2,3,:)= C0(3,3,2,3,2,3,:)+V2;C0(1,1,2,3,3,2,:)= C0(1,1,2,3,3,2,:)-V2*reshape(exp(1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,2,2,3,:)= C0(1,1,3,2,2,3,:)-V2*reshape(exp(-1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);C0(2,2,3,2,3,2,:)= C0(2,2,3,2,3,2,:)+V2;
D0(1,1,2,3,2,3,:)= D0(1,1,2,3,2,3,:)+V2*reshape(exp(1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);D0(3,3,2,3,3,2,:)= D0(3,3,2,3,3,2,:)-V2;D0(2,2,3,2,2,3,:)= D0(2,2,3,2,2,3,:)-V2;D0(1,1,3,2,3,2,:)= D0(1,1,3,2,3,2,:)+V2*reshape(exp(-1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);
P0(4,4,2,3,2,3,:)= P0(4,4,2,3,2,3,:)+V2;P0(4,5,2,3,3,2,:)= P0(4,5,2,3,3,2,:)-V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,4,3,2,2,3,:)= P0(5,4,3,2,2,3,:)-V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);P0(5,5,3,2,3,2,:)= P0(5,5,3,2,3,2,:)+V2;
C0(4,4,2,3,2,3,:)= C0(4,4,2,3,2,3,:)+V2;C0(1,1,2,3,3,2,:)= C0(1,1,2,3,3,2,:)-V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(1,1,3,2,2,3,:)= C0(1,1,3,2,2,3,:)-V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);C0(5,5,3,2,3,2,:)= C0(5,5,3,2,3,2,:)+V2;
D0(1,1,2,3,2,3,:)= D0(1,1,2,3,2,3,:)+V2*reshape(exp(1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);D0(4,4,2,3,3,2,:)= D0(4,4,2,3,3,2,:)-V2;D0(5,5,3,2,2,3,:)= D0(5,5,3,2,2,3,:)-V2;D0(1,1,3,2,3,2,:)= D0(1,1,3,2,3,2,:)+V2*reshape(exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);

%% V3
tem=reshape( 2*V3*(cos(qvecs(:,1))+2*cos(qvecs(:,1)/2).*cos(sqrt(3)*qvecs(:,2)/2)) , [1,1,1,1,1,1,Nq] );
for io=1:Nb
    C0(1,1,io,io,io,io,:)=C0(1,1,io,io,io,io,:)-tem;D0(1,1,io,io,io,io,:)=D0(1,1,io,io,io,io,:)+tem;
    for i=1:6
        P0(n1sites(i),n1sites(i),io,io,io,io,:)=P0(n1sites(i),n1sites(i),io,io,io,io,:)+V3;
        C0(n1sites(i),n1sites(i),io,io,io,io,:)=C0(n1sites(i),n1sites(i),io,io,io,io,:)+V3;
        D0(n1sites(i),n1sites(i),io,io,io,io,:)=D0(n1sites(i),n1sites(i),io,io,io,io,:)-V3;
    end
    
    P0(2,3,io,io,io,io,:)=P0(2,3,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(2,:).'),[1,1,1,1,1,1,Nq]);
    P0(3,2,io,io,io,io,:)=P0(3,2,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(3,:).'),[1,1,1,1,1,1,Nq]);
    P0(4,5,io,io,io,io,:)=P0(4,5,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(4,:).'),[1,1,1,1,1,1,Nq]);
    P0(5,4,io,io,io,io,:)=P0(5,4,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(5,:).'),[1,1,1,1,1,1,Nq]);
    P0(6,7,io,io,io,io,:)=P0(6,7,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(6,:).'),[1,1,1,1,1,1,Nq]);
    P0(7,6,io,io,io,io,:)=P0(7,6,io,io,io,io,:)- reshape(V3*exp(-1j*qvecs*ps(7,:).'),[1,1,1,1,1,1,Nq]);
end

% fprintf('Initialization Finished. Begin TUfRG Flow!\n')

%% vertex_initialization
P=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);C=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);D=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);

%% flow parameters
iter=0;itermax=400;
lambdamin=1e-6;bandwidth=6;
% vstop=3*bandwidth;% vstop=10*bandwidth;% vstop=18;%vstop=50;
vstop=divergence_criterion;vmax=0;

lambdaarr=zeros(itermax,1);pmaxarr=zeros(itermax,1);cmaxarr=zeros(itermax,1);dmaxarr=zeros(itermax,1);kmaxarr=zeros(itermax,1);
% Gid=fcp(lvecs,[0,0]);% Mid=fcp(lvecs,[pi,pi/sqrt(3)]);% BubblesGamma= zeros(2,itermax);% BubblesM= zeros(2,itermax);

if (ledingorderflow_track==1)
    Gid=fcp([0,0],qvecs);
    Mid=fcp([pi,pi/sqrt(3)],qvecs);

    KevalGflow=zeros(itermax,5);
    KevecGflow=zeros(Nf*Nb^2,5,itermax);
    KevalMflow=zeros(itermax,5);
    KevecMflow=zeros(Nf*Nb^2,5,itermax);

    CevalGflow=zeros(itermax,5);
    CevecGflow=zeros(Nf*Nb^2,5,itermax);
    CevalMflow=zeros(itermax,5);
    CevecMflow=zeros(Nf*Nb^2,5,itermax);
    
    PevalGflow=zeros(itermax,5);
    PevecGflow=zeros(Nf*Nb^2,5,itermax);
    % PevalMflow=zeros(itermax,5);
    % PevecMflow=zeros(Nf*Nb^2,5,itermax);
end

if (ansatzorderflow_track==1)
    ansatz;
    % % % Maid=fcp([-pi,-pi/sqrt(3)],qvecs);Mbid=fcp([pi,-pi/sqrt(3)],qvecs);Mcid=fcp([0,2*pi/sqrt(3)],qvecs);
    Maid=fcp([pi,pi/sqrt(3)],qvecs);Mbid=fcp([pi,-pi/sqrt(3)],qvecs);Mcid=fcp([0,2*pi/sqrt(3)],qvecs);
    trackorder={'cdw1','cbo1nna2','cbo2nna2','lco1nna2','lco2nna2','fsc1nn','fsc3nn1','fsc3nn2','fsc4nn1','fsc4nn2','fsc5nn','fsc6nn'};
    % trackorderff={cdw1,cbo1nna2,cbo2nna2,lco1nna2,lco2nna2,fsc1nn,fsc3nn1,fsc3nn2,fsc5nn,fsc6nn};
    trackorderflow=zeros(itermax,length(trackorder));
end

lambda=1.05*bandwidth;dlambda=0.05*lambda;
% fprintf('Bandwidth:%.2f. Initial scale is now:%.4f\n',bandwidth,1.05*bandwidth);

%% Start flow
% fprintf('Star TUfRG flow with mu=%.1f, U=%.1f, nr=%d, nq=%.0f, nk=%.0f!\n',mu,U,Nf,length(qvecs),length(kvecs))
fprintf('Start TUfRG flow with mu=%.2f, U=%.2f, V1=%.2f, V2=%.2f, Nf=%d, Nq=%.0f, Nk=%.0f!\n',mu,U,V1,V2,Nf,length(qvecs),length(kvecs))
if (parallelrun==1 )
    if (runtype ==1)
        c = parcluster('local');

        slurmid = getenv('SLURM_JOB_ID'); %get the unique job ID as seen by Hydra
        % storage_folder=strcat(getenv('SCRATCH'),'/',getenv('SLURM_JOB_ID'));
        % storage_folder = strcat('/mnt/zfsusers/your_username/matlabtemp/', slurmid);
        % mkdir(storage_folder); %make a temporary folder for this session
        storage_folder = strcat('/public/home/jiangkun/zj/matlabtmp/', slurmid);
        mkdir(storage_folder); %make a temporary folder for this session
        c.JobStorageLocation = storage_folder; %and tell the cluster to use it
        
        c.NumWorkers = Ncores;
        parpool(c, c.NumWorkers);
    elseif (runtype ==0) 
        if ( isempty(gcp('nocreate')) )
            parpool('local');
        end
    end
end

while (iter<itermax)&&(lambda>lambdamin)&&(vmax<=vstop)
    iter=iter+1;
    lambda=lambda-dlambda;
    %     fprintf('Iter=%d. Lambda=%.8f. dLambda=%.8f.\n',iter,lambda,dlambda);
    %     fprintf('Calculate Increment.\n')
    %     fprintf('Calculate Bubbles...\n')
    %     fprintf('Calculate Projections.\n')
    %     fprintf('Fourier transform.\n')
    
    buffer=exp(-1j*qvecs*transpose(ps)).*wq;
    tem1=reshape(P,[Nf^2*Nb^4,Nq]);tem2=reshape(C,[Nf^2*Nb^4,Nq]);tem3=reshape(D,[Nf^2*Nb^4,Nq]);
    fv_P=reshape(tem1*buffer,[Nf,Nf,Nb,Nb,Nb,Nb,length(sites)]);fv_C=reshape(tem2*buffer,[Nf,Nf,Nb,Nb,Nb,Nb,length(sites)]);fv_D=reshape(tem3*buffer,[Nf,Nf,Nb,Nb,Nb,Nb,length(sites)]);
    
    PC=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);PD=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);
    CP=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);CD=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);
    DP=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);DC=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);
    
     % fprintf('Project...\n')
    %% projections
%     tic
    for m=1:Nf
        for n=1:Nf
            for mp=1:Nf
                %% C to P projections
                tem=sites(m,:)+sites(n,:)-sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=sites(n,:)-sites(mp,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        PC(m,n,:,:,:,:,iq)=PC(m,n,:,:,:,:,iq)+ fv_C(mp,np,:,:,:,:,tem2)* exp(1j* (ps(tem2,:)*qvecs(iq,:).') );
                    end
                end
                
                %% D to P projections
                tem=sites(m,:)-sites(n,:)-sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=-sites(n,:)-sites(mp,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        PD(m,n,:,:,:,:,iq)=PD(m,n,:,:,:,:,iq)+ fv_D(mp,np,:,:,:,:,tem2)* exp(-1j* (ps(mp,:)*qvecs(iq,:).') );
                    end
                end
                
                %% P to C projections
                tem=sites(m,:)+sites(n,:)-sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=sites(n,:)-sites(mp,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        CP(m,n,:,:,:,:,iq)=CP(m,n,:,:,:,:,iq)+ fv_P(mp,np,:,:,:,:,tem2)* exp(1j* (ps(tem2,:)*qvecs(iq,:).') );
                    end
                end
                
                %% D to C projections
                tem=-sites(m,:)+sites(n,:)+sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=-sites(m,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        CD(m,n,:,:,:,:,iq)=CD(m,n,:,:,:,:,iq)+ fv_D(mp,np,:,:,:,:,tem2)* exp(-1j* (ps(mp,:)*qvecs(iq,:).') );
                    end
                end
                
                %% P to D projections
                tem=-sites(m,:)-sites(n,:)+sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=-sites(m,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        DP(m,n,:,:,:,:,iq)=DP(m,n,:,:,:,:,iq)+ fv_P(mp,np,:,:,:,:,tem2)* exp(-1j* ( (ps(mp,:)-ps(n,:)) * (qvecs(iq,:).') ) );
                    end
                end
                
                %% C to D projections
                tem=-sites(m,:)+sites(n,:)+sites(mp,:);
                divation=abs(sites-tem);
                epsilon=1e-5;
                np=find(divation(:,1)<epsilon & divation(:,2)<epsilon);
                if (length(np)==1) && (np<=Nf)
                    tem1=-sites(m,:);
                    tem2=fcp(sites,tem1);
                    for iq=1:Nq
                        DC(m,n,:,:,:,:,iq)=DC(m,n,:,:,:,:,iq)+ fv_C(mp,np,:,:,:,:,tem2)* exp(-1j* (ps(mp,:)*qvecs(iq,:).') );
                    end
                end
                
            end
        end
    end
%     toc
    
    VP=P0+P+PC+PD;
    VC=C0+C+CP+CD;
    VD=D0+D+DP+DC;
    
     % fprintf('Calculate Bubbles...\n')
    %% caculate the projections of particle-particle and particle-hole bubbles into form-factor basis
    ph=zeros(Nf*Nb^2,Nf*Nb^2,Nq);
    pp=zeros(Nf*Nb^2,Nf*Nb^2,Nq);
    fr=ps(1:Nf,:);
    % parfor iq=1:Nq
    %     fprintf('%d,',iq)
    %     [ph(:,:,iq),pp(:,:,iq)]=susceptibility(qvecs(iq,:),kvecs,fr,lambda);
    % end
    % fprintf('\n')
    parfor iq=1:Nsam
        [ph(:,:,iq),pp(:,:,iq)]=susceptibility_eig(qvecs(iq,:),kvecs,fr,lambda);
    end
    
    %% TUfRG flow for exchange propogators on main patch
    %     fprintf('Calculate main patch....\n')
    dP=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);dC=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);dD=zeros(Nf,Nf,Nb,Nb,Nb,Nb,Nq);
    for iq=1:Nsam
        VPM=reshape(permute(VP(:,:,:,:,:,:,iq),[1,3,4,2,5,6]),[Nf*Nb^2,Nf*Nb^2]);
        VCM=reshape(permute(VC(:,:,:,:,:,:,iq),[1,3,6,2,5,4]),[Nf*Nb^2,Nf*Nb^2]);
        VDM=reshape(permute(VD(:,:,:,:,:,:,iq),[1,3,5,2,6,4]),[Nf*Nb^2,Nf*Nb^2]);
        Xpp=pp(:,:,iq);Xph=ph(:,:,iq);
        
        %% spinful fermionic TUfRG flow equations
%         dPM=   VPM * Xpp * VPM ;
%         dCM=   VCM * Xph * VCM;
%         dDM=  (VCM-VDM) * Xph * VDM +VDM * Xph * (VCM-VDM);
        %% spinless fermionic TUfRG flow equations
        dPM=1/2*VPM * Xpp * VPM ;
        dCM=    VCM * Xph * VCM;
        dDM=  - VDM * Xph * VDM;
        
        dP(:,:,:,:,:,:,iq)=permute(reshape(dPM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,3,5,6]);
        dC(:,:,:,:,:,:,iq)=permute(reshape(dCM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,6,5,3]);
        dD(:,:,:,:,:,:,iq)=permute(reshape(dDM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,6,3,5]);
    end    
    
    P=P+dP*dlambda;
    C=C+dC*dlambda;
    D=D+dD*dlambda;
    
    %% recover all the exchange propagators by point group symmetry
    for isym=2:Nsym
        for o1=1:Nb
            o1p=osymto(o1,isym);
            for o2=1:Nb
                o2p=osymto(o2,isym);
                for o3=1:Nb
                    o3p=osymto(o3,isym);
                    for o4=1:Nb
                        o4p=osymto(o4,isym);
                        for iq=samid
                            iqp=ksymto(iq,isym);
                            for m=1:Nf
                                % mp=rsymto(m,isym);
                                % mp=fcp(ps,ps(rsymto(m,isym),:)+ ps(usymto(o1,isym ),:)-ps(usymto(o2,isym ),:) );
                                mp=rsymtoshift(m,isym,o1,o2);
                                for n=1:Nf
                                    % np=rsymto(n,isym);
                                    % np=fcp(ps,ps(rsymto(n,isym),:)+ ps(usymto(o3,isym ),:)-ps(usymto(o4,isym ),:) );
                                    np=rsymtoshift(n,isym,o3,o4);
                                    
                                    if (mp<=Nf && np<=Nf )
                                        % P(mp,np,o1p,o2p,o3p,o4p,ilp)=exp(-1j*dot(qvecs(ilp,:),ps(usymto(o1,isym ),:)-ps(usymto(o3,isym ),:)))*P(m,n,o1,o2,o3,o4,iq);
                                        % C(mp,np,o1p,o4p,o3p,o2p,ilp)=exp(-1j*dot(qvecs(ilp,:),ps(usymto(o1,isym ),:)-ps(usymto(o3,isym ),:)))*C(m,n,o1,o4,o3,o2,iq);
                                        % D(mp,np,o1p,o4p,o2p,o3p,ilp)=exp(-1j*dot(qvecs(ilp,:),ps(usymto(o1,isym ),:)-ps(usymto(o3,isym ),:)))*D(m,n,o1,o4,o2,o3,iq);
                                        P(mp,np,o1p,o2p,o3p,o4p,iqp)=phasefactormesh(iq,isym,o1,o3)*P(m,n,o1,o2,o3,o4,iq);
                                        C(mp,np,o1p,o4p,o3p,o2p,iqp)=phasefactormesh(iq,isym,o1,o3)*C(m,n,o1,o4,o3,o2,iq);
                                        D(mp,np,o1p,o4p,o2p,o3p,iqp)=phasefactormesh(iq,isym,o1,o3)*D(m,n,o1,o4,o2,o3,iq);
                                    else
                                        % P(m,n,o1,o2,o3,o4,iq)=0;P(m,n,o1,o2,o3,o4,iqp)=0;
                                        % C(m,n,o1,o4,o3,o2,iq)=0;C(m,n,o1,o2,o3,o4,iqp)=0;
                                        % D(m,n,o1,o4,o2,o3,iq)=0;D(m,n,o1,o2,o3,o4,iqp)=0;
                                        P(m,n,o1,o2,o3,o4,:)=0;C(m,n,o1,o4,o3,o2,:)=0;D(m,n,o1,o4,o2,o3,:)=0;
                                    end
                                    
                                end
                            end
                        end
                        
                    end
                end
            end
        end
    end


    %% PHS
    for iq=1:Nq
        PM=reshape(permute(P(:,:,:,:,:,:,iq),[1,3,4,2,5,6]),[Nf*Nb^2,Nf*Nb^2]);
        CM=reshape(permute(C(:,:,:,:,:,:,iq),[1,3,6,2,5,4]),[Nf*Nb^2,Nf*Nb^2]);
        DM=reshape(permute(D(:,:,:,:,:,:,iq),[1,3,5,2,6,4]),[Nf*Nb^2,Nf*Nb^2]);
        
        PM=(PM+PM')/2;CM=(CM+CM')/2;DM=(DM+DM')/2;
        
        P(:,:,:,:,:,:,iq)=permute(reshape(PM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,3,5,6]);
        C(:,:,:,:,:,:,iq)=permute(reshape(CM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,6,5,3]);
        D(:,:,:,:,:,:,iq)=permute(reshape(DM,[Nf,Nb,Nb,Nf,Nb,Nb]),[1,4,2,6,3,5]);
    end
    
    % K=C-2*D;
    K=reshape(permute(C,[1,3,6,2,5,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]) - 2*reshape(permute(D,[1,3,5,2,6,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]);

    %% update
    vmaxp=max(abs(P(:)));vmaxc=max(abs(C(:)));vmaxd=max(abs(D(:)));vmaxk=max(abs(K(:)));

    % fprintf('Pmax=%.8f, Cmax=%.8f, Dmax=%.8f\n',vmaxp,vmaxc,vmaxd)
    % fprintf('Iter=%d, Lambda=%.16f, dLambda=%.16f, Pmax=%.16f, Cmax=%.16f, Kmax=%.16f\n',iter,lambda,dlambda,vmaxp,vmaxc,vmaxk);
    fprintf('Iter=%d, Lambda=%.16f, dLambda=%.16f, Pmax=%.16f, Cmax=%.16f, Dmax=%.16f, Kmax=%.16f\n',iter,lambda,dlambda,vmaxp,vmaxc,vmaxd,vmaxk);
    lambdaarr(iter)=lambda;
    pmaxarr(iter)=vmaxp;
    cmaxarr(iter)=vmaxc;
    dmaxarr(iter)=vmaxd;
    kmaxarr(iter)=vmaxk;
    %     BubblesGamma(:,iter)= [real(ph(1,1,Gid));real(pp(1,1,Gid))];
    %     BubblesM(:,iter)= [real(ph(1,1,Mid));real(pp(1,1,Mid))];
    
    if vmaxp>vstop
        fprintf("P INSTABILITY!!!\n")
    elseif vmaxc>vstop
        fprintf("C INSTABILITY!!!\n")
    elseif vmaxd>vstop
        fprintf("D INSTABILITY!!!\n")
    elseif lambda<lambdamin
        fprintf("No INSTABILITY until Omega_min!!!\n")
    end
    % vmax=max([vmaxp,vmaxc,vmaxd,vmaxk]);
    vmax=max([vmaxp,vmaxc,vmaxd]);
    
    dlambda=min([0.05*lambda, 1/(2*max(abs(dP(:)))) , 1/(2*max(abs(dC(:)))), 1/(2*max(abs(dD(:))))]);
    %     fprintf('Iter=%d,Lambda=%.8f, Pmax=%.8f, Cmax=%.8f, Dmax=%.8f\n',iter,lambda,vmaxp,vmaxc,vmaxd)

%% charge order
    % spin=0;VMq=reshape(permute(C,[1,3,6,2,5,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]) - spin*2*reshape(permute(D,[1,3,5,2,6,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]);
    spin=0;VKMq=reshape(permute(VC,[1,3,6,2,5,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]) - spin*2*reshape(permute(VD,[1,3,5,2,6,4,7]),[Nf*Nb^2,Nf*Nb^2,Nq]);
    VKMG=VKMq(:,:,Gid);VKMM=VKMq(:,:,Mid);
    VKMG=(VKMG+VKMG')/2;VKMM=(VKMM+VKMM')/2;
    %% SC order
    VPMq=reshape(permute(-VP,[1,3,4,2,5,6,7]),[Nf*Nb^2,Nf*Nb^2,Nq]);
    VPMG=VPMq(:,:,Gid);VPMM=VPMq(:,:,Mid);
    VPMG=(VPMG+VPMG')/2;VPMM=(VPMM+VPMM')/2;

    if (ledingorderflow_track==1)
        [evec,eval]=eig(VKMG);eval=diag(eval);
        KevalGflow(iter,:)=eval(end:-1:end-4).';
        KevecGflow(:,:,iter)=evec(:,end:-1:end-4);

        [evec,eval]=eig(VKMM);eval=diag(eval);
        KevalMflow(iter,:)=eval(end:-1:end-4).';
        KevecMflow(:,:,iter)=evec(:,end:-1:end-4);

        [evec,eval]=eig(VPMG);eval=diag(eval);
        PevalGflow(iter,:)=eval(end:-1:end-4).';
        PevecGflow(:,:,iter)=evec(:,end:-1:end-4);

        [evec,eval]=eig(VPMM);eval=diag(eval);
        % PevalMflow(iter,:)=eval(end:-1:end-4).';
        % PevecMflow(:,:,iter)=evec(:,end:-1:end-4);
    end

    if (ansatzorderflow_track==1)
        % trackorder={'cdw1','cbo1nna2','cbo2nna2','lco1nna2','lco2nna2','fsc1nn','fsc3nn1','fsc3nn2','fsc5nn','fsc6nn'};
        VM=VKMG;
        testorder=cdw1;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,1)=real(tem'*VM*tem);

        VM=VKMM;
        testorder=cbo1nna2;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,2)=real(tem'*VM*tem);
        testorder=cbo2nna2;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,3)=real(tem'*VM*tem);
        testorder=lco1nna2;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,4)=real(tem'*VM*tem);
        testorder=lco2nna2;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,5)=real(tem'*VM*tem);


        VM=VPMG;
        % testorder=fsc1nn;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,6)=real(tem'*VM*tem);
        % testorder=fsc3nn1;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,7)=real(tem'*VM*tem);
        % testorder=fsc3nn2;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,8)=real(tem'*VM*tem);
        % testorder=fsc5nn;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,9)=real(tem'*VM*tem);
        % testorder=fsc6nn;tem=reshape(testorder,[Nf*Nb^2,1]);trackorderflow(iter,10)=real(tem'*VM*tem);
        trackorderflow(iter,6:12)=[ctranspose(reshape(fsc1nn,[Nf*Nb^2,1]))*VM*reshape(fsc1nn,[Nf*Nb^2,1]), ...
            ctranspose(reshape(fsc3nn1,[Nf*Nb^2,1]))*VM*reshape(fsc3nn1,[Nf*Nb^2,1]),ctranspose(reshape(fsc3nn2,[Nf*Nb^2,1]))*VM*reshape(fsc3nn2,[Nf*Nb^2,1]),  ...
            ctranspose(reshape(fsc4nn1,[Nf*Nb^2,1]))*VM*reshape(fsc4nn1,[Nf*Nb^2,1]),ctranspose(reshape(fsc4nn2,[Nf*Nb^2,1]))*VM*reshape(fsc4nn2,[Nf*Nb^2,1]),  ...
            ctranspose(reshape(fsc5nn,[Nf*Nb^2,1]))*VM*reshape(fsc5nn,[Nf*Nb^2,1]),ctranspose(reshape(fsc6nn,[Nf*Nb^2,1]))*VM*reshape(fsc6nn,[Nf*Nb^2,1])];
    end

end
delete(gcp('nocreate'));

%% save results
vmaxflow=[lambdaarr,pmaxarr,cmaxarr,dmaxarr,kmaxarr];
vmaxflow=vmaxflow(1:iter,:);

str=['Kagome_spinless_Tflow_',num2str(Nf),'FFs_qmesh=',num2str(Nq),'_kmesh=',num2str(length(kvecs)),'_Vmax=',num2str(vstop),'_mu=',num2str(mu),'_U=',num2str(U),'_V1=',num2str(V1),'_V2=',num2str(V2),'_V3=',num2str(V3),'.mat'];
if (ledingorderflow_track==0 && ansatzorderflow_track==0)
    save(str,'P','C','D','vmaxflow','qvecs');
elseif (ledingorderflow_track==1 && ansatzorderflow_track==0)
    KevalGflow=KevalGflow(1:iter,:);KevecGflow=KevecGflow(:,:,1:iter);
    KevalMflow=KevalMflow(1:iter,:);KevecMflow=KevecMflow(:,:,1:iter);
    PevalGflow=PevalGflow(1:iter,:);PevecGflow=PevecGflow(:,:,1:iter);
    % PevalMflow=PevalMflow(1:iter,:);PevecMflow=PevecMflow(:,:,1:iter);
    save(str,'P','C','D','vmaxflow','qvecs','KevalGflow','KevecGflow','KevalMflow','KevecMflow','PevalGflow','PevecGflow');
elseif (ledingorderflow_track==0 && ansatzorderflow_track==1)    
    trackorderflow=trackorderflow(1:iter,:);
    save(str,'P','C','D','vmaxflow','qvecs','trackorderflow');
elseif (ledingorderflow_track==1 && ansatzorderflow_track==1) 
    KevalGflow=KevalGflow(1:iter,:);KevecGflow=KevecGflow(:,:,1:iter);
    KevalMflow=KevalMflow(1:iter,:);KevecMflow=KevecMflow(:,:,1:iter);
    PevalGflow=PevalGflow(1:iter,:);PevecGflow=PevecGflow(:,:,1:iter);
    % PevalMflow=PevalMflow(1:iter,:);PevecMflow=PevecMflow(:,:,1:iter);
    trackorderflow=trackorderflow(1:iter,:);
    save(str,'P','C','D','vmaxflow','qvecs','KevalGflow','KevecGflow','KevalMflow','KevecMflow','PevalGflow','PevecGflow','trackorderflow');
end

%% save as h5 files
%%% creat h5files
% filename=['Kagome_Tflow_',num2str(Nf),'FFs_qmesh=',num2str(Lq),'_kmesh=',num2str(Lk),'_mu=',num2str(mu),'_U=',num2str(U),'_V1=',num2str(V1),'_V2=',num2str(V2),'.h5'];
% if (exist(filename,'file'))
%     delete(filename)
% end
% 
% h5create(filename,'/Lambda',[length(vmaxflow),1]);
% h5create(filename,'/pmax',[length(vmaxflow),1]);
% h5create(filename,'/cmax',[length(vmaxflow),1]);
% h5create(filename,'/dmax',[length(vmaxflow),1]);
% h5create(filename,'/kmax',[length(vmaxflow),1]);
% 
% h5create(filename,'/p',size(P));
% h5create(filename,'/c',size(C));
% h5create(filename,'/d',size(D));
% 
% h5create(filename,'/coords',size(qvecs));
% % h5create(filename,'/BubblesGamma',size(BubblesGamma));
% % h5create(filename,'/BubblesM',size(BubblesM));
% % h5create(filename,'/evec',[Nq,4]);
% % h5create(filename,'/eval',[4,1]);
% 
% %%% write h5files
% h5write(filename,'/Lambda',vmaxflow(:,1))
% h5write(filename,'/pmax',vmaxflow(:,2))
% h5write(filename,'/cmax',vmaxflow(:,3))
% h5write(filename,'/dmax',vmaxflow(:,4))
% h5write(filename,'/kmax',vmaxflow(:,5))
% 
% h5write(filename,'/p',real(P))
% h5write(filename,'/c',real(C))
% h5write(filename,'/d',real(D))
% 
% h5write(filename,'/coords',qvecs);
% 
% % h5write(filename,'/BubblesGamma',BubblesGamma);
% % h5write(filename,'/BubblesM',BubblesM);
% % h5write(filename,'/evec',real(evec(:,end-3:end)))
% % h5write(filename,'/eval',eval(end-3:end))

%% function
%% functions
function sites=formfactor_grid(s)
N=1+3*s+3*s^2;
sites=zeros(N,2);
a1=[1,0];a2=[0,1];a3=[-1,1];
b1=[-1,1];b2=[-1,0];b3=[0,-1];
sites(1,:)=[0,0];
c=1;
for i=1:s
    for j=0:i-1
        sites(c+1,:)=i*a1+j*b1;
        sites(c+2,:)=-i*a1-j*b1;
        c=c+2;
    end
    for j=0:i-1
        sites(c+1,:)=i*a2+j*b2;
        sites(c+2,:)=-i*a2-j*b2;
        c=c+2;
    end
    for j=0:i-1
        sites(c+1,:)=i*a3+j*b3;
        sites(c+2,:)=-i*a3-j*b3;
        c=c+2;
    end
end
end

function re=fcp(A,B)
d=(A-B).^2;
d=d(:,1)+d(:,2);
[~,re]=min(d);
end

function LaR=rot(La,theta)
[row,~]=size(La);
% LaR=La;
LaR=zeros(size(La));
for i=1:row
    LaR(i,1)=cos(theta)*La(i,1)-sin(theta)*La(i,2);
    LaR(i,2)=sin(theta)*La(i,1)+cos(theta)*La(i,2);
end
end

function LaR=mirror(La,theta)
[row,~]=size(La);
% s=[1/2,sqrt(3)/2; sqrt(3)/2, -1/2];
R=[cos(theta),-sin(theta); sin(theta),cos(theta)];
s=R*[1,0;0,-1]*transpose(R);
LaR=zeros(size(La));
for i=1:row
    LaR(i,1)=s(1,1)*La(i,1)+s(1,2)*La(i,2);
    LaR(i,2)=s(2,1)*La(i,1)+s(2,2)*La(i,2);
end
end

function re = backtoPZ( M )
[row,~]=size(M);
re=zeros(row,2);
count=0;
L=3;
[X,Y]=meshgrid(-L:L);
mnmesh=[X(:),Y(:)];

for i=1:row
    k=M(i,:);
    kx=k(1);ky=k(2);
    if ( kx>=0-1e-10 && kx< 2*pi-1e-10 && ky >= -kx/sqrt(3) -1e-10 && ky < -kx/sqrt(3) + 4*pi/sqrt(3) - 1e-10  )
        count=count+1;
        re(i,:)=k;
    else
        j=1;
%         while (j<=length(mnmesh)&& j~=((2*L+1)^2+1)/2)
        while (j<=length(mnmesh))
            m=mnmesh(j,1);n=mnmesh(j,2);
            if (m==0 && n==0)
                j=j+1;
                continue
            end
%             fprintf("flag, [%d,%d]\n",m,n);
            k0=k+m*4*pi/sqrt(3)*[sqrt(3)/2,-1/2]+n*4*pi/sqrt(3)*[0,1];
            kx=k0(1);ky=k0(2);
            if ( kx>=0-1e-10 && kx< 2*pi-1e-10 && ky >= -kx/sqrt(3) -1e-10 && ky < -kx/sqrt(3) + 4*pi/sqrt(3) - 1e-10 )
%                 fprintf("flag, [%d,%d]\n",m,n);
                count=count+1;
                re(i,:)=k0;
                
                break
            end
            j=j+1;
        end
        
    end
end
if (count<row)
    fprintf("Not all k points are folded back to 1PZ! #=%d\n",count);
end
if (count>row)
    fprintf("More k points are folded back to 1PZ! #=%d\n", count);
end
end

function re=norepeat(La)
[row,~]=size(La);
% repeat=zeros(row,1);
count=0;
epsilon=1e-5;
big=1e4;
for i=1:row
    for j=i+1:row
        if(abs(La(i,1)-La(j,1))<epsilon && abs(La(i,2)-La(j,2))<epsilon)
            count=count+1;
%             repeat(count)=j;
            La(j,1)=big;
        end
    end
end
re=La((La(:,1))<big,:);
end