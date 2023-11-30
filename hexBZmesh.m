function kv = hexBZmesh( ng,wmax,wmin,meshtype,visualization,latticetype,mu)
b=exp(log(wmax/wmin)/ng);

% kmeshnpoints=0;

kmeshgnp=zeros(ng,1);
kmeshgw=zeros(ng,1);
kmeshgarea=zeros(ng,1);
% kmeshgdp=zeros(ng,1);
kmeshgtri=cell(ng,1);
kmeshgp=cell(ng,1);
kmeshginfo=cell(ng,1);

np=6;
kmeshgnp(1)=np;
kmeshgw(1)=wmax;
kmeshginfo{1}=zeros(np,1);
kmeshgp{1}=zeros(np,2);
kmeshgtri{1}=zeros(3,2,np);

vertex=zeros(3,2);
vertex(1,:)=[0,0];vertex(2,:)=[4*pi/3,0];vertex(3,:)=[2*pi/3,2*pi/sqrt(3)];
kmeshgtri{1}(:,:,1)=vertex;
Stot=area(vertex)*6;kmeshgarea(1)=1/6;

for ip=2:np
    kmeshgtri{1}(:,:,ip)=rot(vertex,(ip-1)*pi/3);
end

for ip=1:np
    kmeshgp{1}(ip,:)=(kmeshgtri{1}(1,:,ip)+kmeshgtri{1}(2,:,ip)+kmeshgtri{1}(3,:,ip))/3;
end

ig=1;w=wmax;
while(w>=wmin)
    icount=0;
    np=kmeshgnp(ig);
    for ip=1:np
        vertex=kmeshgtri{ig}(:,:,ip);
        ps=searchpoints_auto(vertex,14);
        if ( any(abs(ek(ps,meshtype,latticetype,mu))<w) )
            kmeshginfo{ig}(ip)=icount+1;
            icount=icount+4;
        end
    end
    
    if(icount==0)
        w=w/b;
        continue
    end
    if(ig+1>ng)
        kmeshginfo{ig}=zeros(np,1);
        break
    end
    
    kmeshgnp(ig+1)=icount;
    kmeshgw(ig+1)=w;
    
    kmeshginfo{ig+1}=zeros(icount,1);
    kmeshgp{ig+1}=zeros(icount,2);
    kmeshgtri{ig+1}=zeros(3,2,icount);
    
    icount=0;
    for ip=1:np
        if (kmeshginfo{ig}(ip)==0)
            continue
        end
        vertex=kmeshgtri{ig}(:,:,ip);
        vertices=newvertices(vertex);
        for i=1:4
            icount=icount+1;
            kmeshgtri{ig+1}(:,:,icount)=vertices(:,:,i);
            kmeshgp{ig+1}(icount,:)=(vertices(1,:,i)+vertices(2,:,i)+vertices(3,:,i))/3;
        end
    end
    kmeshgarea(ig+1)=area(kmeshgtri{ig+1}(:,:,1))/Stot;
    w=w/b;ig=ig+1;
end

icount=0;
for ig=1:ng
    for ip=1:kmeshgnp(ig)
        if(kmeshginfo{ig}(ip))
            icount=icount+1;
        end
    end
end
kmeshnpoints=icount;

nk=kmeshnpoints;
kv=zeros(nk,3);
nk=0;
for ig=1:ng
    for ip=1:kmeshgnp(ig)
        if (kmeshginfo{ig}(ip)==0)
            nk=nk+1;
            kv(nk,1:2)=kmeshgp{ig}(ip,:);
            kv(nk,3)=kmeshgarea(ig);
        end
    end
end

if (visualization)
    figure
    hold on
    % scatter(qv(:,1),qv(:,2),'b.')
    scatter(kv(:,1),kv(:,2),'b.')
    plot(2*pi/sqrt(3)*cos(pi/2+linspace(0,2*pi,6+1)),2*pi/sqrt(3)*sin(pi/2+linspace(0,2*pi,6+1)))
    plot(4*pi/3*cos(0+linspace(0,2*pi,6+1)),4*pi/3*sin(0+linspace(0,2*pi,6+1)))
    axis equal
    divation=0.3;
    axis([-4*pi/3-divation,4*pi/3+divation,-2*pi/sqrt(3)-divation,2*pi/sqrt(3)+divation])
    box on
    set(gca, 'xTick', [-4*pi/3,0,4*pi/3]);
    set(gca,'XTickLabel',{'$$-4\pi/3$$','0','$$4\pi/3$$'},'ticklabelinterpreter','latex')
    set(gca, 'yTick', [-2*pi/sqrt(3),0,2*pi/sqrt(3)]);
    set(gca,'yTickLabel',{'$$-2\pi/\sqrt3$$','0','$$2\pi/\sqrt3$$'},'ticklabelinterpreter','latex')
    set(gca,'fontSize',15, 'fontname' ,'Times','linewidth' ,1 )
    xlabel('$k_x$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ylabel('$k_y$','Interpreter','latex','fontsize',15,'fontname','times new roman','FontWeight','bold','Color','black');
    ax=gca;ax.Position=ax.Position+[0.02,0,0,0];
end

end

 %% local functions used by main fucntion hexBZmesh
function LaR=rot(La,theta)
[row,~]=size(La);
% LaR=La;
LaR=zeros(size(La));
for i=1:row
    LaR(i,1)=cos(theta)*La(i,1)-sin(theta)*La(i,2);
    LaR(i,2)=sin(theta)*La(i,1)+cos(theta)*La(i,2);
end
end

function ps=searchpoints_auto(vertex,n)
a=vertex(3,:)-vertex(2,:);b=vertex(1,:)-vertex(2,:);
[X,Y]=meshgrid(0:1/n:1);
pos=(X+Y<=1+1e-12);X=X(pos);Y=Y(pos);
ps=vertex(2,:)+X*a+Y*b;
end

function re=ek(k,meshtype,latticetype,mu)

if (meshtype=="kmesh")
    if (latticetype=="Triangular")
        re=-2*(cos(k(:,1)) + 2*cos(k(:,1)/2).* cos(sqrt(3)*k(:,2)/2) ) - mu;
    elseif (latticetype=="Honeycomb")
        re=sqrt(3+2*cos(k(:,1))+4*cos(k(:,1)/2).*cos(sqrt(3)*k(:,2)/2)) -mu ;    
    elseif (latticetype=="Kagome")
        re=-1+sqrt(3+2*cos(k(:,1))+4*cos(k(:,1)/2).*cos(sqrt(3)*k(:,2)/2)) -mu ;
    end
    
elseif (meshtype=="qmesh")
%     a=[[1,0]];b=[0,sqrt(3)/2];
    re=(sin(k*[1,0].')).^2+(sin(k*[1/2,sqrt(3)/2].')).^2+(sin(k*[-1/2,sqrt(3)/2].')).^2 ...
        +(sin(k*[-1,0].')).^2+(sin(k*[-1/2,-sqrt(3)/2].')).^2+(sin(k*[1/2,-sqrt(3)/2].')).^2;
% %     re=(sin(k*[1,0].') .* sin(k*[3/4,sqrt(3)/4].') ).^2+(sin(k*[1/2,sqrt(3)/2].') .* sin(k*[0,sqrt(3)/2].') ).^2+(sin(k*[-1/2,sqrt(3)/2].') .* sin(k*[-3/4,sqrt(3)/4].') ).^2 ...
% %         +(sin(k*[-1,0].') .* sin(k*[-3/4,-sqrt(3)/4].') ).^2+(sin(k*[-1/2,-sqrt(3)/2].') .* sin(k*[0,-sqrt(3)/2].') ).^2+(sin(k*[1/2,-sqrt(3)/2].') .* sin(k*[3/4,-sqrt(3)/4].')).^2;
%     re=abs(sin(k*[1,0].') .* sin(k*[3/4,sqrt(3)/4].') )+abs(sin(k*[1/2,sqrt(3)/2].') .* sin(k*[0,sqrt(3)/2].') )+abs(sin(k*[-1/2,sqrt(3)/2].') .* sin(k*[-3/4,sqrt(3)/4].') ) ...
%         +abs(sin(k*[-1,0].') .* sin(k*[-3/4,-sqrt(3)/4].') )+abs(sin(k*[-1/2,-sqrt(3)/2].') .* sin(k*[0,-sqrt(3)/2].') )+abs(sin(k*[1/2,-sqrt(3)/2].') .* sin(k*[3/4,-sqrt(3)/4].'));
end

end

function re=newvertices(vertex)
re=zeros(3,2,4);
v1=vertex(1,:);v2=vertex(2,:);v3=vertex(3,:);
re(:,:,1)=[v1;(v1+v2)/2;(v1+v3)/2];
re(:,:,2)=[v2;(v2+v1)/2;(v2+v3)/2];
re(:,:,3)=[v3;(v3+v1)/2;(v3+v2)/2];
re(:,:,4)=[(v1+v2)/2;(v2+v3)/2;(v1+v3)/2];
end

function re=area(vertex)
a=vertex(2,:)-vertex(1,:);
b=vertex(3,:)-vertex(1,:);
re=a(1)*b(2)-a(2)*b(1);
re=abs(re)/2;
end
