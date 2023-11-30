function ps = setrmesh( shell,visualiztion )
sites=formfactor_grid(shell);
vec1=[1,0];vec2=[1/2,sqrt(3)/2];
% ps=sites(:,1)*vec1+sites(:,2)*vec2;
ps=sites*[vec1;vec2];

%% visulization form factor grid
if (visualiztion)
    nnbonds=fdbrrp(ps,0,1.1);
    figure
    hold on
    axis off
    axis equal
    for p=1:length(nnbonds)
        i=nnbonds(p,1);j=nnbonds(p,2);
        plot([ps(i,1),ps(j,1)],[ps(i,2),ps(j,2)],'k','linewidth',1);
    end
    theta_c=0;
    
    sz=200;
    % scatter(ps(:,1),ps(:,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980])
    if (shell==2)
        sz=600;
        plot(1*cos(theta_c+linspace(0,2*pi,6+1)),1*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,0,1])
        plot(2*cos(theta_c+linspace(0,2*pi,6+1)),2*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,123,255]/255)
        scatter(ps(1,1),ps(1,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[159,0,0]/255)
        scatter(ps(2:7,1),ps(2:7,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[200,0,0]/255)
        scatter(ps(8:19,1),ps(8:19,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,25,25]/255)
    elseif (shell==4)
        plot(1*cos(theta_c+linspace(0,2*pi,6+1)),1*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,0,1])
        plot(2*cos(theta_c+linspace(0,2*pi,6+1)),2*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,123,255]/255)
        plot(3*cos(theta_c+linspace(0,2*pi,6+1)),3*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,182,255]/255)
        plot(4*cos(theta_c+linspace(0,2*pi,6+1)),4*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,235,255]/255)
        scatter(ps(1,1),ps(1,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[159,0,0]/255)
        scatter(ps(2:7,1),ps(2:7,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[200,0,0]/255)
        scatter(ps(8:19,1),ps(8:19,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,25,25]/255)
        scatter(ps(20:37,1),ps(20:37,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,105,105]/255)
        scatter(ps(38:61,1),ps(38:61,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,185,185]/255)
    elseif (shell==1)
        plot(1*cos(theta_c+linspace(0,2*pi,6+1)),1*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,0,1])
        scatter(ps(1,1),ps(1,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[159,0,0]/255)
        scatter(ps(2:7,1),ps(2:7,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[200,0,0]/255)
    elseif (shell==3)
        plot(1*cos(theta_c+linspace(0,2*pi,6+1)),1*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,0,1])
        plot(2*cos(theta_c+linspace(0,2*pi,6+1)),2*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,123,255]/255)
        plot(3*cos(theta_c+linspace(0,2*pi,6+1)),3*sin(theta_c+linspace(0,2*pi,6+1)),'linewidth',2,'color',[0,182,255]/255)
        scatter(ps(1,1),ps(1,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[159,0,0]/255)
        scatter(ps(2:7,1),ps(2:7,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[200,0,0]/255)
        scatter(ps(8:19,1),ps(8:19,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,25,25]/255)
        scatter(ps(20:37,1),ps(20:37,2),sz,'MarkerEdgeColor','k','MarkerFaceColor',[255,105,105]/255)
    end
    
    for i=1:length(ps)
        if (shell==2)
            text(ps(i,1),ps(i,2),num2str(i),'FontSize',20,'Color','w','HorizontalAlignment','center')
        else
            text(ps(i,1),ps(i,2),num2str(i),'FontSize',10,'Color','w','HorizontalAlignment','center')
        end
    end
end

end

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

function re=fdbrrp(P,r,rp)
N=length(P);
re=zeros(N*(N-1)/2,2);
count=0;
for i=1:N
    for j=i+1:N
        if( r<sqrt( (P(i,1)-P(j,1))^2+ (P(i,2)-P(j,2))^2) && sqrt( (P(i,1)-P(j,1))^2+ (P(i,2)-P(j,2))^2)<rp )
            count=count+1;
            re(count,:)=[i,j];
        end
    end
end
re=re(1:count,:);
end