function kv = setkmesh( ng,wmax,wmin,visualization ,latticetype,mu)
kv=hexBZmesh(ng,wmax,wmin,"kmesh",0,latticetype,mu);
%% symmetrization

%% visualization
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

