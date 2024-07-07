function plot_3d_truss(nele,coord,conn, eid_flag)
if nargin==3
    eid_flag=0;
end

figure();
hold on
for k=1:nele
    X = coord(conn(k,:),1);
    Y = coord(conn(k,:),2);
    Z = coord(conn(k,:),3);
    line(X,Y,Z,'Color','k','LineWidth',2); % Original configuration of truss
    %line(X,Y,Z,'Color',[0 0 0 0.3],'LineWidth',3); % Original configuration of truss
    plot3(X,Y,Z,'o','MarkerSize',7,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
    if eid_flag==1
        cx =mean(X)+0.1;
        cy =mean(Y)+0.1;
        cz =mean(Z)+0.1;
        text(cx,cy,cz,num2str(k,'%d'));
       
        
    end
end
box(gca,'off');
axis(gca,'auto');
grid(gca,'off');
axis equal
az = -45;
el = 20;
view(az, el);
end