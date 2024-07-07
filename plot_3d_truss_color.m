function plot_3d_truss_color(fixity, nele,coord,conn,ng, group,eid_flag)
if nargin==5
    eid_flag = 0;
end

color_array =["r","g","b","c","m","y","k"];

figure();
hold on
nnodes = size(fixity,2);
for j=1:nnodes
    if isnan(fixity(1, j))
        plot3(coord(j,1),coord(j,2),coord(j,3),"o",'MarkerSize',8,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
    else
        plot3(coord(j,1),coord(j,2),coord(j,3),'diamond','MarkerSize',8,'MarkerEdgeColor','r', 'LineWidth', 2)
    end
end

for k=1:nele
    color = [];
    for r=1:ng
        if ~isempty(intersect(group(r).gid,k))
            if (r > numel(color_array))
                break;
            end
            color = color_array(r);
        end
    end
    if isempty(color)
        color = rand(1,3);
    end
    X = coord(conn(k,:),1);
    Y = coord(conn(k,:),2);
    Z = coord(conn(k,:),3);
    
    line(X,Y,Z,'Color',color,'LineWidth',3); % Original configuration of truss

    if eid_flag==1
        cx =mean(X)+0.1;
        cy =mean(Y)+0.1;
        cz =mean(Z)+0.1;
        text(cx,cy,cz,num2str(k,'%d'));
    end
   
end

flag = 1;
if flag==1
    box(gca,'off');
    axis(gca,'auto');
    Ax = gca;
    Ax.XAxis.Visible = 'off';
    Ax.YAxis.Visible = 'off';
    Ax.ZAxis.Visible = 'off';
    Ax.XGrid = 'off';
    Ax.YGrid = 'off';
    Ax.ZGrid = 'off';
    Ax.Color = 'none';
end
grid(gca,'off');
axis equal
az = -45;
el = 20;
view(az, el);
end