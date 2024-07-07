function plot_pts3d(coord)
n = size(coord,1);
figure;
hold on
axis equal
for k=1:n
    plot3(coord(k,1),coord(k,2),coord(k,3),'k*')
    text(coord(k,1)+0.1,coord(k,2)+0.1,coord(k,3)+0.1,num2str(k,'%d'));
end
drawnow;
%axis tight;
%axis off;
end