function [nele, nnodes, coord, conn, fixity] = truss_data_3d_spacetruss_Cnv(N, r, h, imp_flag)
if nargin ==3
    imp_flag = 0; % Imperfection flag
end
nnodes = N*2 + 1;
nele = N + N + N*2;
coord = zeros(nnodes,3);
t = 360/N;
th = t/2;
r2 = 2*r*cosd(th); 
for k =1:N
    coord(k,:) = [r2*cosd(t*(k-1)+th),r2*sind(t*(k-1)+th), 0];
    coord(k+N,:) = [r*cosd(t*(k-1)),r*sind(t*(k-1)), h];
end
coord(end,:) = [0, 0, 2*h];
if imp_flag==1
    %coord(end,:) = [-0.15*h, 0.1*h, 2*h]; % Add imperfections at the top node
    coord(end,:) = [-0.001*h, 0.008*h, 2*h]; % Add imperfections at the top node
end
%
idx1 = 1:N;
idx2 = idx1+N;
idx3 = nnodes*ones(1,N);
conn=[];
temp = [idx2', idx3'];
conn = [conn;temp];
temp = [idx1', [idx1+N]'];
conn = [conn;temp];
temp = [idx1', [idx1+N+1]'];
conn = [conn;temp];
conn(end,2) = N+1;
temp = [idx2', [idx2+1]'];
conn = [conn;temp];
conn(end,2) = N+1;
nele = size(conn,1);
%% Boundary
fixity = [];
fixity = NaN*ones(size(coord));
fixity = fixity';
fixity(:,1:N) = 0.0;
end