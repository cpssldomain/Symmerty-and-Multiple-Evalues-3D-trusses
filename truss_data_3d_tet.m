function [nele, nnodes, coord, conn, fixity] = truss_data_3d_tet(r,imp_flag)
if nargin ==1
    imp_flag = 0; % Imperfection flag
end
nnodes = 4 + 4;
nele = 10;
coord = zeros(nnodes,3);
a= sqrt(8/9);
b= sqrt(2/3);
c= sqrt(2/9);
n1 = [a,0,-1/3]; n1 = n1/norm(n1);
n2 = [-c,b,-1/3]; n2 = n2/norm(n2);
n3 = [-c,-b,-1/3]; n3 = n3/norm(n3);
n4 = [0,0,1];
coord(1,:) = r*n1;
coord(2,:) = r*n2;
coord(3,:) = r*n3;
coord(4,:) = r*n4;
coord(5,:) = 2*r*n1;
coord(6,:) = 2*r*n2;
coord(7,:) = 2*r*n3;
coord(8,:) = 2*r*n4;

if imp_flag==1
    impvec = [1.1,-1.2,1.25]; 
    coord(1,:) = impvec.*coord(1,:); % Add imperfections at one node
end
%
%idx = 1:nnodes;
conn = [1,2;
        2,3;
        3,1
        1,4;
        2,4;
        3,4;
        1,5
        2,6;
        3,7;
        4,8;];
%% Boundary
fixity = NaN*ones(size(coord));
fixity = fixity';
fixity(:,5:8) = 0.0;
end