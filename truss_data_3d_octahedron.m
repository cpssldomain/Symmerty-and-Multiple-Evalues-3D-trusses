function [nele, nnodes, coord, conn, fixity] = truss_data_3d_octahedron(r)
%% https://en.wikipedia.org/wiki/Octahedron
close all;
clc;
nnodes = 6 + 6;
%nele = 10;
coord = zeros(nnodes,3);
coord(1,:) = r*[1,0,0];
coord(2,:) = r*[-1,0,0];
coord(3,:) = r*[0,1,0];
coord(4,:) = r*[0,-1,0];
coord(5,:) = r*[0,0,1];
coord(6,:) = r*[0,0,-1];
coord(7:end,:)= 1.5*r*coord(1:6,:);
%
%idx = 1:nnodes;
conn = [1,3;
        1,4;
        1,5;
        1,6;
        2,3;
        2,4;
        2,5;
        2,6;
        3,5;
        3,6;
        4,5;
        4,6;];
idx1= 1:6;idx2=idx1+6;
conn2 = [idx1', idx2'];
conn = [conn;conn2];
nele = size(conn,1);
%% Boundary
fixity = NaN*ones(size(coord));
fixity = fixity';
fixity(:,7:end) = 0.0;
end