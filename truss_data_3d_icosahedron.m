function [nele, nnodes, coord, conn, fixity] = truss_data_3d_icosahedron(r)
close all;
clc;
nnodes = 12 + 12; % 20 vertices
nele  = 30 + 12; % 30 edges + 12 extra 
coord = zeros(nnodes,3);
g = (1+sqrt(5))*0.5;
%gi=  1/g;
% edge length = 2.0
coord(1,:) = [ 0, 1, g];
coord(2,:) = [ 0, 1,-g];
coord(3,:) = [ 0,-1, g];
coord(4,:) = [ 0,-1,-g];
coord(5,:) = [ 1, g, 0];
coord(6,:) = [ 1,-g, 0];
coord(7,:) = [-1, g, 0];
coord(8,:) = [-1,-g, 0];
coord(9,:) =  [ g, 0, 1];
coord(10,:) = [-g, 0, 1];
coord(11,:) = [ g, 0,-1];
coord(12,:) = [-g, 0,-1];
coord(13:end,:)= 1.5*coord(1:12,:);
%
%idx = 1:nnodes;
%conn=[];
conn =  [1,3;
         1,5;
         1,7;
         1,9;
         1,10;
         2,5;
         2,7
         2,4;
         2,11;
         2,12
         3,6;
         3,8;
         3,9
         3,10;
         4,6;
         4,8;
         4,11;
         4,12;
         5,7;
         5,9;
         5,11;
         6,8;
         6,9;
         6,11;
         7,10;
         7,12;
         8,10;
         8,12;
         9,11;
         10,12;];
idx1 = 1:12; idx2=idx1+12;
conn2 = [idx1', idx2'];
conn = [conn;conn2];
%nele = size(conn,1);
%% Boundary
%fixity = [];
fixity = NaN*ones(size(coord));
fixity = fixity';
fixity(:,13:end) = 0.0;
end

