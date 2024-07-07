function [nele, nnodes, coord, conn, fixity] = truss_data_3d_dodecahedron(r)
close all;
clc;
nnodes = 20 + 20; % 20 vertices
nele  = 30 + 20; % 30 edges + 20 extra 
coord = zeros(nnodes,3);
g = (1+sqrt(5))*0.5;
gi=  1/g;
coord(1,:) = [ 1, 1, 1];
coord(2,:) = [-1, 1, 1];
coord(3,:) = [ 1,-1, 1];
coord(4,:) = [ 1, 1,-1];
coord(5,:) = [-1,-1, 1];
coord(6,:) = [ 1,-1,-1];
coord(7,:) = [-1, 1,-1];
coord(8,:) = [-1,-1,-1];
coord(9,:) =  [ 0, g, gi];
coord(10,:) = [ 0,-g, gi];
coord(11,:) = [ 0, g,-gi];
coord(12,:) = [ 0,-g,-gi];
coord(13,:) = [ gi, 0, g];
coord(14,:) = [ gi, 0,-g];
coord(15,:) = [-gi, 0, g];
coord(16,:) = [-gi, 0,-g];
coord(17,:) = [ g, gi, 0,];
coord(18,:) = [ g,-gi, 0,];
coord(19,:) = [-g, gi, 0,];
coord(20,:) = [-g,-gi, 0,];
coord(21:end,:)= 1.5*coord(1:20,:);
%
%idx = 1:nnodes;
%conn=[];
conn =  [1,9;
         1,17;
         1,13;
         2,9;
         2,15;
         2,19;
         3,10
         3,18;
         3,13;
         4,11
         4,14;
         4,17;
         5,10
         5,20;
         5,15;
         6,12
         6,14;
         6,18;
         7,11
         7,16;
         7,19;
         8,12
         8,16;
         8,20;
         9,11;
         10,12;
         13,15;
         14,16;
         17,18;
         19,20;];
idx1= 1:20;idx2=idx1+20;
conn2 = [idx1', idx2'];
conn = [conn;conn2];
%nele = size(conn,1);
%% Boundary
%fixity = [];
fixity = NaN*ones(size(coord));
fixity = fixity';
fixity(:,21:40) = 0.0;
end