function [dKff_dA, dMff_dA] = truss_derivatives_3d(nele, nnodes, fixity, conn, coord, E, rho, eleid)
%% Book keeping
node_dof = 3;  % DOF per node (u,v)
nele_dof = node_dof*2;  % DOF per element
ndof = node_dof*nnodes; % Total # of DOF in the problem
% ID array/matrix: Labels of DOFs
ID = reshape(1:ndof,node_dof,nnodes);
% LM array/matrix: DOF labels for elements
LM = zeros(nele_dof,nele);
for k=1:nele
    LM(1:node_dof,k) = ID(:,conn(k,1));
    LM(node_dof+1:end,k) = ID(:,conn(k,2));
end
%%  Global Data
dKg_dA = zeros(ndof,ndof); % Global stiffness matrix
dMg_dA = zeros(ndof,ndof); % Global mass matrix
%% Boundary conditions
all_dof = ID(:);   % All DOFS
free_dof = all_dof(isnan(fixity)); % Free DOFs
%% Assemble the global stiffness matrix (Kg) and global load vector Po
% Here we use the stiffness function we have developed
% Loop over elements
k = eleid;
coordI = coord(conn(k,1),:); % I end coords
coordJ = coord(conn(k,2),:); % J end coords
[Ke, Me] = truss_element_3d(E(k),rho(k), coordI,coordJ,1);
dKg_dA(LM(:,k),LM(:,k))= dKg_dA(LM(:,k),LM(:,k)) + Ke;
dMg_dA(LM(:,k),LM(:,k))= dMg_dA(LM(:,k),LM(:,k)) + Me;
dKg_dA = 0.5*(dKg_dA+dKg_dA');
dMg_dA = 0.5*(dMg_dA+dMg_dA');
dKff_dA = dKg_dA(free_dof,free_dof);
dMff_dA = dMg_dA(free_dof,free_dof);
end

%%
function [dKe_dA, dMe_dA ] = truss_element_3d(E,rho, coordI, coordJ, type)
% Calculate length
L = sqrt((coordI(1)-coordJ(1))^2 + (coordI(2)-coordJ(2))^2 +...
   (coordI(3)-coordJ(3))^2 );
% calculate tranformation matrix
lamx = (coordJ(1)-coordI(1))/L;
lamy = (coordJ(2)-coordI(2))/L;
lamz = (coordJ(3)-coordI(3))/L;
T = [lamx, lamy,lamz, 0, 0, 0;
        0, 0, 0, lamx, lamy,lamz];
% 2x2 stiffness matrix in LCS
K_local = (E/L)*[1,-1;
             -1,1];
% Apply tranformations from LCS to GCS
% 4x4 stiffness matrix in GCS
dKe_dA = T'*K_local*T;

% Mass Matrix
mt = L*rho; %% Total mass
if type==1 % consistent
    M_local = mt/6*[2,1;
                    1,2];
else
    % Lumped
    M_local = mt/2*eye(2);
end
dMe_dA = T'*M_local*T;
end
