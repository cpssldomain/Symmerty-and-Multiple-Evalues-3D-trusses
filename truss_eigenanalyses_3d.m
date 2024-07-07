function [emodes, evals] = truss_eigenanalyses_3d(nele, nnodes, fixity, conn, coord, A, E, rho, type)
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
U = zeros(ndof,1);     % Displacement vector
Kg = zeros(ndof,ndof); % Global stiffness matrix
Mg = zeros(ndof,ndof); % Global mass matrix
%% Boundary conditions
all_dof = ID(:);   % All DOFS
free_dof = all_dof(isnan(fixity)); % Free DOFs
supp_dof = setdiff(all_dof,free_dof);  % Supported DOFs
U(supp_dof) = fixity(supp_dof); % Applied displacement
%P(free_dof) = concen(free_dof);
%% Assemble the global stiffness matrix (Kg) and global load vector Po
% Here we use the stiffness function we have developed
% Loop over elements
for k=1:nele
    coordI = coord(conn(k,1),:); % I end coords
    coordJ = coord(conn(k,2),:); % J end coords
    [Ke, Me] = truss_element_3d(E(k), A(k),rho(k), coordI,coordJ,1);
    Kg(LM(:,k),LM(:,k))= Kg(LM(:,k),LM(:,k)) + Ke;
    Mg(LM(:,k),LM(:,k))= Mg(LM(:,k),LM(:,k)) + Me;
    %Po(LM(:,k))= Po(LM(:,k)) + Poe;
end
Kg = 0.5*(Kg+Kg');
Mg = 0.5*(Mg+Mg');
Kff = Kg(free_dof,free_dof);
Mff = Mg(free_dof,free_dof);
if type==1
    [emodes, evals] = eig(Kff);
else
    [emodes, evals] = eig(Kff, Mff,'chol');
end
% sort results
evals  = diag(evals);
[evals, idx] = sort(evals);
emodes = emodes(:, idx);
end


%%
function [Ke, Me ] = truss_element_3d(E, A, rho, coordI, coordJ, type)
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
K_local = (E*A/L)*[1,-1;
             -1,1];
% Apply tranformations from LCS to GCS
% 4x4 stiffness matrix in GCS
Ke = T'*K_local*T;

% Mass Matrix
mt = L*A*rho; %% Total mass
if type==1 % consistent
    M_local = mt/6*[2,1;
                    1,2];
else
    % Lumped
    M_local = mt/2*eye(2);
end
Me = T'*M_local*T;
end
