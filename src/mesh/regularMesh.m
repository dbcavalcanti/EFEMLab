function [Node,ELEM] = regularMesh(Lx,Ly,Nx,Ny)
% -------------------------------------------------------------------------
% This function generates a regular finite element mesh with quadrilateral
% elements.
%
% Input:
%   Lx: length in the x-direction
%   Ly: lenght in the y-direction
%   Nx: number of elements in the x-direction
%   Ny: number of elements in the y-direction
%
% Output:
%   NODE: matrix with the nodes coordinates
%   ELEM: matrix with the elements connectivity
% -------------------------------------------------------------------------

% Generate the NODE matrix
[X,Y]= meshgrid(linspace(0,Lx,Nx+1),linspace(0,Ly,Ny+1));
Node= [reshape(X,numel(X),1) reshape(Y,numel(Y),1)];

% Generate the elements 
k= 0; ELEM= zeros(Nx*Ny,4);
% Numbering elements first along x and then along y
for j=1:Ny
    for i=1:Nx  
        k = k+1;
        n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
        ELEM(k,:) = [n1 n2 n2+1 n1+1];
    end
end

end