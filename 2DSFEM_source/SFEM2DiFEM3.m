function[Lambda,N1,N2,N3] = SFEM2DiFEM3(Surf_deg, PDE_deg, Geo, NumRefinements, EigTarget, NumEigs, varargin)
%Dependent on iFEM for generating mesh data, surface projections, and quadrature. 
%Write code to perturb inner mesh points slightly!!!!!
%{
clear all
EigTarget = 2;  %Target Eigenvalue. Finds 3 closest to this.
NumEigs = 3;
Geo = 1;
Surf_deg = 2;
PDE_deg = 2;
NumRefinements = 2;
nargin = 7;
%}
varargin{1}
if nargin > 6;
    perturb = 'p';
else
    perturb = 'X';
end



%Generate Mesh Data
%==========================================================================
[node,elem,surfacedata] = GenerateMesh(Geo,NumRefinements);
N1 = mean(sqrt(sum((node(elem(:,1),:)- node(elem(:,2),:)).^2,2)));
N2 = mean(sqrt(sum((node(elem(:,2),:)- node(elem(:,3),:)).^2,2)));
N3 = mean(sqrt(sum((node(elem(:,1),:)- node(elem(:,3),:)).^2,2)));

%node = RadialPerturb(node,N1,Surf_deg);
%{
if (Surf_deg == 1)
    GeoElem = elem;
    GeoElem = GeoElem(:,[3,1,2]);
    GeoNode = node;
else
    [GeoElem,GeoNode] = HigherElems2D(Surf_deg, elem, node);
    GeoElem = ElemNumberRearrange2D(GeoElem,Surf_deg);   %Rearrange HigherOrder Elem indices to work with BasisBuild2D
end

GeoNode = surfacedata.project(GeoNode); %project geometry nodes onto surface based on iFEM code
%}
%figure
%scatter3(GeoNode(:,1),GeoNode(:,2),GeoNode(:,3),'r')
%trimesh(elem,node(:,1),node(:,2),node(:,3),0*node(:,3),'facealpha',0.5)
%hold on

%trisurf(elem,node(:,1),node(:,2),node(:,3),0*node(:,3))
%==========================================================================


%Generate Higher Order Elements
%==========================================================================
if (PDE_deg == 1)
    HigherOrderElem = elem;
    HigherOrderElem = HigherOrderElem(:,[3,1,2]);
else
    HigherOrderElem = HigherElems2D(PDE_deg, elem, node);
    HigherOrderElem = ElemNumberRearrange2D(HigherOrderElem,PDE_deg);   %Rearrange HigherOrder Elem indices to work with BasisBuild2D
end

NumHigherOrderNodes = max(max(HigherOrderElem))
[Phi,Phi1] = BasisBuild2D(PDE_deg);   %PDE Lagrange Basis
%==========================================================================



%Quadrature
%==========================================================================
%Set Quadrature Rule (Max order: 9)
Quad_degree = min(2*PDE_deg,2*Surf_deg) + 1;
[lamb,Qw] = quadpts(Quad_degree);   %Using iFEM quadrature rules that have max order of 9

%Calculate Quadrature points on reference element in Euclidean coordinates
%from Barycentric coordinates
for p = 1:length(Qw)
    Qp(p,:) = lamb(p,1)*[0,0] + lamb(p,2)*[1,0] + lamb(p,3)*[0,1];
end
%==========================================================================



%Basis Function Evaluation at quadrature points
%==========================================================================
PhiQuads = zeros(PDE_deg+1,length(Qw));
Phis1Quads = zeros(PDE_deg+1,length(Qw));
Phit1Quads = zeros(length(Phi),length(Qw));
for i = 1:length(Phi)
   
    PhiQuads(i,:) = polyval2D(Phi{i},Qp(:,1),Qp(:,2));
    Phis1Quads(i,:) = polyval2D(Phi1{i,1},Qp(:,1),Qp(:,2));
    Phit1Quads(i,:) = polyval2D(Phi1{i,2},Qp(:,1),Qp(:,2));

end
%==========================================================================



%Pre-allocate for sparse matrix construction
%==========================================================================
ntriplets = 0;
Total_Entries = size(HigherOrderElem,2)*length(Qw)*length(Phi)^2;
I = zeros(Total_Entries,1);
J = zeros(Total_Entries,1);
XA = zeros(Total_Entries,1);
XM = zeros(Total_Entries,1);
%==========================================================================



%Loop over elements
%==========================================================================
NumPerElem = length(Qw)*length(Phi)^2;
if strcmp(varargin{1}, 'Perturb_Basis_Locations')
    
    [Bary,GeometryNodes] = LagrangeNodes2D(Surf_deg);
    
    %Perturb Locations of Lagrange Nodes on Unit Cell and Build Basis From Nodes
    %======================================================================
    GeometryNodes = PerturbInteriorLagrangeNodes(Surf_deg,GeometryNodes);
    [L,L1] = CustomBasisBuild2D(Surf_deg,GeometryNodes);
    %======================================================================
    
    %[L,L1,GeometryNodes] = BasisBuild2D(Surf_deg,'p',h);      %Geometry Lagrange Basis
    for j = 1:size(HigherOrderElem,1)
        
        %Build Lift Gradient
        %======================================================================
        Projk1 = TriangleLift(L1,elem(j,:),node,'p',GeometryNodes,surfacedata);
        %======================================================================

        %System Assembly
        %======================================================================
        %{
        [I_elem,J_elem,XA_elem,XM_elem] = AssembleLocal(HigherOrderElem(j,:),PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp);
        I(1 + (j-1)*length(Qw)*length(Phi)^2:j*length(Qw)*length(Phi)^2) = I_elem;
        J(1 + (j-1)*length(Qw)*length(Phi)^2:j*length(Qw)*length(Phi)^2) = J_elem;
        XA(1 + (j-1)*length(Qw)*length(Phi)^2:j*length(Qw)*length(Phi)^2) = XA_elem;
        XM(1 + (j-1)*length(Qw)*length(Phi)^2:j*length(Qw)*length(Phi)^2) = XM_elem;
        %}
        [I(1 + (j-1)*NumPerElem:j*NumPerElem), J(1 + (j-1)*NumPerElem:j*NumPerElem),...
            XA(1 + (j-1)*NumPerElem:j*NumPerElem), XM(1 + (j-1)*NumPerElem:j*NumPerElem)]...
            = AssembleLocal(HigherOrderElem(j,:),PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp);        
        
    end
    
elseif strcmp(varargin{1}, 'Perturb_Mapping_Radially')
    if (Surf_deg == 1)
        GeoElem = elem;
        GeoElem = GeoElem(:,[3,1,2]);
        GeoNode = node;
    else
        [GeoElem,GeoNode] = HigherElems2D(Surf_deg, elem, node);
        GeoElem = ElemNumberRearrange2D(GeoElem,Surf_deg);   %Rearrange HigherOrder Elem indices to work with BasisBuild2D
    end
    GeoNode = surfacedata.project(GeoNode); %project geometry nodes onto surface based on iFEM code
    
    
    
    %Perturb Points off Surface Radially
    %==========================================================================
    GeoNode = RadialPerturb(GeoNode,N1,Surf_deg); %Putting a radial error of h^{k+1} in the lifts radially
    %==========================================================================
    
    
    [L,L1] = BasisBuild2D(Surf_deg);      %Geometry Lagrange Basis    
    for j = 1:size(HigherOrderElem,1)
        %Build Lift Gradient
        %======================================================================
        Projk1 = TriangleLift(L1,GeoElem(j,:),GeoNode);
        %======================================================================
        
        %System Assembly
        %======================================================================
        [I(1 + (j-1)*NumPerElem:j*NumPerElem), J(1 + (j-1)*NumPerElem:j*NumPerElem),...
            XA(1 + (j-1)*NumPerElem:j*NumPerElem), XM(1 + (j-1)*NumPerElem:j*NumPerElem)]...
            = AssembleLocal(HigherOrderElem(j,:),PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp);
        
    end
    
elseif strcmp(varargin{1}, 'Perturb_Mapping_Ball')
    if (Surf_deg == 1)
        GeoElem = elem;
        GeoElem = GeoElem(:,[3,1,2]);
        GeoNode = node;
    else
        [GeoElem,GeoNode] = HigherElems2D(Surf_deg, elem, node);
        GeoElem = ElemNumberRearrange2D(GeoElem,Surf_deg);   %Rearrange HigherOrder Elem indices to work with BasisBuild2D
    end
    GeoNode = surfacedata.project(GeoNode); %project geometry nodes onto surface based on iFEM code
    
    
    %Perturb each point off/on surface in any direction within a ball
    %==========================================================================
    GeoNode = BallPerturb(GeoNode,N1,Surf_deg); %Putting a random error of h^{k+1} at all the lifts in any direction within a ball
    %==========================================================================
    
    
    [L,L1] = BasisBuild2D(Surf_deg);      %Geometry Lagrange Basis    
    for j = 1:size(HigherOrderElem,1)
        %Build Lift Gradient
        %======================================================================
        Projk1 = TriangleLift(L1,GeoElem(j,:),GeoNode);
        %======================================================================
        
        %System Assembly
        %======================================================================
        [I(1 + (j-1)*NumPerElem:j*NumPerElem), J(1 + (j-1)*NumPerElem:j*NumPerElem),...
            XA(1 + (j-1)*NumPerElem:j*NumPerElem), XM(1 + (j-1)*NumPerElem:j*NumPerElem)]...
            = AssembleLocal(HigherOrderElem(j,:),PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp);
        
    end
    
else
    if (Surf_deg == 1)
        GeoElem = elem;
        GeoElem = GeoElem(:,[3,1,2]);
        GeoNode = node;
    else
        [GeoElem,GeoNode] = HigherElems2D(Surf_deg, elem, node);
        GeoElem = ElemNumberRearrange2D(GeoElem,Surf_deg);   %Rearrange HigherOrder Elem indices to work with BasisBuild2D
    end
    GeoNode = surfacedata.project(GeoNode); %project geometry nodes onto surface based on iFEM code

    [L,L1] = BasisBuild2D(Surf_deg);      %Geometry Lagrange Basis    
    for j = 1:size(HigherOrderElem,1)
        %Build Lift Gradient
        %======================================================================
        Projk1 = TriangleLift(L1,GeoElem(j,:),GeoNode);
        %======================================================================
        
        %System Assembly
        %======================================================================
        [I(1 + (j-1)*NumPerElem:j*NumPerElem), J(1 + (j-1)*NumPerElem:j*NumPerElem),...
            XA(1 + (j-1)*NumPerElem:j*NumPerElem), XM(1 + (j-1)*NumPerElem:j*NumPerElem)]...
            = AssembleLocal(HigherOrderElem(j,:),PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp);
        
    end        
end
%==========================================================================



%Create Sparse System Matrices
%==========================================================================
n = max(max(HigherOrderElem));
clearvars -except I J XA XM n EigTarget NumEigs N1 N2 N3 %GeoNode
A = sparse (I,J,XA,n,n) ;
M = sparse (I,J,XM,n,n) ;

%{
See link below for reason local to global mapping is done the way it is:
http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/
%}
%==========================================================================


%Clear variables to free memory (Essential to avoid OUT OF MEMORY)
%==========================================================================
clearvars -except A M EigTarget NumEigs N1 N2 N3 %GeoNode
%==========================================================================

%Calculate Eigenvalues
%==========================================================================
opts.tol = 1e-12;
Lambda = eigs(A,M,NumEigs,EigTarget,opts);%length(A));
%Lambda = eigs(A,M,20,'sm');%length(A));
%==========================================================================


