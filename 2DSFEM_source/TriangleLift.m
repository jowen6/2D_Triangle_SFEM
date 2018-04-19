function [Projk1] = TriangleLift(L1,Elem,Node,varargin)
%Input single triangle element, entire node list, and Geometry basis
%functions


%Build Goemetric Approximation
%==========================================================================
    %{
    Projkx = zeros(length(L{1})); %x-component of poly interp of geometry
    Projky = zeros(length(L{1})); %y-component of poly interp of geometry
    Projkz = zeros(length(L{1})); %y-component of poly interp of geometry
    %}

    Projkxs1 = zeros(length(L1{1}));    %s Derivative of x-component of poly interp of geometry
    Projkys1 = zeros(length(L1{1}));    %s Derivative of y-component of poly interp of geometry
    Projkzs1 = zeros(length(L1{1}));    %s Derivative of z-component of poly interp of geometry
    
    Projkxt1 = zeros(length(L1{1}));    %t Derivative of x-component of poly interp of geometry
    Projkyt1 = zeros(length(L1{1}));    %t Derivative of y-component of poly interp of geometry
    Projkzt1 = zeros(length(L1{1}));    %t Derivative of z-component of poly interp of geometry    
    
    
if (nargin>3 && varargin{1} == 'p')    %Build perturbed basis
    GeometryNodes = varargin{2};
    surfacedata = varargin{3};
    %Direction vectors for 2D triangle in 3D
    V1 = Node(Elem(3),:) - Node(Elem(1),:);   
    V2 = Node(Elem(2),:) - Node(Elem(1),:);
    %{
          t  
          ^\
          | \
          |  \
        V2|   \
          |    \
          ------>s
            V1
    %}
    %Parameterization of Linear geometry element by unit triangle in (s,t)-coordinates
    Fx = @(s,t) s*V1(1) + t*V2(1) + Node(Elem(1),1);
    Fy = @(s,t) s*V1(2) + t*V2(2) + Node(Elem(1),2);
    Fz = @(s,t) s*V1(3) + t*V2(3) + Node(Elem(1),3);

    for n = 1:size(L1,1)  %Loop over interpolation points of element
            %scatter3(Fx(GeometryNodes(n,1),GeometryNodes(n,2)),Fy(GeometryNodes(n,1),GeometryNodes(n,2)),Fz(GeometryNodes(n,1),GeometryNodes(n,2)),'r')
            %projection onto surface based on iFEM code
            P = surfacedata.project([Fx(GeometryNodes(n,1),GeometryNodes(n,2)), Fy(GeometryNodes(n,1),GeometryNodes(n,2)), Fz(GeometryNodes(n,1),GeometryNodes(n,2))]);

            %Construct polynomial interpolation of geometry
            %{
            Projkx = Projkx + L{n}*P(1);
            Projky = Projky + L{n}*P(2);
            Projkz = Projkz + L{n}*P(3);
            %}
            
            %Derivatives of interpolation of geometry for building metric
            Projkxs1 = Projkxs1 + L1{n,1}*P(1);
            Projkys1 = Projkys1 + L1{n,1}*P(2);
            Projkzs1 = Projkzs1 + L1{n,1}*P(3);
            
            Projkxt1 = Projkxt1 + L1{n,2}*P(1);
            Projkyt1 = Projkyt1 + L1{n,2}*P(2);            
            Projkzt1 = Projkzt1 + L1{n,2}*P(3); 
    end
    
    
    
else    %Build unperturbed basis
 
    for n = 1:size(L1,1)  %Loop over interpolation points of element

            %Construct polynomial interpolation of geometry
            %{
            Projkx = Projkx + L{n}*Node(Elem(n),1);
            Projky = Projky + L{n}*Node(Elem(n),2);
            Projkz = Projkz + L{n}*Node(Elem(n),3);
            %}
            
            %Derivatives of interpolation of geometry for building metric
            Projkxs1 = Projkxs1 + L1{n,1}*Node(Elem(n),1);
            Projkys1 = Projkys1 + L1{n,1}*Node(Elem(n),2);
            Projkzs1 = Projkzs1 + L1{n,1}*Node(Elem(n),3);
            
            Projkxt1 = Projkxt1 + L1{n,2}*Node(Elem(n),1);
            Projkyt1 = Projkyt1 + L1{n,2}*Node(Elem(n),2);            
            Projkzt1 = Projkzt1 + L1{n,2}*Node(Elem(n),3); 
    end

end    

    Projk1{1,1} = Projkxs1;
    Projk1{2,1} = Projkys1;
    Projk1{3,1} = Projkzs1;
    Projk1{1,2} = Projkxt1;
    Projk1{2,2} = Projkyt1;
    Projk1{3,2} = Projkzt1;
    
    
    