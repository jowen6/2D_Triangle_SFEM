function [node] = BallPerturb(node,h,Surf_deg)
%Take in a list of nodes and the h value for the mesh. Perturb nodes in a
%Gaussian distributed random direction by O(h^{k+1}) 

%Only for sphere for now

RandPerturbs = 2*(rand(size(node,1),3) -.5)*h^(Surf_deg+1);%Entire interval

%RandPerturbs = (1.25*rand(size(node,1),3)-0.25)*h^(Surf_deg+1);%Positive bias
%RandPerturbs = (rand(size(node,1),3))*h^(Surf_deg+1);%Positive bias
%RandPerturbs = 2*(round(rand(size(node,1),3))-0.5)*h^(Surf_deg+1); %Binary

%RandPerturbs = rand(size(node,1),1)*h^(Surf_deg+1);

%RandPerturbs = ones(size(node,1),1)*h^(Surf_deg+1);

%max(RandPerturbs)
%RandPerturbs = 3*RandPerturbs; 
%NodePerturbs =repmat(RandPerturbs,[1,3]).* node./repmat(sqrt(sum(node.*node,2)),[1,3]); %points perturbation in radial direction

node = node + RandPerturbs;

%mean(sqrt(sum(node.^2,2)))