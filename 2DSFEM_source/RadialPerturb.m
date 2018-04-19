function [node] = RadialPerturb(node,h,Surf_deg)
%Take in a list of nodes and the h value for the mesh. Perturb nodes along
%the normal to the surface by O(h^{k+1}) randomly

%Only for sphere for now
%RandPerturbs = 2*(rand(size(node,1),1) -.5)*h^(Surf_deg+1);
%RandPerturbs = rand(size(node,1),1)*h^(Surf_deg+1);
%RandPerturbs = 2*(round(rand(size(node,1),1))-0.5)*h^(Surf_deg+1);
RandPerturbs = (2*rand(size(node,1),1)-0.5)*h^(Surf_deg+1);
%RandPerturbs = ones(size(node,1),1)*h^(Surf_deg+1);
%max(RandPerturbs)
%RandPerturbs = 3*RandPerturbs; 
NodePerturbs =repmat(RandPerturbs,[1,3]).* node./repmat(sqrt(sum(node.*node,2)),[1,3]);

node = node + NodePerturbs;

%mean(sqrt(sum(node.^2,2)))