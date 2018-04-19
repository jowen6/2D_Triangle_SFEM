function [HigherOrderElem,HigherOrderNode] = HigherElems2D(degree, elem, node, varargin)
%{
clear all
[node,elem] = spheremesh(1);
degree = 3;
perturb = 'p';
%}

if nargin>3 
    perturb = varargin{1};
else
    perturb = 'X';
end

%Edges length(elem) apart correspond to same triangle
edge = [elem(:,1),elem(:,2);elem(:,2),elem(:,3);elem(:,3),elem(:,1)];

[edgeSortLR,I1] = sort(edge,2);
[edgeSortUD,I2] = sortrows(edgeSortLR); %lexicographic ordering based on row values


%Exploit lack of boundary which implies every edge has a partner. Edges come in pairs 
%Create matrix of interior edge nodes
NewEdgeIndices = length(node)+1:length(node) + (degree-1)*length(edgeSortUD)/2;
NewEdgeIndices = reshape(NewEdgeIndices,[(degree-1),length(edgeSortUD)/2])';

%Concatenate interior edge node matrix with edge vertex nodes
HigherOrderEdge1 = [edgeSortUD(1:2:end,1), NewEdgeIndices, edgeSortUD(1:2:end,2)];


%Create list of new edge node coordinates
NumberInteriorNodes = (degree+1)*(degree+2)/2 - 3*degree; %n*(n+1)/2 total nodes minus number edge nodes
HigherOrderNode = [node;zeros((degree-1)*length(edgeSortUD)/2 + NumberInteriorNodes*length(elem),3)];
Step = repmat((1/degree)*[1:degree-1]',[1,3]);

%Edge Interior Nodes
if (perturb == 'p')
    for i = 1:size(HigherOrderEdge1,1)
        
        
        
        Slope1 = node(HigherOrderEdge1(i,end),:) - node(HigherOrderEdge1(i,1),:);
        Slope = repmat(Slope1,[degree-1,1]);
        
        perturbation = 2*(rand(degree-1,1)-.5)* norm(Slope1)^(degree+1);
        perturbation = repmat(perturbation,[1,3]);
        
        
        
        
        
        Intercept1 = node(HigherOrderEdge1(i,1),:);
        Intercept = repmat(Intercept1,[degree-1,1]);
        HigherOrderNode(length(node) + 1 + (i-1)*(degree-1):length(node) + i*(degree-1),:) = Step.*Slope + Intercept + perturbation.*Slope/norm(Slope);
        
    end
else
    
    for i = 1:size(HigherOrderEdge1,1)
        
        
        
        Slope1 = node(HigherOrderEdge1(i,end),:) - node(HigherOrderEdge1(i,1),:);
        Slope = repmat(Slope1,[degree-1,1]);
        Intercept1 = node(HigherOrderEdge1(i,1),:);
        Intercept = repmat(Intercept1,[degree-1,1]);
        HigherOrderNode(length(node) + 1 + (i-1)*(degree-1):length(node) + i*(degree-1),:) = Step.*Slope + Intercept;
        
    end
end    
    
    

%Map higher order edges back to original ordering as edge 
HigherOrderEdge = zeros(2*length(HigherOrderEdge1),degree+1);
HigherOrderEdge(1:2:2*length(HigherOrderEdge1),:) = HigherOrderEdge1;
HigherOrderEdge(2:2:2*length(HigherOrderEdge1),:) = HigherOrderEdge1;
HigherOrderEdge(I2,:) = HigherOrderEdge;
I3 = find(I1(:,1)==2);
HigherOrderEdge(I3,:) = fliplr(HigherOrderEdge(I3,:));

%Create higher order elem and add triangle interior nodes
HigherOrderElem = [HigherOrderEdge(1:length(elem),1:size(HigherOrderEdge,2)-1), HigherOrderEdge(length(elem)+1:2*length(elem),1:size(HigherOrderEdge,2)-1), HigherOrderEdge(2*length(elem)+1:3*length(elem),1:size(HigherOrderEdge,2)-1)];
NumberInteriorNodes = (degree+1)*(degree+2)/2 - 3*degree; %n*(n+1)/2 total nodes minus number edge nodes


InteriorNodes = reshape([length(node) + (degree-1)*length(edgeSortUD)/2 + 1 : length(node) + (degree-1)*length(edgeSortUD)/2 + NumberInteriorNodes*length(elem)],[NumberInteriorNodes,length(HigherOrderElem)])';


HigherOrderElem = [HigherOrderElem,InteriorNodes];

%Triangle Interior Nodes
if (perturb == 'p')
    if (NumberInteriorNodes > 0)
        for i = 1:size(HigherOrderElem,1)
            Slope1 = node(HigherOrderElem(i,1),:) - node(HigherOrderElem(i,1 + degree),:);
            n = 0;
            index = size(HigherOrderElem,2) - NumberInteriorNodes + 1;
            for j = degree-2:-1:1
                n = n+1;
                for k = 1:j
                    perturbation1 = 2*(rand(2,1)-.5)* norm(Slope1)^(degree+1);
                    B = HigherOrderNode(HigherOrderElem(i,1),:);
                    V1 = (HigherOrderNode(HigherOrderElem(i,n+1),:) - B);
                    V2 = (HigherOrderNode(HigherOrderElem(i,2*degree + 1 + n + k),:) - B);
                    HigherOrderNode(HigherOrderElem(i,index),:) =  B + V1 + V2 ...
                        + perturbation1(1)*V1/norm(V1) + perturbation1(2)*V2/norm(V2);

                    index = index + 1;
                end
            end
        end
    end

else
    if (NumberInteriorNodes > 0)
        for i = 1:size(HigherOrderElem,1)
            n = 0;
            index = size(HigherOrderElem,2) - NumberInteriorNodes + 1;
            for j = degree-2:-1:1
                n = n+1;
                for k = 1:j
                    B = HigherOrderNode(HigherOrderElem(i,1),:);
                    V1 = (HigherOrderNode(HigherOrderElem(i,n+1),:) - B);
                    V2 = (HigherOrderNode(HigherOrderElem(i,2*degree + 1 + n + k),:) - B);
                    HigherOrderNode(HigherOrderElem(i,index),:) =  B + V1 + V2;
                    
                    index = index + 1;
                end
            end
        end
    end
end


%{
%scatter3(HigherOrderNode(length(node)+1:end,1),HigherOrderNode(length(node)+1:end,2),HigherOrderNode(length(node)+1:end,3))
%hold on
%scatter3(HigherOrderNode(end-80+1,1),HigherOrderNode(end-80+1,2),HigherOrderNode(end-80+1,3),'r')
Interior_index = length(node) + (degree-1)*length(edgeSortUD)/2 + 1;
figure
scatter3(HigherOrderNode(length(node)+1:Interior_index-1,1),HigherOrderNode(length(node)+1:Interior_index-1,2),HigherOrderNode(length(node)+1:Interior_index-1,3))
%scatter3(HigherOrderNode(Interior_index:end,1),HigherOrderNode(Interior_index:end,2),HigherOrderNode(Interior_index:end,3),'r')
hold on
trimesh(elem,node(:,1),node(:,2),node(:,3),0*node(:,3),'facealpha',0.5)
%}



