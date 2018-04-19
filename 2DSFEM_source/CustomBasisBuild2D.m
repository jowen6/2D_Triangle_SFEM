function [Phi,Phi1] = CustomBasisBuild2D(degree,Nodes)
%Creates Lagrange basis functions of a given degree centered at Nodes 

%{
clear all
degree = 3;
[Bary,Nodes] = LagrangeNodes2D(degree);
%}

NumberNodes = (degree+1)*(degree+2)/2;
B = MonomialBasisBuild(degree);
U = zeros(NumberNodes);
for i = 1:NumberNodes
    for j = 1:NumberNodes
        
        U(i,j) = polyval2D(B{j},Nodes(i,1),Nodes(i,2));
        
    end
end
Phi = cell(NumberNodes,1);
Phi1 = cell(NumberNodes,2);
for k = 1:NumberNodes
    b = zeros(NumberNodes,1);
    b(k) = 1;
    x = U\b;
    Phi{k} = zeros(length(B{1}));
    for j =1:NumberNodes
        Phi{k} = Phi{k} + x(j)*B{j}; %Need to make dimensions of B matrices all the same.
    end
    Phi{k} = Phi{k}.*(abs(Phi{k})>1e-13);
    [Phi1{k,1},Phi1{k,2}] = polyder2D(Phi{k});                          %Gradient of Phi    
end