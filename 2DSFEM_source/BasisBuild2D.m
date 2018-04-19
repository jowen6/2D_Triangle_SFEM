function [Phi,Phi1,nodeCoord] = BasisBuild2D(degree,varargin)
%Build Lagrange basis with varargin option to perturb interior triangle
%nodes 
%To perturb: BasisBuild2D(degree,'p',h) 
%h is mesh size

%{
clear all;
degree = 3;
perturb = 'p';
%}

if nargin>1
    perturb = varargin{1};
    h = varargin{2};
else
    perturb = 'X';
end

L1 = [0 , 1; 0, 0]; %x
L2 = [0, 0; 1, 0]; %y
L3 = [1, -1;-1, 0]; %1 - x - y
Identity = [1, 0; 0, 0]; %1

[Bary, nodeCoord] = LagrangeNodes2D(degree);
%The Phi's are indexed like a matrix. The first phi is top left the last
%phi is bottom right
%{
For degree 2 the numbering looks like this:

       Phi{1}O
             |\
             | \
             |  \
       Phi{2}O   O Phi{4}
             |    \
             |     \
       Phi{3}O---O--O Phi{6}
                 ^
                 |
                Phi{5}

%}
Phi = cell(length(Bary),1);
Phi1 = cell(length(Bary),2);
%Build basis functions by multiplying barycentric functions
for i = 1:length(Bary)
    Phi{i} = Identity;
    for j = 1:size(Bary{i},1)
        if Bary{i}(j,1) == 1
            Phi{i} = conv2(Phi{i},L1 - Bary{i}(j,2)*Identity);   
        elseif Bary{i}(j,1) == 2
            Phi{i} = conv2(Phi{i},L2 - Bary{i}(j,2)*Identity); 
        else
            Phi{i} = conv2(Phi{i},L3 - Bary{i}(j,2)*Identity); 
        end
    
    end
    Phi{i} = rot90(Phi{i},2);
    Phi{i} = Phi{i}/polyval2D(Phi{i},nodeCoord(i,1),nodeCoord(i,2));    %Normalize Phi
    [Phi1{i,1},Phi1{i,2}] = polyder2D(Phi{i});                          %Gradient of Phi
end



if (perturb == 'p')     %Build Basis with perturbed interior functions
    
    NumberInteriorNodes = (degree-2)*(degree-1)/2;
    %Random perturbation
    %R = 2*(rand(NumberInteriorNodes,2) -.5)*(h^(degree+1))/(h); %roughly h^{k+1} on unit triangle: (1/4)^degree+1
    R = [0.001*ones(NumberInteriorNodes,1), 0*ones(NumberInteriorNodes),1];
    %R = R*0;
    if (NumberInteriorNodes > 0)
        n = 0;
        index = 1;
        I1 = 0;
        for j = degree-2:-1:1
            n = n+2;
            for k = 1:j
                
                nodeCoord(degree+n+k + I1,1) = nodeCoord(degree+n+k + I1,1) + R(index,1);
                nodeCoord(degree+n+k + I1,2) = nodeCoord(degree+n+k + I1,2) + R(index,2);
                Phi{degree+n+k + I1} = Phi{degree+n+k + I1}/polyval2D(Phi{degree+n+k + I1},nodeCoord(degree+n+k + I1,1),nodeCoord(degree+n+k + I1,2));    %Normalize Phi at perturbed point
                [Phi1{degree+n+k + I1,1},Phi1{degree+n+k + I1,2}] = polyder2D(Phi{degree+n+k + I1});                          %Gradient of Phi at perturbed point
                
                index = index +1;
            end
            I1 = I1+j;
        end
    end
end
%{
x = linspace(0,1,101);
[X,Y] = meshgrid(x,x);
Z = polyval2D(Phi{1},X(:),Y(:));
Z = Z.*(X(:)+Y(:)<1)';
Z = reshape(Z,[size(X,1),size(X,2)]);
[Dx,Dy] = polyder2D(Phi{1});
Z1 = polyval2D(Dy,X(:),Y(:));
Z1 = Z1.*(X(:)+Y(:)<1)';
Z1 = reshape(Z1,[size(X,1),size(X,2)]);

mesh(X,Y,Z)
figure
mesh(X,Y,Z1)
%}