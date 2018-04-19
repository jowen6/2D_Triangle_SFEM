function [Bary,nodeCoord] = LagrangeNodes2D(degree)
%The nodeCoord is the set of coordinates for the nodes on the unit
%triangle starting at the top left and ending at the bottom right going down each column. 

%Bary = [Lagrange basis function, corresponding shift]

%L1 = x
%L2 = y
%L3 = 1 - x - y

%clear all
%degree = 3;


n = degree+1;
N = tril(ones(n));
[Irow,Icol] = ind2sub([size(N,1),size(N,2)],find(N)); %matrix indices of nodes

Nx = zeros(n);
Ny = zeros(n);
x = linspace(0,1,n);
for i = 1:n
    Nx(i:end,i) = x(i); 
    Ny(i,1:i) = x(end-i+1);
end

for i = 1:length(Irow)
    
    node = [Irow(i),Icol(i)];
    nodeCoord(i,:) = [Ny(Irow(i),Icol(i)),Nx(Irow(i),Icol(i))];
    %scatter3(Nx(:),Ny(:),N(:))
    N(node(1),node(2)) = 0;

    for j = 1:degree
        RowSums = sum(N,2);
        ColumnSums = sum(N,1);
        for k = 1:n
            DiagonalSums(k) = sum(diag(N,1-k));
        end
        Side_sizes = [RowSums, ColumnSums', DiagonalSums'];
   


        Side_sizes(node(1),1) = 0;
        Side_sizes(node(2),2) = 0;
        Side_sizes(node(1)-node(2) + 1,3) = 0;
    
        [M,I] = max(Side_sizes);
        [MM,J] = max(M);
        
        %Bary = [Barycentric Function, Shift]
        
        if (J == 1)
            Bary{i}(j,:) = [J,(n-I(J))/(n-1)];
            N(I(J),:) = 0;
        elseif (J == 2)
            Bary{i}(j,:) = [J,(I(J)-1)/(n-1)];
            N(:,I(J)) = 0;
        else
            Bary{i}(j,:) = [J,(I(J)-1)/(n-1)];
            N = N - diag(diag(N,1-I(J)),1-I(J));    %zero out subdiagonal
        end
    
        %Bary{i}(j,:) = [J,I(J)];
    end
    
    N = tril(ones(n));
    
end

