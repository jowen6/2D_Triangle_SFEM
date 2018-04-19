function [Dx,Dy] = polyder2D(A)
%Differentiate a 2D matrix representation of the basis function.
%clear all
%A = ones(3);
%A = [3 , 0, 0; 0, 0, 1; 0, 0, 0];   %If A is n x n matrix, then A_{i,j} is coefficient of x^{n-i}y^{n-j}
der = [length(A)-1:-1:0];
[Dery,Derx] = meshgrid(der,der);    %Derivative coefficients for a polynomial
Dx = circshift(A.*Derx,1,1);    
Dy = circshift(A.*Dery,1,2);