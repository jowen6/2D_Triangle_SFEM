function [p] = polyval2D(A,x,y)
%A should be square
for i = 1:length(A)
   
    py(i,:) = polyval(A(i,:),x);  %Evaluate at x coordinates
    
end

for j = 1:length(x)
   
    p(j) = polyval(py(:,j),y(j));
    
end