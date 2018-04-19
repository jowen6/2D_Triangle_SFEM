function [nodeCoord] = PerturbInteriorLagrangeNodes(degree,nodeCoord)
%Perturbs the Lagrange nodes on the interior of the triangle
%Nodes on edges are not modified due to the need to have neighboring
%triangles remain compatible

NumberInteriorNodes = (degree-2)*(degree-1)/2;
%Random perturbation
R = 2*(rand(NumberInteriorNodes,2) -.5)*.2;%(h^(degree+1))/(h); %roughly h^{k+1} on unit triangle: (1/4)^degree+1
%R = [0*ones(NumberInteriorNodes,1), -0.1*ones(NumberInteriorNodes,1)];
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
            index = index +1;
        end
        I1 = I1+j;
    end
end
