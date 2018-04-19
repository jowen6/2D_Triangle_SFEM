function [B] = MonomialBasisBuild(degree)

%clear all
%degree = 2;

Lx = [0 , 1; 0, 0]; %x
Ly = [0, 0; 1, 0]; %y
Identity = [1, 0; 0, 0]; %1

B = cell(degree*(degree+1)/2,1);
n = 1; 
for k = 0:degree
    for j = 0:k
        I = Identity;
        for i = 0:j
            if i>0
                I = conv2(I,Lx);
            end
        end

        for i = 0:k-j
            if i>0
                I = conv2(I,Ly);
            end
        end
        B{n} = I;
        n = n+1;
    end
end

N = length(B{end});

for i = 1:(degree+1)*(degree+2)/2 - degree-1
    for j = 1:N-length(B{i})
        B{i} = conv2(B{i},Identity);
    end
end

for i = 1:length(B)    
    B{i} = rot90(B{i},2);
end





