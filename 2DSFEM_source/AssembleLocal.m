function [I_elem,J_elem,XA_elem,XM_elem] = AssembleLocal(HigherOrderElem,PhiQuads,Phis1Quads,Phit1Quads,Projk1,Qw,Qp)
    %System Assembly
    %======================================================================
    %Build the Riemannian Metric evaluated at quadrature points
    dChi1 = [polyval2D(Projk1{1,1},Qp(:,1),Qp(:,2));polyval2D(Projk1{1,2},Qp(:,1),Qp(:,2))];
    dChi2 = [polyval2D(Projk1{2,1},Qp(:,1),Qp(:,2));polyval2D(Projk1{2,2},Qp(:,1),Qp(:,2))];
    dChi3 = [polyval2D(Projk1{3,1},Qp(:,1),Qp(:,2));polyval2D(Projk1{3,2},Qp(:,1),Qp(:,2))];
    Total_Entries = length(Qw)*size(PhiQuads,1)^2;
    I_elem = zeros(Total_Entries,1);
    J_elem = zeros(Total_Entries,1);
    XA_elem = zeros(Total_Entries,1);
    XM_elem = zeros(Total_Entries,1);
    %Construct Local Stiffness and Mass matrix
    ntriplets = 0;
    for m = 1:length(Qw)
        Gc = [dChi1(:,m),dChi2(:,m),dChi3(:,m)]*[dChi1(:,m),dChi2(:,m),dChi3(:,m)]'; %Riemannian metric
        GcI = inv(Gc);
        W = sqrt(det(Gc));
        for i = 1:size(PhiQuads,1)
            for k = 1:size(PhiQuads,1)

                GcInverse = GcI*[Phis1Quads(k,m),Phit1Quads(k,m)]';

                ntriplets = ntriplets + 1;
                I_elem(ntriplets) = HigherOrderElem(i);
                J_elem(ntriplets) = HigherOrderElem(k);
                XA_elem(ntriplets) = Qw(m)*([Phis1Quads(i,m),Phit1Quads(i,m)]*GcInverse*W);
                XM_elem(ntriplets) = Qw(m)*PhiQuads(i,m)*PhiQuads(k,m)*W;
            end
        end
    end
    %======================================================================
