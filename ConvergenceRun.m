%Shape Options
%{
    case 1  %SphereJ
        
    case 2  %TorusJ
          
    case 3  %ellipsoidJ
    
    case 6  %HeartJ

%}

%Perturbation Options
%{
    Perturb_Basis_Locations:     Perturbs locations of basis nodes on unit cell but leaves interpolation points on surface

    Perturb_Mapping_Radially:    Perturbs interpolation points off surface in normal direction

    Perturb_Mapping_Ball:    Perturbs interpolation points in any direction off surface

    Unperturbed:     No perturbation
 
%}
clear all
EigTarget = 2;  %Target Eigenvalue. Finds 3 closest to this.
NumEigs = 3;    
Geo = 1;
Surf_deg = 2;
PDE_deg = 3;
for NumRefinements = 1:6
    NumRefinements
    [Lambda,N1,N2,N3] = SFEM2DiFEM3(Surf_deg, PDE_deg, Geo, NumRefinements, EigTarget, NumEigs, 'Perturb_Mapping_Radially'); %Main function
    L(:,NumRefinements) = Lambda;
    N11(NumRefinements) = N1;
    N22(NumRefinements) = N2;
    N33(NumRefinements) = N3;
end

h1 = N11;
figure
if (Geo == 1)
    p = polyfit(log(h1(end-3:end)),log(abs(L(1,end-3:end)-EigTarget)),1);
    p
    fsize = 18;
    loglog(h1,abs(L(1,:)-EigTarget),'-o',h1,abs(L(2,:)-EigTarget),'-o',h1,abs(L(3,:)-EigTarget),'-o',h1,h1.^3,'--',h1,h1.^4,'--',h1,h1.^5,'--',h1,h1.^6,'--',h1,h1.^7,'--',h1,h1.^8,'--')
    axis([10^(-2) 1 10^(-7) 1])
    set(gca,'FontSize',fsize)
    
    [hleg1, hobj1]=legend('EigErr1','EigErr2','EigErr3','$h^3$','$h^4$','$h^5$','$h^6$','$h^7$','$h^8$')
    h=legend;
    set(h, 'interpreter', 'latex','Fontsize',12);
    %set(hleg1,'units','normalized','position',[0 0 .38 .1])

else
    loglog(h1(1:end-1),abs(L(1,1:end-1)-L(1,end)),'-o',h1(1:end-1),abs(L(2,1:end-1)-L(2,end)),'-x',h1(1:end-1),abs(L(2,1:end-1)-L(2,end)),'-*',h1(1:end-1),h1(1:end-1).^3,'--',h1(1:end-1),h1(1:end-1).^4,'--',h1(1:end-1),h1(1:end-1).^5,'--',h1(1:end-1),h1(1:end-1).^6,'--',h1(1:end-1),h1(1:end-1).^7,'--',h1(1:end-1),h1(1:end-1).^8,'--')
    
    legend('EigErr1','EigErr2','EigErr3','h^3','h^4','h^5','h^6','h^7','h^8')
end
   
