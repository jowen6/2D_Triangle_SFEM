fsize = 18;
    loglog(Er3k2Unperturbed(2,:),Er3k2Unperturbed(1,:),'-ok',Er3k2Cperturbed(2,:),Er3k2Cperturbed(1,:),'-sk',Er3k2Operturbed(2,:),Er3k2Operturbed(1,:),'-+k',Er3k2Unperturbed(2,:),Er3k2Unperturbed(2,:).^3,'-.k',Er3k2Unperturbed(2,:),Er3k2Unperturbed(2,:).^4,'--k')
    axis([.2*10^(-1) 1 10^(-8) 1])
    set(gca,'FontSize',fsize)
    
    [hleg1, hobj1]=legend('Unperturbed','Centered Perturbed','Biased Perturbed','$h^3$','$h^4$');
    h=legend;
    set(h, 'interpreter', 'latex','Fontsize',12);
    
    xlabel('$h$','interpreter', 'latex','Fontsize',fsize);
    ylabel('$|\lambda-\Lambda|$','interpreter', 'latex','Fontsize',fsize);