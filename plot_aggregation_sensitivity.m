function plot_aggregation_sensitivity(ele_vec, dgdx_cdm_uncluster, dgdx_uncluster, agg_type)
    figure()
    hold on
    
    plot(ele_vec, dgdx_uncluster,'-r','linewidth',1.5)
    plot(ele_vec, dgdx_cdm_uncluster,':k','linewidth',1.5)
    %xlabel('All design variable number $e$','Interpreter','latex')
    xlabel('Symmetric design variable number $e_{sym}$','Interpreter','latex')
    if agg_type==1  % p-norm
        ylabel('$d|\!|\mbox{\boldmath $\lambda$}|\!|_p/dx_e$', 'Interpreter','latex')
    elseif agg_type==2  %ks function
        ylabel('$dKS(\mbox{\boldmath $\lambda$})/dx_e$', 'Interpreter','latex')
    elseif agg_type==3  %gen function
        ylabel('$dh(\mbox{\boldmath $\lambda$})/dx_{e_{sym}}$', 'Interpreter','latex')
        %ylabel('$dh(\mbox{\boldmath $\lambda$})/dx_{e_sym}$', 'Interpreter','latex')
    elseif agg_type==4  %gen function of cluster mean 
        %name = ['$', 'd','(\mbox{\boldmath\bar{\Lambda}})', '\bar{\Lambda}', '/','dx_e$'];
        name = ['$', 'df','(\bar{\mathbf{\Lambda}})', '/','dx_e$'];
  
        ylabel(name,'Interpreter','latex')
        
    end
    
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca,'fontname','times')  % Set it to times
    grid on
    legend('Ana. sensitivities','CDM',  'Location', 'northeast')
    
    hold off
    
end