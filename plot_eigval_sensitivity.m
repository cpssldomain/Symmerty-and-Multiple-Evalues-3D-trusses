
function plot_eigval_sensitivity(e_ord, ele_vec, dfdx_cdm_uncluster, dfdx_uncluster)
    figure()
    hold on
    
    plot(ele_vec, dfdx_uncluster,'-r','linewidth',1.5)
    plot(ele_vec, dfdx_cdm_uncluster,':k','linewidth',1.5)
    name = ['$', 'd\', 'lambda_','{',num2str(e_ord),'}', '/','dx_e$'];
    xlabel('All design variable number $e$','Interpreter','latex')
    ylabel(name,'Interpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca,'fontname','times')  % Set it to times
    grid on
    legend('Ana. sensitivities','CDM',  'Location', 'northwest')
    
    hold off
    
end


