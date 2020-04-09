


function [] = save_to_tif(length, width, name, laberFontSize, axisFontSize)

    figWidth = width;

    figHeight = length;

    set(gca,'FontName','Helvetica','FontSize', laberFontSize, 'LineWidth', 1);
    
    ax = gca;
    
    ax.FontSize = axisFontSize;

    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);
    
    print(gcf, [name, '.tif'],'-r600','-dtiff');
    
end

