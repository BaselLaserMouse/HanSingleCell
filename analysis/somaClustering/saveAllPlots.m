 [n,c,abrv]=brainAreaNames.visualAreas;

 for ii=2:13, 
    [~,data{ii}]=clusterPos(cleanCells,ii);
    drawnow
    set(gcf,'PaperPosition',[0,0,25,25],'units','inches')

    print('-dpng', fullfile('images',abrv{ii}))
    print('-depsc', fullfile('images',abrv{ii}))

end
