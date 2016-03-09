function myPlot(params,spectrum)
  
  delete(findobj(gca,'color','red'));
  
  if spectrum.col == 1
    data = bw(params,[0:1:length(spectrum.one)]);
    plot([0:1:length(spectrum.one)],spectrum.one,'r-');
  else
    data = bw(params,spectrum.one);
    plot(spectrum.one, data,'r-');
  end
  