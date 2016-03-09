function setParamSliders(params, val)
  
%%%
% Put parameters of peak 1 on sliders and slider text, set ranges
%%%
  
% B0
  set(findobj('Tag','B0Size'), 'String',params(1));
  set(findobj('Tag','B0Slider'),'value',params(1));
  set(findobj('Tag','B0Slider'),'min',params(1)-2*abs(params(1)));
  set(findobj('Tag','B0Slider'),'max',params(1)+2*abs(params(1)));
  set(findobj('Tag','minB0'),'String',params(1)-2*abs(params(1)));
  set(findobj('Tag','maxB0'),'String',params(1)+2*abs(params(1)));
  
  % B1
  set(findobj('Tag','B1Size'), 'String',params(2));
  set(findobj('Tag','B1Slider'),'value',params(2));
  set(findobj('Tag','B1Slider'),'min',params(2)-2*abs(params(2)));
  set(findobj('Tag','B1Slider'),'max',params(2)+2*abs(params(2)));
  set(findobj('Tag','minB1'),'String',params(2)-2*abs(params(2)));
  set(findobj('Tag','maxB1'),'String',params(2)+2*abs(params(2)));
  
  
  % C
  set(findobj('Tag','CSize'),'String',params(4*(val-1)+3));
  set(findobj('Tag','CSlider'), 'value',params(4*(val-1)+3));
  set(findobj('Tag','CSlider'),'min',params(4*(val-1)+3)-2*abs(params(4*(val-1)+3)));
  set(findobj('Tag','CSlider'),'max',params(4*(val-1)+3)+0.5*abs(params(4*(val-1)+3)));
  set(findobj('Tag','minC'),'String',params(4*(val-1)+3)-2*abs(params(4*(val-1)+3)));
  set(findobj('Tag','maxC'),'String',params(4*(val-1)+3)+0.5*abs(params(4*(val-1)+3)));
  
  % D
  set(findobj('Tag','DSize'), 'String',params(4*(val-1)+4));
  set(findobj('Tag','DSlider'),'value',params(4*(val-1)+4));
  set(findobj('Tag','DSlider'),'min',params(4*(val-1)+4)-abs(params(4*(val-1)+4)));
  set(findobj('Tag','DSlider'),'max',params(4*(val-1)+4)+abs(params(4*(val-1)+4)));
  set(findobj('Tag','minD'),'String',params(4*(val-1)+4)-abs(params(4*(val-1)+4)));
  set(findobj('Tag','maxD'),'String',params(4*(val-1)+4)+abs(params(4*(val-1)+4)));
  
  % Gam
  set(findobj('Tag','GamSize'), 'String',params(4*(val-1)+5));
  set(findobj('Tag','GamSlider'),'value',params(4*(val-1)+5));
  set(findobj('Tag','GamSlider'),'min',0);
  set(findobj('Tag','GamSlider'),'max',params(4*(val-1)+5)+0.5*abs(params(4*(val-1)+5)));
  set(findobj('Tag','minGam'),'String','0');
  set(findobj('Tag','maxGam'),'String',params(4*(val-1)+5)+0.5*abs(params(4*(val-1)+5)));
  
  % Freq
  set(findobj('Tag','FreqSize'), 'String',params(4*(val-1)+6));
  set(findobj('Tag','FreqSlider'),'value',params(4*(val-1)+6));
  set(findobj('Tag','FreqSlider'),'min',params(4*(val-1)+6)-abs(params(4*(val-1)+6))/100);
  set(findobj('Tag','FreqSlider'),'max',params(4*(val-1)+6)+abs(params(4*(val-1)+6))/100);
  set(findobj('Tag','minFreq'),'String',params(4*(val-1)+6)-abs(params(4*(val-1)+6))/100);
  set(findobj('Tag','maxFreq'),'String',params(4*(val-1)+6)+abs(params(4*(val-1)+6))/100);
  
  
