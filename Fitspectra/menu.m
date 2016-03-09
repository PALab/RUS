function menu(Action);
  
  
  if nargin == 0
    
  else
    
    switch Action
      
     case 'toCWI'
      if strcmp(get(gcf,'tag'),'cwiFig')
      else
       close(gcbf);
       cwi;    
      end
      
     case 'toFIT'
      if strcmp(get(gcf,'tag'),'fitspectraFig')
      else
       close(gcbf);
       fitspectra;    
      end
      
     case 'toXYLine'
      if strcmp(get(gcf,'tag'),'XYLineFig')
      else
	close(gcbf);
	xyLineScan;
      end
      
     case 'toXYPlane'
      if strcmp(get(gcf,'tag'),'XYPlaneFig')
      else
	close(gcbf);
	xyPlaneScan;
      end
      
     case 'toCircle'
      if strcmp(get(gcf,'tag'),'circleFig')
      else
	close(gcbf);
	CircleScan;
      end
      
    end
  end