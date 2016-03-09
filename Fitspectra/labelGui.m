function labelGui(Action)
  
  
  if nargin == 0
    
    figpos = get(gcf,'Position');
    bXSize = 100;
    bYSize = 20;
    miniXsize = 20;
    miniYsize = 20;
    pos1 = get(findobj('tag','axes1'),'Position');
    xpos1 = figpos(3)*( (pos1(1)+pos1(3)) + ...
			(1-(pos1(1)+pos1(3)))/2 )-bXSize/2;
    ypos1 = figpos(4)*( (pos1(2)+pos1(4)) ) - bYSize;
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UI Menu for setting axes label text 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AO = uimenu('Label','Axes Options');
    uimenu(AO,'Label','Set axes1 X-label', ...
	   'Callback','labelGui(''setXlabel1'')');
    uimenu(AO,'Label','Set axes1 Y-label', ...
	   'Callback','labelGui(''setYlabel1'')');
    uimenu(AO,'Label','Set axes1 Title', ...
	   'Callback','labelGui(''setTitle1'')');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start label/title uicontrols
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%  
    % X-label uicontrols 1
    %%%%%%%%%%%%%%%%%%%%%% 
    uicontrol('Style', 'checkbox', ...
	      'Position', [xpos1 ypos1 bXSize bYSize], ...
	      'String', 'X-label: off', ...
	      'Tag', 'xButt1', ...
	      'Value', 0, ...
	      'Callback', 'labelGui(''hideX1'')');
    uicontrol('Style', 'edit', ...
	      'Position', [xpos1 ypos1 bXSize bYSize], ...
	      'String', '', ...
	      'Tag', 'editXlabel1', ...
	      'Visible', 'off', ...
	      'Callback', 'labelGui(''changeXlabel1'')');
   
    %%%%%%%%%%%%%%%%%%%%%%   
    % Y-label uicontrols 1
    %%%%%%%%%%%%%%%%%%%%%%
    uicontrol('Style', 'checkbox', ...
	      'Position', [xpos1 ypos1-1.5*bYSize bXSize bYSize], ...
	      'String', 'Y-label: off', ...
	      'Tag', 'yButt1', ...
	      'Value', 0, ...
	      'Callback', 'labelGui(''hideY1'')');
    uicontrol('Style', 'edit', ...
	      'Position', [xpos1 ypos1-1.5*bYSize bXSize bYSize], ...
	      'String', '', ...
	      'Tag', 'editYlabel1', ...
	      'Visible', 'off', ...
	      'Callback', 'labelGui(''changeYlabel1'')');
    
    %%%%%%%%%%%%%%%%%%%%   
    % Title uicontrols 1
    %%%%%%%%%%%%%%%%%%%%
    uicontrol('Style', 'checkbox', ...
	      'Position', [xpos1 ypos1-3.0*bYSize bXSize bYSize], ...
	      'String', 'Title: off', ...
	      'Tag', 'titleButt1', ...
	      'Value', 0, ...
	      'Callback', 'labelGui(''hideTitle1'')');
    uicontrol('Style', 'edit', ...
	      'Position', [xpos1 ypos1-3.0*bYSize bXSize bYSize], ...
	      'String', '', ...
	      'Tag', 'editTitle1', ...
	      'Visible', 'off', ...
	      'Callback', 'labelGui(''changeTitle1'')');
    
    % Do for CWI only
    if strcmp(get(gcf,'tag'),'cwiFig')
      pos2 = get(findobj('tag','axes2'),'Position');
      xpos2 = figpos(3)*( (pos2(1)+pos2(3)) + ...
			  (1-(pos1(1)+pos2(3)))/2 )-bXSize/2;
      ypos2 = figpos(4)*( (pos2(2)+pos2(4)) ) - bYSize;
      
      uimenu(AO,'Label','Set axes2 X-label', ...
	     'Callback','labelGui(''setXlabel2'')','Separator','on');
      uimenu(AO,'Label','Set axes2 Y-label', ...
	     'Callback','labelGui(''setYlabel2'')');
      uimenu(AO,'Label','Set axes2 Title', ...
	     'Callback','labelGui(''setTitle2'')');
      
      %%%%%%%%%%%%%%%%%%%%%% 
      % X-label uicontrols 2
      %%%%%%%%%%%%%%%%%%%%%% 
      uicontrol('Style', 'checkbox', ...
		'Position', [xpos2 ypos2 bXSize bYSize], ...
		'String', 'X-label: off', ...
		'Tag', 'xButt2', ...
		'Value', 0, ...
		'Callback', 'labelGui(''hideX2'')');
      uicontrol('Style', 'edit', ...
		'Position', [xpos2 ypos2 bXSize bYSize], ...
		'String', '', ...
		'Tag', 'editXlabel2', ...
		'Visible', 'off', ...
		'Callback', 'labelGui(''changeXlabel2'')');
      
      %%%%%%%%%%%%%%%%%%%%%%
      % Y-label uicontrols 2
      %%%%%%%%%%%%%%%%%%%%%%
      uicontrol('Style', 'checkbox', ...
		'Position', [xpos2 ypos2-1.5*bYSize bXSize bYSize], ...
		'String', 'Y-label: off', ...
		'Tag', 'yButt2', ...
		'Value', 0, ...
		'Callback', 'labelGui(''hideY2'')');
      uicontrol('Style', 'edit', ...
		'Position', [xpos2 ypos2-1.5*bYSize bXSize bYSize], ...
		'String', '', ...
		'Tag', 'editYlabel2', ...
		'Visible', 'off', ...
		'Callback', 'labelGui(''changeYlabel2'')');
      
      %%%%%%%%%%%%%%%%%%%%
      % Title uicontrols 2
      %%%%%%%%%%%%%%%%%%%%
      uicontrol('Style', 'checkbox', ...
		'Position', [xpos2 ypos2-3.0*bYSize bXSize bYSize], ...
		'String', 'Title: off', ...
		'Tag', 'titleButt2', ...
		'Value', 0, ...
		'Callback', 'labelGui(''hideTitle2'')');
      uicontrol('Style', 'edit', ...
		'Position', [xpos2 ypos2-3.0*bYSize bXSize bYSize], ...
		'String', '', ...
		'Tag', 'editTitle2', ...
		'Visible', 'off', ...
		'Callback', 'labelGui(''changeTitle2'')');  
    end
    
    
  else
    
    switch Action
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start axes/title Actions
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%
      % X-axes1
      %%%%%%%%%
     case 'hideX1'       % 
      h=findobj('tag','axes1');axes(h);
      xval = get(findobj('Tag', 'xButt1'),'Value');
      if xval == 0 
	xlabel('');
	set(findobj('Tag', 'xButt1'),'String','X-label: off');
      elseif xval == 1 & strcmp(get(findobj('Tag','editXlabel1'),'String'),'')
	xlabel('Frequency (KHz)');
	set(findobj('Tag', 'xButt1'),'String','X-label: on');
      else
	xlabel(get(findobj('Tag','editXlabel1'),'String'));
	set(findobj('Tag', 'xButt1'),'String','X-label: on');
      end
      
     case 'setXlabel1'
      if strcmp(get(findobj('tag','editXlabel1'),'string'),'') == 0 
	set(findobj('Tag','xButt1'),'visible','off');
	set(findobj('Tag','editXlabel1'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editXlabel1'),'String');
	set(findobj('Tag','editXlabel1'),'String',oldLabel);
	set(findobj('Tag','editXlabel1'),'Visible','on');
      end
      set(findobj('Tag','xButt1'),'value',1);
      
     case 'changeXlabel1'
      h=findobj('tag','axes1');axes(h);
      newXlabel = get(findobj('Tag','editXlabel1'),'String');
      xlabel(newXlabel);
      set(findobj('Tag','editXlabel1'),'visible','off');
      set(findobj('Tag','xButt1'),'String', 'X-label: on');
      set(findobj('Tag','xButt1'),'visible','on');
      
      %%%%%%%%%
      % X-axes2
      %%%%%%%%%
     case 'hideX2'       % 
      h=findobj('tag','axes2');axes(h);
      xval = get(findobj('Tag', 'xButt2'),'Value');
      if xval == 0 
	xlabel('');
	set(findobj('Tag', 'xButt2'),'String','X-label: off');
      elseif xval == 1 & strcmp(get(findobj('Tag','editXlabel2'),'String'),'')
	xlabel('time (seconds)');
	set(findobj('Tag', 'xButt2'),'String','X-label: on');
      else
	xlabel(get(findobj('Tag','editXlabel2'),'String'));
	set(findobj('Tag', 'xButt2'),'String','X-label: on');
      end
      
     case 'setXlabel2'
      if strcmp(get(findobj('tag','editXlabel2'),'string'),'') == 0 
	set(findobj('Tag','xButt2'),'visible','off');
	set(findobj('Tag','editXlabel2'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editXlabel2'),'String');
	set(findobj('Tag','editXlabel2'),'String',oldLabel);
	set(findobj('Tag','editXlabel2'),'Visible','on');
      end
      
     case 'changeXlabel2'
      h=findobj('tag','axes2');axes(h);
      newXlabel = get(findobj('Tag','editXlabel2'),'String');
      xlabel(newXlabel);
      set(findobj('Tag','editXlabel2'),'visible','off');
      set(findobj('Tag','xButt2'),'String', 'X-label: on');
      set(findobj('Tag','xButt2'),'visible','on');
      
      %%%%%%%%%
      % Y-axes1
      %%%%%%%%%
     case 'hideY1'
      h=findobj('tag','axes1');axes(h);
      yval = get(findobj('Tag', 'yButt1'),'Value');
      if yval == 0 
	ylabel('');
       set(findobj('Tag', 'yButt1'),'String','Y-label: off');
      elseif yval == 1 & strcmp(get(findobj('Tag','editYlabel1'),'String'),'')
	ylabel('amplitude');
	set(findobj('Tag', 'yButt1'),'String','Y-label: on');
      else
	ylabel(get(findobj('Tag','editYlabel1'),'String'));
	set(findobj('Tag', 'yButt1'),'String','Y-label: on');
      end
      
     case 'setYlabel1'
      if strcmp(get(findobj('tag','editYlabel1'),'string'),'') == 0 
	set(findobj('Tag','yButt1'),'visible','off');
	set(findobj('Tag','editYlabel1'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editYlabel1'),'String');
	set(findobj('Tag','editYlabel1'),'String',oldLabel);
	set(findobj('Tag','editYlabel1'),'Visible','on');
      end
      set(findobj('Tag','yButt1'),'value',1);
 
     case 'changeYlabel1'
      h=findobj('tag','axes1');axes(h);
      newYlabel = get(findobj('Tag','editYlabel1'),'String');
      ylabel(newYlabel);
      set(findobj('Tag','editYlabel1'),'visible','off');
      set(findobj('Tag','yButt1'),'String', 'Y-label: on');
      set(findobj('Tag','yButt1'),'visible','on');

      
      %%%%%%%%%
      % Y-axes2
      %%%%%%%%%
     case 'hideY2'
      h=findobj('tag','axes2');axes(h);
      yval = get(findobj('Tag', 'yButt2'),'Value');
      if yval == 0 
	ylabel('');
       set(findobj('Tag', 'yButt2'),'String','Y-label: off');
      elseif yval == 1 & strcmp(get(findobj('Tag','editYlabel2'),'String'),'')
	ylabel('dV/V');
	set(findobj('Tag', 'yButt2'),'String','Y-label: on');
      else
	ylabel(get(findobj('Tag','editYlabel2'),'String'));
	set(findobj('Tag', 'yButt2'),'String','Y-label: on');
      end
      
     case 'setYlabel2'
      if strcmp(get(findobj('tag','editYlabel2'),'string'),'') == 0 
	set(findobj('Tag','yButt2'),'visible','off');
	set(findobj('Tag','editYlabel2'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editYlabel2'),'String');
	set(findobj('Tag','editYlabel2'),'String',oldLabel);
	set(findobj('Tag','editYlabel2'),'Visible','on');
      end
      set(findobj('Tag','yButt2'),'value',1);
      
     case 'changeYlabel2'
      h=findobj('tag','axes2');axes(h);
      newYlabel = get(findobj('Tag','editYlabel2'),'String');
      ylabel(newYlabel);
      set(findobj('Tag','editYlabel2'),'visible','off');
      set(findobj('Tag','yButt2'),'String', 'Y-label: on');
      set(findobj('Tag','yButt2'),'visible','on');
      
      %%%%%%%%%
      % Title1
      %%%%%%%%%
     case 'hideTitle1'       % 
      h=findobj('tag','axes1');axes(h);
      tval = get(findobj('Tag', 'titleButt1'),'Value');
      if tval == 0 
	title('');
	set(findobj('Tag', 'titleButt1'),'String','Title: off');
      elseif tval == 1 & strcmp(get(findobj('Tag','editTitle1'),'String'),'')
	title('Amplitude versus Frequency Plot');
	set(findobj('Tag', 'titleButt1'),'String','Title: on');
      else
	title(get(findobj('Tag','editTitle1'),'String'));
	set(findobj('Tag', 'titleButt1'),'String','Title: on');
      end
      
     case 'setTitle1'
      if strcmp(get(findobj('tag','editTitle1'),'string'),'') == 0 
	set(findobj('Tag','titleButt1'),'visible','off');
	set(findobj('Tag','editTitle1'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editTitle1'),'String');
	set(findobj('Tag','editTitle1'),'String',oldLabel);
	set(findobj('Tag','editTitle1'),'Visible','on');
      end
      set(findobj('Tag','titleButt1'),'value',1);

     case 'changeTitle1'
      h=findobj('tag','axes1');axes(h);
      newTitle = get(findobj('Tag','editTitle1'),'String');
      title(newTitle);
      set(findobj('Tag','editTitle1'),'visible','off');
      set(findobj('Tag','titleButt1'),'String', 'Title: on');
      set(findobj('Tag','titleButt1'),'visible','on');
      
      
      %%%%%%%%%
      % Title2
      %%%%%%%%%
     case 'hideTitle2'       % 
      h=findobj('tag','axes2');axes(h);
      tval = get(findobj('Tag', 'titleButt2'),'Value');
      if tval == 0 
	title('');
	set(findobj('Tag', 'titleButt2'),'String','Title: off');
      elseif tval == 1 & strcmp(get(findobj('Tag','editTitle2'),'String'),'')
	title('Spectrum Plot');
	set(findobj('Tag', 'titleButt2'),'String','Title: on');
      else
	title(get(findobj('Tag','editTitle2'),'String'));
	set(findobj('Tag', 'titleButt2'),'String','Title: on');
      end
      
     case 'setTitle2'
      if strcmp(get(findobj('tag','editTitle2'),'string'),'') == 0 
	set(findobj('Tag','titleButt2'),'visible','off');
	set(findobj('Tag','editTitle2'),'visible','on');
      else
	oldLabel = get(findobj('Tag','editTitle2'),'String');
	set(findobj('Tag','editTitle2'),'String',oldLabel);
	set(findobj('Tag','editTitle2'),'Visible','on');
      end
      set(findobj('Tag','titleButt2'),'value',1);
      
     case 'changeTitle2'
      h=findobj('tag','axes2');axes(h);
      newTitle = get(findobj('Tag','editTitle2'),'String');
      title(newTitle);
      set(findobj('Tag','editTitle2'),'visible','off');
      set(findobj('Tag','titleButt2'),'String', 'Title: on');
      set(findobj('Tag','titleButt2'),'visible','on');      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % end axes/title Actions
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     end
  end
