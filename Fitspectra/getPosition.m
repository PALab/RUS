function position = getPosition(array,value)
  
  diff_value = abs(array-value);
  
  
  position = 1;
  for i=2:length(array)
    if diff_value(i)<diff_value(position)
      position = i;
    end
  end
 
  
  clear diff_value