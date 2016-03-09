% fit a Breit-Wigner model to a measured amplitude spectrum
function newparams = fit(params,mag,freq)
  
  NP = (length(params)-2)/4;
  l = 1:NP;                   
  peaks = params(4*(l-1)+6);
  %nf=length(freq);          
    
  %test bounds
  LB = [-inf -inf];
  UB = [inf inf];
  for i=1:length(l)
    LB = [LB -inf -inf -inf (peaks(i)-.01*peaks(i))];
    UB = [UB inf inf inf (peaks(i)+.01*peaks(i))];
  end
  
   [newparams,resnorm, residual, exitflag, output,options] = ...
      lsqcurvefit('bw',params,freq',mag', LB, UB);
  %fprintf(1,'The residual is: %f\n',residual);
  
  
 
  






