function initialparams = getInitialParams(varargin)
  %cell2mat(varargin(1))
  %cell2mat(varargin(2))
  %cell2mat(varargin(3))
  if nargin == 2
    freqs = [0:1:length(cell2mat(varargin(2)))];
    amps = cell2mat(varargin(2));
  else
    freqs = cell2mat(varargin(2));
    amps = cell2mat(varargin(3));
  end
  x = cell2mat(varargin(1));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %B0 = mean([amps(1:(length(amps)/10)); ...
%	     amps(length(amps)-length(amps)/10+1:length(amps))]);
  B0 = mean([amps(1:10); ...
	     amps(length(amps)-10:length(amps))]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  B1 = (amps(length(amps))-amps(1))/(max(freqs)-min(freqs));
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Gam and D
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%%% Peak parameters, number determined by length of x.
  
  % Find the value in the frequency array that is nearest to the
  %  frequency the user has chosen in peaks array.
  for i=1:length(x)
    position(i) = getPosition(freqs,x(i));
  end
  
  % For and while loops to find the actual position in the frequency and
  %  amplitude arrays of the selected peak's maximum amplitude.
  for j=1:length(position)
    i=position(j);
    if amps(i)<amps(i-1)
      while amps(i)<amps(i-1)
	i=i-1;
      end
    else
      while amps(i)<amps(i+1)
	i=i+1;
      end
    end
    cp(j) = i; %cp -> center of peak, highest point
    %fprintf('cp is: %i\n',cp(j));
  end
  
  
  % The initial full width at half max of each peak is found by taking
  %  the(x,y) position of the maximum value of the amplitude of the peak,
  %  and calculating the point (x,y/2), or in this case (f,a/2) where f
  %  is the frequency of the maximum value of the peak, and a/2 is half
  %  the maximum amplitude.  
  % A search is now done to the left of f to find the point in the data
  %  set which is closest in geometric distance to (f,a/2), calculated as
  %  distance^2 = (f-f(i))^2 + (a/2-a(i))^2.
  % Note that calculating the sqrt(distance) isn't necessary.
  % The search is done to the right also, and the result is two points,
  %  lstop and rstop, where rstop > lstop, and width = rstop-lstop.
  
  % The initial skew it then calculated as:
  %  skew = log((f(cp)-f(lstop))/(f(rstop)-f(cp)))
  % Skew is then scaled by 1/1000 for reasons unknown to me to get a 
  %  good initial value.
  fmax = max(freqs);
  for j=1:length(cp)
    amax = amps(cp(j));
    amps = fmax*amps ./amax;  
    ld = (amps(cp(j))-amps(cp(j))/2);
    rd = ld;                                     
    for i=cp(j):-1:1
      ldist = sqrt((freqs(cp(j))-freqs(i))^2+...
		   (amps(cp(j))/2-amps(i))^2);
      if ldist <= ld 
    	ld = ldist;
    	lstop(j)=i;
      end
    end
    
    for i=cp(j):length(freqs)
      rdist = sqrt((freqs(cp(j))-freqs(i))^2+...
		   (amps(cp(j))/2-amps(i))^2);
      if rdist <= rd
    	rd = rdist;
    	rstop(j)=i;
      end
    end
 
    Gam(j)=freqs(rstop(j))-freqs(lstop(j));
    
    
    if freqs(rstop(j))-freqs(cp(j)) == 0 
      if freqs(cp(j))-freqs(lstop(j)) == 0
	D(j) = 0.0001;
      else
	%D(j) = log(freqs(cp(j))-freqs(lstop(j)))/10000;
	D(j) = 0.0001;
      end
    else
      if freqs(cp(j))-freqs(lstop(j)) == 0
	D(j) = 0.0001;
      else
	%D(j) = log((freqs(cp(j))-freqs(lstop(j)))/...
	%    (freqs(rstop(j))-freqs(cp(j))))/10;
	D(j) = 0.0001;
      end
    end
  end
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % C
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This gives a good approximate value for amplitude parameters:
  %   2^-(n-1) = Gam
  %   Abar = 4^n
  %   A = amax/Abar     where amax is maximum amplitude of each peak.
  
  for i=1:length(position)
    n(i) = -log2(Gam(i)) + 1;
    Abar(i) = 4^n(i);
    A(i)=amax/Abar(i);
    %fprintf('A is %f\n',A);
    %A(i)=.1;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Create initialparams list
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initialparams(1) = B0;
  initialparams(2) = B1;
  for i=1:length(position)
    initialparams(4*(i-1) + 3) = A(i);
    initialparams(4*(i-1) + 4) = D(i);
    initialparams(4*(i-1) + 5) = Gam(i);
    initialparams(4*(i-1) + 6) = freqs(cp(1,i));
  end
  %initialparams
  
  % Cleanup
  clear freqs amps B0 B1 position i j cp Gam D n Abar A
