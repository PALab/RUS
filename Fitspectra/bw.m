function A = bw(x,f);
  % given an array of frequencies and BW parameters, this function returns
  % and array of amplitude values
    
  N = (length(x)-2)/4;
  M = length(f);
  l = 1:N;

  % these are the parameters of a Breit-Wigner representation of the spectrum
  % cf Bertlesen p. 70 for notation.  The parameters are all unpacked
  % from a single vector x.
  B0 = x(1);
  B1 = x(2);
  C = x(4*(l-1)+3);
  D = x(4*(l-1)+4);
  Gam = x(4*(l-1)+5);
  peak = x(4*(l-1)+6);
  
  fmax = max(peak);
  fmin = min(peak);
  f0 = (fmax + fmin)/2;
  
  % A is the amplitude as a function of frequency
  A(1:M) = 0;
  for l = 1:M % loop over number of input frequencies
  for k = 1:N % loop over assumed number of peaks
  A(l) = A(l)+(C(k) + D(k)*(f(l) - f0))/((f(l) - peak(k))^2 + .25 * Gam(k)^2);
  end
  A(l) = A(l) + B0+ B1 * (f(l) - f0);
  end 
  
% references:
% G. Breit, E. Wigner, Capture of Slow Neutrons, Phys. Rev. 49 (1936) 519. 
% G. Breit, Theory of Resonance Reactions, in: Handbuch der Physik
% XLI/1 Springer, Berlin, Heidelberg, 1959. 
