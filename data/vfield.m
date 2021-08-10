function [p,v]=vfield(r0, s, r, z)

  p = v = zeros(size(r)) ;

  rmax = 4 ;
  ng = 192 ;
  nt = 96 ;

  [t, w] = grule(ng) ;
  
  for i=1:length(r(:))
    if ( r(i) == 0 ) tmax = pi/2 ; else tmax = 3*pi/2 ; endif

    th = linspace(-pi/2, tmax, nt) ;
    for j=1:length(th)
      C = cos(th(j)) ;
      if ( C < 0 ) rr = -r(i)/C ; else rr = rmax ; endif

      rho = 0.5*rr*(1+t) ;
      r1 = r(i) + rho*C ;
      z1 = z(i) + rho*sin(th(j)) ;

      R1 = sqrt((z(i)-z1).^2 + (r(i)-r1).^2) ;
      R2 = sqrt((z(i)-z1).^2 + (r(i)+r1).^2) ;
      lm = (R2 - R1)./(R2 + R1) ;

      lm(find(isnan(lm))) = 0 ;
      
      om = exp(-((r1-r0).^2 + z1.^2)/s^2) ;

      [K,E] = ellipke(lm.^2) ;

      f = (R2 + R1).*(K - E).*om/2/pi ;

      p(i) += sum(f.*w.*rho)*0.5*rr*diff(th(1:2)) ;
      
    endfor

  endfor
