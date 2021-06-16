function mkring(file, s, tol)

  ng = 17 ;
  dx = 1.5 ;
  r0 = 1 ;
  
  [x,y,z] = meshgrid(linspace(-dx, dx, ng)) ;

  x = x(:) ; y = y(:) ; z = z(:) ;

  r = sqrt(x.^2 + y.^2) ;
  th = atan2(y, x) ;
  R = sqrt((r-r0).^2 + z.^2) ;
  w0 = exp(-R.^2/s^2) ;

  w = [-w0.*sin(th) w0.*cos(th) 0*w0] ;

  fid = fopen(file, "w") ;

  fprintf(fid, "%d", length(x)) ;
  fprintf(fid, " %d\n", 3) ;

  dat = [x y z w]' ;

  fprintf(fid, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", dat) ;
  
  fclose(fid) ;
  
