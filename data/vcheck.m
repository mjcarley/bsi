load velocity.dat

x = velocity(:, 1:3) ;
uc = velocity(:, 4:6) ;

load velocity0.dat
u0 = velocity0(:, 4:6) ;


##ii = find(x(:,2) == 0) ;
