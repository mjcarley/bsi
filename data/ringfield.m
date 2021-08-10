nr = 65 ;
nz = 33 ;

dz = 1 ;

[r,z] = meshgrid(linspace(0,1.5,nr), linspace(-dz, dz, nz)) ;

p = vfield(1, 0.3, r, z) ;

dz = diff(z(1:2)) ;
dr = diff(r'(1:2)) ;

[pr,pz]=gradient(p, dr, dz) ;
uz = pr./r ;
ur = -pz./r ;
