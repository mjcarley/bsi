load nodes.dat

load nodes1.dat

idx = [] ;

for i=1:size(nodes,1)
  ii = find(nodes1(:,1)==nodes(i,1) &
	    nodes1(:,2)==nodes(i,2) &
	    nodes1(:,3)==nodes(i,3) ) ;
  if ( any(abs(nodes1(ii,4:6) - nodes(i, 4:6)) > 1e-3) )
    break ;
  endif
  
#  if ( length(ii) ~= 1 )
#    error("") ;
#  endif
#  idx = [idx; ii] ;
endfor

if 0

[xn, i] = sort(nodes(:,1)) ;
nodes = nodes(i,:) ;
[xn, i] = sort(nodes(:,2)) ;
nodes = nodes(i,:) ;
[xn, i] = sort(nodes(:,3)) ;
nodes = nodes(i,:) ;

[xn, i] = sort(nodes1(:,1)) ;
nodes1 = nodes1(i,:) ;
[xn, i] = sort(nodes1(:,2)) ;
nodes1 = nodes1(i,:) ;
[xn, i] = sort(nodes1(:,3)) ;
nodes1 = nodes1(i,:) ;
endif
