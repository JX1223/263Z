use <./parametric-loft/loft.scad>;
use <./parametric-loft/pcurve.scad>;
include <./airfoil.scad>;

  
loft([
	[],
    for (i=[0:1:len(airfoil)-1]) indexed_parametrize(airfoil[i])
    ,
	[],
]);
