
	integer*4 max_pnts			!max # of E' points.
	parameter  (max_pnts = 150)

	integer*4 max_sets			!max # of data sets.
	parameter (max_sets = 8)

	real*8 param1(40), param2(40), current_slope, current_offset
	real*8 innt, innp, fermip(197), epsn(197), aa, bb, cc, dd, ee, ff
	real*8 param3(40), param4(40)

	common /fits/ param1, param2, current_slope, current_offset,
     &   innp, innt, fermip, epsn, aa, bb, cc, dd, ee, ff, param3,
     &   param4
