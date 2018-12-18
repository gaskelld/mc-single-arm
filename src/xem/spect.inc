C ========================== BEGIN OF INCLUDE FILE =============================

	real*8 spect_sign
        parameter (spect_sign = +1.0)   !theta_effective=theta_spect+theta_tar
					!			    ^
	real*8 can_radius, t_can_side	!thickness (cm)
	real*8 t_can_front,t_can_back	!thickness (cm)
	real*8 t_exitwin		!thickness (cm)
	real*8 t_air			!thickness (cm)
	real*8 t_kevlar			!thickness (cm)
	real*8 t_mylar			!thickness (cm)

c	parameter (t_can_front=.0028*2.54)	!2.8 mil crytarget entrance window.
	parameter (t_can_front=.005*2.54)	!5 mil crytarget entrance window. (tuna can)
	parameter (t_can_back=.005*2.54)	!5 mil crytarget exit window.
***CHECK THIS!****
	parameter (t_can_side=.005*2.54)	!5 mil? crytarget walls.
	parameter (can_radius=1.325*2.54)       !radius of cell.
	parameter (t_exitwin=.016*2.54)		!16 mil chamber exit window.
        parameter (t_air=18+40.17)		!~18 cm air gap (rolf:6-8").
        parameter (t_kevlar=.015*2.54)		!15 mil kevlar entrance window
        parameter (t_mylar=.005*2.54)		!5 mil mylar entrance window

C ========================== END OF INCLUDE FILE ===============================
