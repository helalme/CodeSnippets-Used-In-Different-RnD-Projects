C This UMAT subroutine is a modified version of Huang's original UMAT. Necessary modification is done for
C all related functions. Flow rule, slip systems, hardening laws, hardening stresses, and related all variables and functions are modified accordingly.
C 14 slip system is applicable for Al-rich TiAl single crystals at high temperature  

      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C	  Function for setting displacement rate
      INCLUDE 'ABA_PARAM.INC'
C      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3),TIME(2),COORDS(3)
      U(1)=-0.01*TIME(1)
	  U(2)=0.0
	  U(3)=0.0
      RETURN
      END 

C---------------------------------------------------------------
      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
C_HC_ necessary to put total number of slip systems
      PARAMETER (ND=14)
C_HC_ necessary to put total number of slip systems
	  
C-----  The parameter ND determines the dimensions of the arrays in 
C     this subroutine.  The current choice 150 is a upper bound for a 
C     cubic crystal with up to three sets of slip systems activated.  
C     Users may reduce the parameter ND to any number as long as larger
C     than or equal to the total number of slip systems in all sets.  
C     For example, if {110}<111> is the only set of slip system 
C     potentially activated, ND could be taken as twelve (12).  
c
      include 'aba_param.inc'
c
      CHARACTER*8 CMNAME
      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), 
     2          SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND), 
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND), 
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3), 
     5          FSLIP(ND), DFDZSP(ND), DDEMSD(6,ND), 
     6          DDGDDE(ND,6), 
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND), 
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3)
	 
!      DIMENSION FSLIP1(ND), STRES1(6), GAMMA1(ND), TAUSP1(ND), 
!     2          GSLP1(ND), SPNOR1(3,ND), SPDIR1(3,ND), DDSDE1(6,6),
!     3          DSOLD(6), DGAMOD(ND), DTAUOD(ND), DGSPOD(ND), 
!     4          DSPNRO(3,ND), DSPDRO(3,ND), 
!     5          DHDGDG(ND,ND)	 
	 
      DIMENSION DelX(ND),DelR(ND),H(ND,ND)
	  
C-----  DDEMSD -- double dot product of the elastic moduli tensor with 
C                 the slip deformation tensor plus, only for finite 
C                 rotation, the dot product of slip spin tensor with 
C                 the stress
C-----  DDGDDE -- derivatice of the shear strain increments in slip 
C                 systems w.r.t. the increment of strains
C-----  DELATS -- strain-increments associated with lattice stretching
C                 DELATS(1) - DELATS(3) -- normal strain increments
C                 DELATS(4) - DELATS(6) -- engineering shear strain 
C                                          increments
C-----  DSPIN  -- spin-increments associated with the material element
C                 DSPIN(1) -- component 12 of the spin tensor
C                 DSPIN(2) -- component 31 of the spin tensor
C                 DSPIN(3) -- component 23 of the spin tensor
C-----  DVGRAD -- increments of deformation gradient in the current 
C                 state, i.e. velocity gradient times the increment of 
C                 time
C       SDV(NSTATV)= total nos of slip systems in all set
C       SDV(NSTATV-1)= total nos of slip systems in the 3rd set
C       SDV(NSTATV-2)= total nos of slip systems in the 2nd set
C       SDV(NSTATV-3)= total nos of slip systems in the 1st set
C       SDV(NSTATV-4)= Sum(abs(initial shear strain of SS) +
C                      Sum(abs(Dgamma(t1))) + Sum(abs(Dgamma(t2)))+
C                      Sum(abs(Dgamma(t3))) + ........ till end of time step
C       SDV(10*NSLPTL+1)~ SDV(11*NSLPTL)= for kinematic hardening variable X
	  
C-----  Elastic matrix in local cubic crystal system: DLOCAL
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.
         END DO
      END DO

      CHECK=0.
      DO J=10,21
         CHECK=CHECK+ABS(PROPS(J))
      END DO

      IF (CHECK.EQ.0.) THEN
         DO J=4,9
            CHECK=CHECK+ABS(PROPS(J))
         END DO

         IF (CHECK.EQ.0.) THEN

            IF (PROPS(3).EQ.0.) THEN

C-----  Isotropic material
               GSHEAR=PROPS(1)/2./(1.+PROPS(2))
               E11=2.*GSHEAR*(1.-PROPS(2))/(1.-2.*PROPS(2))
               E12=2.*GSHEAR*PROPS(2)/(1.-2.*PROPS(2))

               DO J=1,3
                  DLOCAL(J,J)=E11

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=E12
                  END DO

                  DLOCAL(J+3,J+3)=GSHEAR
               END DO

            ELSE

C-----  Cubic material
               DO J=1,3
                  DLOCAL(J,J)=PROPS(1)

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=PROPS(2)
                  END DO

                  DLOCAL(J+3,J+3)=PROPS(3)
               END DO

            END IF

         ELSE

C-----  Orthotropic metarial
            DLOCAL(1,1)=PROPS(1)
            DLOCAL(1,2)=PROPS(2)
            DLOCAL(2,1)=PROPS(2)
            DLOCAL(2,2)=PROPS(3)

            DLOCAL(1,3)=PROPS(4)
            DLOCAL(3,1)=PROPS(4)
            DLOCAL(2,3)=PROPS(5)
            DLOCAL(3,2)=PROPS(5)
            DLOCAL(3,3)=PROPS(6)

            DLOCAL(4,4)=PROPS(7)
            DLOCAL(5,5)=PROPS(8)
            DLOCAL(6,6)=PROPS(9)

         END IF

      ELSE

C-----  General anisotropic material
         ID=0
         DO J=1,6
            DO I=1,J
               ID=ID+1
               DLOCAL(I,J)=PROPS(ID)
               DLOCAL(J,I)=DLOCAL(I,J)
            END DO
         END DO
      END IF

C-----  Rotation matrix: ROTATE, i.e. direction cosines of [100], [010]
C     and [001] of a cubic crystal in global system
C
      CALL ROTATION (PROPS(57), ROTATE)

C-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL 
C     to global elastic matrix D
C
      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

C-----  Elastic matrix in global system: D
C     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose
C
      DO J=1,6
         DO I=1,6
            D(I,J)=0.
         END DO
      END DO

      DO J=1,6
         DO I=1,J

            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO

            D(J,I)=D(I,J)

         END DO
      END DO

C-----  Total number of sets of slip systems: NSET
      NSET=NINT(PROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         STOP
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*) 
     2     '***ERROR - more than three sets of slip systems'
         STOP
      END IF

C-----  Implicit integration parameter: THETA
      THETA=PROPS(145)

C-----  Finite deformation ?
C-----  NLGEOM = 0,   small deformation theory
C       otherwise, theory of finite rotation and finite strain, Users 
C     must declare "NLGEOM" in the input file, at the *STEP card
C
      IF (PROPS(146).EQ.0.) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF

C-----  Iteration?
C-----  ITRATN = 0, no iteration
C       otherwise, iteration (solving increments of stresses and 
C     solution dependent state variables)
C
      IF (PROPS(153).EQ.0.) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF

      ITRMAX=NINT(PROPS(154))
      GAMERR=PROPS(155)

      NITRTN=-1
 
!	  DO I=1,NTENS
!         DSOLD(I)=0.
!      END DO
!
!     DO J=1,ND
!         DGAMOD(J)=0.
!         DTAUOD(J)=0.
!        DGSPOD(J)=0.
!         DO I=1,3
!            DSPNRO(I,J)=0.
!            DSPDRO(I,J)=0.
!         END DO
!      END DO
	  
	  
C-----  Increment of spin associated with the material element: DSPIN
C     (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO

            TERM(J,J)=TERM(J,J)+1.D0
            TRM0(J,J)=TRM0(J,J)-1.D0
         END DO

         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP)

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)

      END IF

C-----  Increment of dilatational strain: DEV
      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO

C-----  Parameter NITRTN: number of iterations
C       NITRTN = 0 --- no-iteration solution
C
      NITRTN=NITRTN+1
      
C-----  Check whether the current stress state is the initial state
      IF (STATEV(1).EQ.0.) THEN

C-----  Initial state
C
C-----  Generating the following parameters and variables at initial 
C     state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Unit vectors in initial slip directions SLPDIR
C          Unit normals to initial slip planes SLPNOR
C
         NSLPTL=0
         DO I=1,NSET
            ISPNOR(1)=NINT(PROPS(25+8*I))
            ISPNOR(2)=NINT(PROPS(26+8*I))
            ISPNOR(3)=NINT(PROPS(27+8*I))

            ISPDIR(1)=NINT(PROPS(28+8*I))
            ISPDIR(2)=NINT(PROPS(29+8*I))
            ISPDIR(3)=NINT(PROPS(30+8*I))

            CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), SLPDIR(1,NSLPTL+1), 
     2                    SLPNOR(1,NSLPTL+1), ROTATE)

            NSLPTL=NSLPTL+NSLIP(I)
         END DO

         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*) 
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
         END IF

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

C-----  Initial value of state variables: unit normal to a slip plane 
C     and unit vector in a slip direction
C
         STATEV(NSTATV)=FLOAT(NSLPTL)
         DO I=1,NSET
            STATEV(NSTATV-4+I)=FLOAT(NSLIP(I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=SLPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=SLPDIR(I,J)
            END DO
         END DO

C-----  Initial value of the current strength for all slip systems
C
!         CALL GSLPINIT (STATEV(1), NSLIP, NSLPTL, NSET, PROPS(97))
         
C-----  Initial value of shear strain in slip systems
         DO I=1,NSLPTL
            STATEV(NSLPTL+I)=0.
            STATEV(9*NSLPTL+I)=0.
         END DO
CFIX--  Initial value of cumulative shear strain in each slip systems         
         STATEV(NSTATV-4)=0.
         
C-----  Initial value of the isotropic hardening stress r for all slip systems
         CALL rSlipInit (STATEV(1), NSLIP, NSLPTL, NSET, PROPS(97),
     2                   STATEV(NSLPTL+1),H)
         
C       Initial value of the kinametic hardening stress X         
         CALL Xslip(STATEV(9*NSLPTL+1),STATEV(10*NSLPTL+1), NSLIP, 
     2              NSLPTL, NSET, PROPS(97))

C-----  Initial value of the resolved shear stress in slip systems
         DO I=1,NSLPTL
            TERM1=0.

            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO

            STATEV(2*NSLPTL+I)=TERM1
         END DO

      ELSE
      
C-----  Current stress state
C
C-----  Copying from the array of state variables STATVE the following
C          parameters and variables at current stress state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Current slip directions SLPDIR
C          Normals to current slip planes SLPNOR
C
         NSLPTL=NINT(STATEV(NSTATV))
         DO I=1,NSET
            NSLIP(I)=NINT(STATEV(NSTATV-4+I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               SLPNOR(I,J)=STATEV(IDNOR)

               IDDIR=IDDIR+1
               SLPDIR(I,J)=STATEV(IDDIR)
            END DO
         END DO

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

      END IF

C-----  Slip spin tensor: SLPSPN (only needed for finite rotation)
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF

C-----  Double dot product of elastic moduli tensor with the slip 
C     deformation tensor (Schmid factors) plus, only for finite 
C     rotation, the dot product of slip spin tensor with the stress: 
C     DDEMSD
C
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF

C-----  Shear strain-rate in a slip system at the start of increment: 
C     FSLIP, and its derivative: DFDZSP
C
      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL STRAINRATE (STATEV(NSLPTL+ID), STATEV(2*NSLPTL+ID), 
     2                    STATEV(ID), NSLIP(I), FSLIP(ID), DFDZSP(ID), 
     3                    PROPS(65+8*I),STATEV(10*NSLPTL+ID))
      END DO
      

C-----  Kinematic and isotropic-hardening laws
      DO I=1,NSET
         ID=1
         IF (I.GT.1) ID=ID+NSLIP(I-1)      
         CALL DelKinematicHard(STATEV(9*NSLPTL+ID), NSLIP(I), FSLIP(ID),
     3                    PROPS(81+16*I),DelX,ID,NSLPTL)
      END DO  

      
      CALL DelIsotroicHard(STATEV(NSLPTL+1),NSET, NSLIP,H(1,1), 
     3                    PROPS(97),DelR,NSLPTL)


C-----  LU decomposition to solve the increment of shear strain in a 
C     slip system
C
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
!         TAUSLP=STATEV(2*NSLPTL+I)
!         RSLIP=STATEV(I)
!         Xslp=STATEV(10*NSLPTL+I)
!         Z=abs(TAUSLP-Xslp)-RSLIP
         TermForTau=TERM1*DFDZSP(I)
         TermForX=TERM1*DFDZSP(I)
         TermForR=-TERM1*DFDZSP(I)


         DO J=1,NSLPTL
            TERM4=0.
            DO K=1,6
               TERM4=TERM4+DDEMSD(K,I)*SLPDEF(K,J)
            END DO
            WORKST(I,J)=TermForTau*TERM4
C_HC  term for kinematic and isotropic hardening            
            IF (I.EQ.J) WORKST(I,J)=WORKST(I,J)-TermForX*DelX(I)
     2                             -TermForR*DelR(I)

         END DO

         WORKST(I,I)=WORKST(I,I)+1.
      END DO

      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP)

C-----  Increment of shear strain in a slip system: DGAMMA
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
!            TAUSLP=STATEV(2*NSLPTL+I)
!            RSLIP=STATEV(I)
!            Xslp=STATEV(10*NSLPTL+I)            
!            Z=abs(TAUSLP-Xslp)-RSLIP
            TERM2=TERM1*DFDZSP(I)

            DGAMMA(I)=0.
            DO J=1,NDI
               DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
            END DO
            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF
            DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME
      END DO

      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)


C-----  Update the shear strain in a slip system: STATEV(NSLPTL+1) - 
C     STATEV(2*NSLPTL)
C
      DO I=1,NSLPTL
         STATEV(NSLPTL+I)=STATEV(NSLPTL+I)+DGAMMA(I)
      END DO
      
C-----  Total cumulative shear strains on all slip systems (sum of the 
C       absolute values of shear strains in all slip systems)
CFIX--  Total cumulative shear strains on each slip system (sum of the 
CFIX    absolute values of shear strains in each individual slip system)
C
      DO I=1,NSLPTL
         STATEV(NSTATV-4)=STATEV(NSTATV-4)+ABS(DGAMMA(I))
         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))
      END DO
      
C-----  Updating DelX and DelR using DGAMMA
!      DO I=1,NSET
!         ID=1
!         IF (I.GT.1) ID=ID+NSLIP(I-1)      
!         CALL DelKinematicHard(STATEV(9*NSLPTL+ID), NSLIP(I), FSLIP(ID),
!     3                    PROPS(81+16*I),DelX,ID,NSLPTL)
!      END DO  
!      
!      CALL DelIsotroicHard(STATEV(NSLPTL+1),NSET, NSLIP,H(1,1), 
!     3                    PROPS(97),DelR,NSLPTL)      
      
      
      
C_HC updating stateV for X and r     
      DO I=1,NSLPTL
         STATEV(10*NSLPTL+I)=STATEV(10*NSLPTL+I)+DelX(I)*DGAMMA(I)
         STATEV(I)=STATEV(I)+DelR(I)*DGAMMA(I)
      END DO      

C-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

C-----  Increment of deformation gradient associated with lattice 
C     stretching in the current state, i.e. the velocity gradient 
C     (associated with lattice stretching) times the increment of time:
C     DVGRAD (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.
                  ELSE
                     TERM1=-1.
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                  END DO
               END IF

            END DO
         END DO

      END IF

C-----  Increment of resolved shear stress in a slip system: DTAUSP
      DO I=1,NSLPTL
         DTAUSP(I)=0.
         DO J=1,6
            DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
         END DO
      END DO

C-----  Update the resolved shear stress in a slip system: 
C     STATEV(2*NSLPTL+1) - STATEV(3*NSLPTL)
C
      DO I=1,NSLPTL
         STATEV(2*NSLPTL+I)=STATEV(2*NSLPTL+I)+DTAUSP(I)
      END DO

C-----  Increment of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.
         END DO
      ELSE
         DO I=1,NTENS
            DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF

      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR

            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
            END DO

         END DO
      END IF

C-----  Update the stress: STRESS
      DO I=1,NTENS
         STRESS(I)=STRESS(I)+DSTRES(I)
      END DO

C-----  Increment of normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.
               DSPDIR(I,J)=0.

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO

C-----  Update the normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)
            END DO
         END DO

      END IF

C-----  Derivative of shear strain increment in a slip system w.r.t. 
C     strain increment: DDGDDE
C
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
!            TAUSLP=STATEV(2*NSLPTL+J)
!            RSLIP=STATEV(J)
!            Xslp=STATEV(10*NSLPTL+I)            
!            Z=abs((TAUSLP-Xslp))-RSLIP            
            TERM2=TERM1*DFDZSP(J)
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO

C-----  Derivative of stress increment w.r.t. strain increment, i.e. 
C     Jacobian matrix
C
C-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF

C-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2                                DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2                            DDEMSD(I,K)*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2                            DDEMSD(J+3,K)*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF

C=========================================================================
C====== This commented section is not necessary for this UMAT 
C-----  Iteration ?
!      IF (ITRATN.NE.0) THEN

C-----  Save solutions (without iteration):
C            Shear strain-rate in a slip system FSLIP1
C            Current strength in a slip system GSLP1
C            Shear strain in a slip system GAMMA1
C            Resolved shear stress in a slip system TAUSP1
C            Normal to a slip plane SPNOR1
C            Slip direction SPDIR1
C            Stress STRES1
C            Jacobian matrix DDSDE1
C
!         IF (NITRTN.EQ.0) THEN
!
!            IDNOR=3*NSLPTL
!            IDDIR=6*NSLPTL
!            DO J=1,NSLPTL
!               FSLIP1(J)=FSLIP(J)
!               GSLP1(J)=STATEV(J)
!               GAMMA1(J)=STATEV(NSLPTL+J)
!               TAUSP1(J)=STATEV(2*NSLPTL+J)
!               DO I=1,3
!                  IDNOR=IDNOR+1
!                  SPNOR1(I,J)=STATEV(IDNOR)
!
!                  IDDIR=IDDIR+1
!                  SPDIR1(I,J)=STATEV(IDDIR)
!               END DO
!            END DO

!           DO J=1,NTENS
!               STRES1(J)=STRESS(J)
!               DO I=1,NTENS
!                  DDSDE1(I,J)=DDSDDE(I,J)
!               END DO
!            END DO

!         END IF

C-----  Increments of stress DSOLD, and solution dependent state 
C     variables DGAMOD, DTAUOD, DGSPOD, DSPNRO, DSPDRO (for the next 
C     iteration)
C
!        DO I=1,NTENS
!            DSOLD(I)=DSTRES(I)
!         END DO

!         DO J=1,NSLPTL
!            DGAMOD(J)=DGAMMA(J)
!            DTAUOD(J)=DTAUSP(J)
!            DGSPOD(J)=DGSLIP(J)
!            DO I=1,3
!               DSPNRO(I,J)=DSPNOR(I,J)
!               DSPDRO(I,J)=DSPDIR(I,J)
!            END DO
!         END DO

C-----  Check if the iteration solution converges
!         IDBACK=0
!         ID=0
!         DO I=1,NSET
!            DO J=1,NSLIP(I)
!               ID=ID+1
!               X=STATEV(2*NSLPTL+ID)/STATEV(ID)
!               RESIDU=THETA*DTIME*F(X,PROPS(65+8*I))+DTIME*(1.0-THETA)*
!     2                FSLIP1(ID)-DGAMMA(ID)
!               IF (ABS(RESIDU).GT.GAMERR) IDBACK=1
!            END DO
!         END DO

!         IF (IDBACK.NE.0.AND.NITRTN.LT.ITRMAX) THEN
C-----  Iteration: arrays for iteration
CFIXA
!            CALL ITERATION (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                      STATEV(1), STATEV(9*NSLPTL+1), 
     3                      STATEV(10*NSLPTL+1), NSLPTL, 
     4                      NSET, NSLIP, ND, PROPS(97), DGAMOD,
     5                      DHDGDG)
CFIXB

!            GO TO 1000

!         ELSE IF (NITRTN.GE.ITRMAX) THEN
C-----  Solution not converge within maximum number of iteration (the 
C     solution without iteration will be used)
C
!            DO J=1,NTENS
!               STRESS(J)=STRES1(J)
!               DO I=1,NTENS
!                  DDSDDE(I,J)=DDSDE1(I,J)
!               END DO
!            END DO

 !           IDNOR=3*NSLPTL
!            IDDIR=6*NSLPTL
!            DO J=1,NSLPTL
!               STATEV(J)=GSLP1(J)
!               STATEV(NSLPTL+J)=GAMMA1(J)
!               STATEV(2*NSLPTL+J)=TAUSP1(J)

!               DO I=1,3
!                  IDNOR=IDNOR+1
!                  STATEV(IDNOR)=SPNOR1(I,J)

!                  IDDIR=IDDIR+1
!                  STATEV(IDDIR)=SPDIR1(I,J)
!               END DO
!            END DO
!
!         END IF

!      END IF

C-----  Total cumulative shear strains on all slip systems (sum of the 
C       absolute values of shear strains in all slip systems)
CFIX--  Total cumulative shear strains on each slip system (sum of the 
CFIX    absolute values of shear strains in each individual slip system)
C
!      DO I=1,NSLPTL
CFIXA
!         STATEV(10*NSLPTL+1)=STATEV(10*NSLPTL+1)+ABS(DGAMMA(I))
!         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))
CFIXB
!      END DO
	  
CHelal_ Updating SLPDEF
C-----  Updating SLPDEF (Schmid factors)
!      DO J=1,NSLPTL
!         SLPDEF(1,J)=STATEV(6*NSLPTL+3*J-2)*STATEV(3*NSLPTL+3*J-2)
!         SLPDEF(2,J)=STATEV(6*NSLPTL+3*J-1)*STATEV(3*NSLPTL+3*J-1)
!         SLPDEF(3,J)=STATEV(6*NSLPTL+3*J)*STATEV(3*NSLPTL+3*J)
!         SLPDEF(4,J)=STATEV(6*NSLPTL+3*J-2)*STATEV(3*NSLPTL+3*J-1)+
!     1                  STATEV(6*NSLPTL+3*J-1)*STATEV(3*NSLPTL+3*J-2)
!         SLPDEF(5,J)=STATEV(6*NSLPTL+3*J-2)*STATEV(3*NSLPTL+3*J)+
!     2                  STATEV(6*NSLPTL+3*J)*STATEV(3*NSLPTL+3*J-2)
!         SLPDEF(6,J)=STATEV(6*NSLPTL+3*J-1)*STATEV(3*NSLPTL+3*J)+
!     3                  STATEV(6*NSLPTL+3*J)*STATEV(3*NSLPTL+3*J-1)
!      END DO
CHelal_ plastic strain calculate kore store korar jonno SDV(122-127)
!     DO I=1,6
!        STATEV(10*NSLPTL+1+I)=0.0
!        DO J=1,NSLPTL
!         STATEV(10*NSLPTL+1+I)= STATEV(10*NSLPTL+1+I)+
!     1                          DGAMMA(J)*SLPDEF(I,J)
!        END DO 
!        STATEV(10*NSLPTL+1+I)=STATEV(10*NSLPTL+1+I)+ eps_PrevStep(I)
!      END DO
C_Helal_ necessary to store plastic strain of the previous time step
!     DO I=1,6
!         eps_PrevStep(I)=statev(ND*10+1+I)
!      END DO
C_Helal_ necessary to store plastic strain of the previous time step
CHelal_ plastic strain calculate kore store korar jonno SDV(122-127)

C=========================================================================
      RETURN
      END
C----------------------------------------------------------------------

      SUBROUTINE ROTATION (PROP, ROTATE)

C-----  This subroutine calculates the rotation matrix, i.e. the 
C     direction cosines of cubic crystal [100], [010] and [001] 
C     directions in global system

C-----  The rotation matrix is stored in the array ROTATE.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3) 

C-----  Subroutines:
C
C       CROSS  -- cross product of two vectors
C
C       LUDCMP -- LU decomposition
C
C       LUBKSB -- linear equation solver based on LU decomposition 
C                 method (must call LUDCMP first)


C-----  PROP -- constants characterizing the crystal orientation 
C               (INPUT)
C
C            PROP(1) - PROP(3) -- direction of the first vector in 
C                                 local cubic crystal system
C            PROP(4) - PROP(6) -- direction of the first vector in 
C                                 global system
C
C            PROP(9) - PROP(11)-- direction of the second vector in 
C                                 local cubic crystal system
C            PROP(12)- PROP(14)-- direction of the second vector in 
C                                 global system
C
C-----  ROTATE -- rotation matrix (OUTPUT):
C
C            ROTATE(i,1) -- direction cosines of direction [1 0 0] in 
C                           local cubic crystal system
C            ROTATE(i,2) -- direction cosines of direction [0 1 0] in 
C                           local cubic crystal system
C            ROTATE(i,3) -- direction cosines of direction [0 0 1] in 
C                           local cubic crystal system
	  
	  
C-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)

C-----  LU decomposition of TERM1
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP)

C-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.
            ELSE
               TERM2(I,J)=0.
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

C-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)

C-----  Check: the angle between first and second vector in local and 
C     global systems must be the same.  The relative difference must be
C     less than 0.1%.
C
      IF (ABS(ANGLE1/ANGLE2-1.).GT.0.001) THEN 
         WRITE (6,*) 
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF

C-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END
C----------------------------------------------------------------------

           SUBROUTINE CROSS (A, B, C, ANGLE)

C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*) 
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF

           RETURN
           END
C----------------------------------------------------------------------


      SUBROUTINE SLIPSYS (ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR, 
     2                    ROTATE)

C-----  This subroutine generates all slip systems in the same set for 
C     a CUBIC crystal.  For other crystals (e.g., HCP, Tetragonal, 
C     Orthotropic, ...), it has to be modified to include the effect of
C     crystal aspect ratio.

C-----  Denote s as a slip direction and m as normal to a slip plane.  
C     In a cubic crystal, (s,-m), (-s,m) and (-s,-m) are NOT considered
C     independent of (s,m).

C-----  Subroutines:  LINE1 and LINE
	 
C-----  Variables:
C
C     ISPDIR -- a typical slip direction in this set of slip systems 
C               (integer)  (INPUT)
C     ISPNOR -- a typical normal to slip plane in this set of slip 
C               systems (integer)  (INPUT)
C     NSLIP  -- number of independent slip systems in this set 
C               (OUTPUT)
C     SLPDIR -- unit vectors of all slip directions  (OUTPUT)
C     SLPNOR -- unit normals to all slip planes  (OUTPUT)
C     ROTATE -- rotation matrix (INPUT)
C          ROTATE(i,1) -- direction cosines of [100] in global system
C          ROTATE(i,2) -- direction cosines of [010] in global system
C          ROTATE(i,3) -- direction cosines of [001] in global system
C
C     NSPDIR -- number of all possible slip directions in this set
C     NSPNOR -- number of all possible slip planes in this set
C     IWKDIR -- all possible slip directions (integer)
C     IWKNOR -- all possible slip planes (integer)


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISPDIR(3), ISPNOR(3), SLPDIR(3,12), SLPNOR(3,12), 
     *          ROTATE(3,3), IWKDIR(3,12), IWKNOR(3,12), TERM(3)

      NSLIP=0
      NSPDIR=0
      NSPNOR=0
      
C_HC modified for 14 slip system in (111)[110] and  (001)[110]  
C-----  Generating all possible slip directions in this set
      I=IABS(ISPNOR(1))+IABS(ISPNOR(2))+IABS(ISPNOR(3))

      IF (I.EQ.3.) THEN
C        ---(111)[110] type
          NSPDIR=12
          NSPNOR=12
          RMODIR=sqrt(2.0)
          RMONOR=sqrt(3.0)
          
          IWKDIR(1,1)=0.
          IWKDIR(2,1)=-1.
          IWKDIR(3,1)=1.
          IWKDIR(1,2)=1.
          IWKDIR(2,2)=0.
          IWKDIR(3,2)=-1.
          IWKDIR(1,3)=-1.
          IWKDIR(2,3)=1.
          IWKDIR(3,3)=0.
          IWKDIR(1,4)=1.
          IWKDIR(2,4)=0.
          IWKDIR(3,4)=1.
          IWKDIR(1,5)=1.
          IWKDIR(2,5)=1.
          IWKDIR(3,5)=0.
          IWKDIR(1,6)=0.
          IWKDIR(2,6)=-1.
          IWKDIR(3,6)=1.
          IWKDIR(1,7)=0.
          IWKDIR(2,7)=1.
          IWKDIR(3,7)=1.
          IWKDIR(1,8)=1.
          IWKDIR(2,8)=1.
          IWKDIR(3,8)=0.
          IWKDIR(1,9)=1.
          IWKDIR(2,9)=0.
          IWKDIR(3,9)=-1.
          IWKDIR(1,10)=0.
          IWKDIR(2,10)=1.
          IWKDIR(3,10)=1.
          IWKDIR(1,11)=1.
          IWKDIR(2,11)=0.
          IWKDIR(3,11)=1.
          IWKDIR(1,12)=-1.
          IWKDIR(2,12)=1.
          IWKDIR(3,12)=0.
          
          IWKNOR(1:3,1:3)=1.0
          IWKNOR(1,4:6)=-1.0
          IWKNOR(2:3,4:6)=1.0
          IWKNOR(1,7:9)=1.0
          IWKNOR(2,7:9)=-1.0  
          IWKNOR(3,7:9)=1.0 
          IWKNOR(1:2,10:12)=1.0
          IWKNOR(3,10:12)=-1.0

      ELSE 
        IF (I.EQ.1.) THEN
C        ---(001)[110] type
          NSPDIR=2
          NSPNOR=2
          RMODIR=sqrt(2.0)
          RMONOR=sqrt(1.0)
          
          IWKDIR(1,1)=1.
          IWKDIR(2,1)=1.
          IWKDIR(1,2)=1.
          IWKDIR(2,2)=-1.

          IWKNOR(3,1:2)=1.

        END IF
      END IF

C-----  Printing all slip systems in this set
C
C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in local cubic crystal system
C
      WRITE (6,*) '          '
      WRITE (6,*) ' #          Slip plane          Slip direction'
  
         DO I=1,NSPDIR
               NSLIP=NSLIP+1
               DO K=1,3
                  SLPDIR(K,NSLIP)=IWKDIR(K,I)/RMODIR
                  SLPNOR(K,NSLIP)=IWKNOR(K,I)/RMONOR
               END DO

               WRITE (6,10) NSLIP, 
     2                      (IWKNOR(K,I),K=1,3), (IWKDIR(K,I),K=1,3)

         END DO

10    FORMAT(1X,I2,9X,'(',3(1X,I2),1X,')',10X,'[',3(1X,I2),1X,']')

      WRITE (6,*) 'Number of slip systems in this set = ',NSLIP
      WRITE (6,*) '          '

      IWKDIR=0.0
      IWKNOR=0.0
C_HC modified for 14 slip system in (111)[110] and  (001)[110]    
          
      IF (NSLIP.EQ.0) THEN
         WRITE (6,*) 
     *      'There is no slip direction normal to the slip planes!'
         STOP

      ELSE

C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in global system
C
         DO J=1,NSLIP
            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPDIR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPDIR(I,J)=TERM(I)
            END DO

            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPNOR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPNOR(I,J)=TERM(I)
            END DO
         END DO

      END IF

      RETURN
      END
C----------------------------------------------------------------------

      SUBROUTINE rSlipInit (rSLIP, NSLIP, NSLPTL, NSET, PROP,Gamma,H)

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL rSLIP0
      DIMENSION rSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET),Gamma(NSLPTL),
     2           H(NSLPTL,NSLPTL),Term0(NSLPTL),Term1(NSLPTL)
      
      DO I=1,NSLPTL
         b=PROP(7,1)
         if (I.GT.NSLIP(1)) b=PROP(7,2)
         Term0(I)=1-exp(-b*abs(Gamma(I)))
         DO J=1,NSLPTL
            if (I.EQ.J) then
               H(I,J)=1.0
            else if(I.GT.NSLIP(1).OR.J.GT.NSLIP(1)) then 
                  H(I,J)=PROP(10,1)
                  if (I.GT.NSLIP(1).AND.J.GT.NSLIP(1)) H(I,J)=PROP(9,1) 
            else
                  H(I,J)=PROP(9,1)
            end if
         END DO
      END DO
      
      DO I=1,NSLPTL
         Term1(I)=0.0
        DO J=1,NSLPTL
          Term1(I)= Term1(I)+ H(I,J)* Term0(J)
        END DO
      END DO
                  
      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            rSLIP(ID)=rSLIP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET,
     2                       Term1(ID))
         END DO
      END DO

      RETURN
      END
C----------------------------------------------------------------------

       REAL*8 FUNCTION rSLIP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET,Term1)

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION NSLIP(NSET), PROP(16)

       rSLIP0=PROP(5)+PROP(6)*Term1

       RETURN
       END
C----------------------------------------------------------------------

      SUBROUTINE Xslip(vSS,XSLIP0, NSLIP, NSLPTL, NSET, PROP)

C----- Kinematic hardening variable X could either be the same for all slip systems 
C     in each set, or be different from set to set, e.g.
C     <110>{111} and <110>{100}.
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL XSLP0
      DIMENSION XSLIP0(NSLPTL), NSLIP(NSET),PROP(16,NSET),vSS(NSLPTL)

      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            XSLIP0(ID)=XSLP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET,
     2                 vSS(ID))
         END DO
      END DO

      RETURN
      END
C----------------------------------------------------------------------

           REAL*8 FUNCTION XSLP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET,vSlp)
C-----     User-supplied function subprogram given the initial value of
C         Kinematic hardening variable at initial state

           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION NSLIP(NSET), PROP(16)

           phiV=PROP(13)+(1-PROP(13))*exp(-PROP(14)*abs(vSlp))
           Term1=phiV*DSIGN(1.D0,vSlp)/PROP(12)
           Term2=PROP(14)*(1-PROP(13))*exp(-PROP(14)*abs(vSlp))
           Term3=Term2/(PROP(12)*(PROP(12)-PROP(14)*DSIGN(1.D0,vSlp)))
           ASLP0=Term1+Term3
           XSLP0=PROP(11)*ASLP0
           
           RETURN
           END
C----------------------------------------------------------------------

      SUBROUTINE STRAINRATE (GAMMA, TAUSLP, Rslip, NSLIP, FSLIP, 
     2                       DFDZSP, PROP,Xslip)

C-----  This subroutine calculates the shear strain-rate in each slip 
C     system for a rate-dependent single crystal.  The POWER LAW 
C     relation between shear strain-rate and resolved shear stress 
C     proposed by Hutchinson, Pan and Rice, is used here.

C-----  The power law exponents are assumed the same for all slip 
C     systems in each set, though they could be different from set to 
C     set, e.g. <110>{111} and <110>{100}.  The strain-rate coefficient
C     in front of the power law form are also assumed the same for all 
C     slip systems in each set. 

C-----  Users who want to use their own constitutive relation may 
C     change the function subprograms F and its derivative DFDX, 
C     where F is the strain hardening law, dGAMMA/dt = F(X), 
C     X=TAUSLP/GSLIP.  The parameters characterizing F are passed into 
C     F and DFDX through array PROP.

C-----  Function subprograms:
C
C       F    -- User-supplied function subprogram which gives shear 
C               strain-rate for each slip system based on current 
C               values of resolved shear stress and current strength
C
C       DFDX -- User-supplied function subprogram dF/dX, where x is the
C               ratio of resolved shear stress over current strength

C-----  Variables:
C
C     GAMMA  -- shear strain in each slip system at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in each slip system (INPUT)
C     GSLIP  -- current strength (INPUT)
C     NSLIP  -- number of slip systems in this set (INPUT)
C
C     FSLIP  -- current value of F for each slip system (OUTPUT)
C     DFDXSP -- current value of DFDX for each slip system (OUTPUT)
C
C     PROP   -- material constants characterizing the strain hardening 
C               law (INPUT)
C
C               For the current power law strain hardening law 
C               PROP(1) -- power law hardening exponent
C               PROP(1) = infinity corresponds to a rate-independent 
C               material
C               PROP(2) -- coefficient in front of power law hardening


      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F, DFDZ
      DIMENSION GAMMA(NSLIP), TAUSLP(NSLIP), Rslip(NSLIP), 
     2          FSLIP(NSLIP), DFDZSP(NSLIP), PROP(8),Xslip(NSLIP)

      DO I=1,NSLIP
         sgn=DSIGN(1.D0,(TAUSLP(I)-Xslip(I)))
         Z=abs(TAUSLP(I)-Xslip(I))-Rslip(I)
         FSLIP(I)=F(Z,PROP,sgn)
         DFDZSP(I)=DFDZ(Z,PROP,sgn)
      END DO

      RETURN
      END
C----------------------------------------------------------------------

           REAL*8 FUNCTION F(Z,PROP,sgn)

C-----     User-supplied function subprogram which gives shear 
C        strain-rate for each slip system based on current values of 
C        resolved shear stress and current strength

           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)
           
           if (Z.GE.0.) then
              F=(Z/PROP(2))**PROP(1)*sgn
           else
              F=0.0
           end if   

           RETURN
           END
C----------------------------------------------------------------------

           REAL*8 FUNCTION DFDZ(Z,PROP,sgn)

C-----     User-supplied function subprogram dF/dX, where x is the 
C        ratio of resolved shear stress over current strength

           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)
           
           if (Z.GE.0.) then
               DFDZ=(PROP(1)/PROP(2))*(Z/PROP(2))**(PROP(1)-1.)
c			   *sgn
           else
              DFDZ=0.0
           end if
           	
	     RETURN
           END
C----------------------------------------------------------------------


      SUBROUTINE DelKinematicHard(GmmaVslip,NSLIP, FSLP,
     3                    PROP, DelX,ID,NSLPTL)
   
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION GmmaVslip(NSLIP),PROP(16),
     2               DelX(NSLPTL),FSLP(NSLIP)
     
        Term1=PROP(11)*PROP(14)*(1-PROP(13))
        DO I=1,NSLIP 
           DelX(ID)=0.0      
           term2=PROP(14)-PROP(12)*DSIGN(1.D0,FSLP(I))
           DelX(ID)=Term1*exp(-PROP(14)*GmmaVslip(I))/term2
           ID=ID+1
        End do   
     
      RETURN
      END
C---------------------------------------------------------------------- 

      SUBROUTINE DelIsotroicHard(Gmma,NSET, NSLP,H, 
     3                    PROP,DelR,NSLPTL)   

     
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Gmma(NSLPTL),PROP(48),H(NSLPTL,NSLPTL),Qb(NSLPTL),
     2               DelR(NSLPTL),NSLP(NSET),Term0(NSLPTL)
          
      DO I=1,NSLPTL
         DO J=1,NSLPTL
            if (I.EQ.J) then
               H(I,J)=1.0
            else if(I.GT.NSLP(1).OR.J.GT.NSLP(1)) then 
                  H(I,J)=PROP(10)
                  if (I.GT.NSLP(1).AND.J.GT.NSLP(1)) H(I,J)=PROP(9) 
            else
                  H(I,J)=PROP(9)
            end if
         END DO
      END DO

      ID=0
      DO I=1,NSET
        DO J=1,NSLP(I)
         ID=ID+1
         b=PROP(16*I-9)
         Term0(ID)=1.0-exp(-b*abs(Gmma(ID)))
         Qb(ID)=PROP(16*I-10)*b
        END DO 
      END DO
      
      DO I=1,NSLPTL
         DelR(I)=0.0
        DO J=1,NSLPTL
          DelR(I)= DelR(I)+ H(I,J)* Term0(J)
        END DO
        DelR(I)=DelR(I)*Qb(I)
      END DO      
      
      RETURN
      END
C----------------------------------------------------------------------


      SUBROUTINE LUDCMP(A, N, NP, INDX, D)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.
      DO I=1,N
         AAMAX=0.

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

      RETURN
      END
C----------------------------------------------------------------------


      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C========= Following functions are not necessary for this UMAT formulation

!           SUBROUTINE LINE (I1, I2, I3, IARRAY)

C-----  Generating all possible slip directions <lmn> (or slip planes 
C     {lmn}) for a cubic crystal, where l,m,n are not zeros.

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
!           DIMENSION IARRAY(3,4)

!           DO J=1,4
!              IARRAY(1,J)=I1
!              IARRAY(2,J)=I2
!              IARRAY(3,J)=I3
!           END DO

!           DO I=1,3
!              DO J=1,4
!                 IF (J.EQ.I+1) IARRAY(I,J)=-IARRAY(I,J)
!              END DO
!           END DO

!           RETURN
!          END


C-----------------------------------


!           SUBROUTINE LINE1 (J1, J2, IARRAY, ID)

C-----  Generating all possible slip directions <0mn> (or slip planes 
C     {0mn}) for a cubic crystal, where m,n are not zeros and m does 
C     not equal n.

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
!           DIMENSION IARRAY(3,2)

!          IARRAY(ID,1)=0
!          IARRAY(ID,2)=0

!           ID1=ID+1
!           IF (ID1.GT.3) ID1=ID1-3
!           IARRAY(ID1,1)=J1
!           IARRAY(ID1,2)=J1

!           ID2=ID+2
!           IF (ID2.GT.3) ID2=ID2-3
!           IARRAY(ID2,1)=J2
!           IARRAY(ID2,2)=-J2
!  
!           RETURN
!           END



C----------------------------------------------------------------------

CFIXA
!      SUBROUTINE LATENTHARDEN (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
!     2                         NSLIP, NSLPTL, NSET, H, PROP, ND)
CFIXB

C-----  This subroutine calculates the current self- and latent-
C     hardening moduli for all slip systems in a rate-dependent single 
C     crystal.  Two kinds of hardening law are used here.  The first 
C     law, proposed by Asaro, and Pierce et al, assumes a HYPER SECANT 
C     relation between self- and latent-hardening moduli and overall 
C     shear strain.  The Bauschinger effect has been neglected.  The 
C     second is Bassani's hardening law, which gives an explicit 
C     expression of slip interactions between slip systems.  The 
C     classical three stage hardening for FCC single crystal could be 
C     simulated.

C-----  The hardening coefficients are assumed the same for all slip 
C     systems in each set, though they could be different from set to 
C     set, e.g. <110>{111} and <110>{100}.

C-----  Users who want to use their own self- and latent-hardening law 
C     may change the function subprograms HSELF (self hardening) and 
C     HLATNT (latent hardening).  The parameters characterizing these 
C     hardening laws are passed into HSELF and HLATNT through array 
C     PROP.


C-----  Function subprograms:
C
C       HSELF  -- User-supplied self-hardening function in a slip 
C                 system
C
C       HLATNT -- User-supplied latent-hardening function

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip system 
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems 
C               (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C
C     H      -- current value of self- and latent-hardening moduli 
C               (OUTPUT)
C               H(i,i) -- self-hardening modulus of the ith slip system
C                         (no sum over i)
C               H(i,j) -- latent-hardening molulus of the ith slip 
C                         system due to a slip in the jth slip system 
C                         (i not equal j)
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in 
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set and jth set (i not equal j) 
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip 
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C     ND     -- leading dimension of arrays defined in subroutine UMAT 
C               (INPUT) 


C-----  Use single precision on cray
C
!      IMPLICIT REAL*8 (A-H,O-Z)
!      EXTERNAL HSELF, HLATNT
CFIXA
!      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
!     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
!     3          H(ND,NSLPTL)
CFIXB

!      CHECK=0.
!      DO I=1,NSET
!         DO J=4,8
!            CHECK=CHECK+ABS(PROP(J,I))
!         END DO
!      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

!      ISELF=0
!      DO I=1,NSET
!         ISET=I
!         DO J=1,NSLIP(I)
!            ISELF=ISELF+1

!            DO LATENT=1,NSLPTL
!               IF (LATENT.EQ.ISELF) THEN
CFIXA
!                  H(LATENT,ISELF)=HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,
!     2                                  NSET,NSLIP,PROP(1,I),CHECK,
!     3                                  ISELF,ISET)
CFIXB
!               ELSE
CFIXA
!                  H(LATENT,ISELF)=HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,
!     2                                   NSET,NSLIP,PROP(1,I),CHECK,
!     3                                   ISELF,ISET,LATENT)
CFIXB

!               END IF
!            END DO

!         END DO
!      END DO

!      RETURN
!      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
!           REAL*8 FUNCTION HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                           NSLIP,PROP,CHECK,ISELF,ISET)
CFIXB

C-----     User-supplied self-hardening function in a slip system

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
!           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
!     2               GMSLTL(NSLPTL)
CFIXB

!           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
!              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              HSELF=PROP(1)*TERM2**2

!           ELSE

C-----  Bassani's hardening law
CFIXA
!              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

!              ID=0
!              G=1.
!              DO I=1,NSET
!                 IF (I.EQ.ISET) THEN
!                    GAMMA0=PROP(5)
!                    FAB=PROP(7)
!                 ELSE
!                    GAMMA0=PROP(6)
!                    FAB=PROP(8)
!                 END IF

!                 DO J=1,NSLIP(I)
!                    ID=ID+1
!                    IF (ID.NE.ISELF) THEN
CFIXA
!		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
!		    END IF

!                 END DO
!              END DO

!              HSELF=F*G

!           END IF

!           RETURN
!           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
!           REAL*8 FUNCTION HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT)
CFIXB

C-----     User-supplied latent-hardening function

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
!           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
!     2               GMSLTL(NSLPTL)
CFIXB

!           Q=PROP(10)
!           d=abs(LATENT-ISELF)
!           if (d.LE.2) then
!             IF (mod(ISELF,3).EQ.1.AND.LATENT.GT.ISELF) then
!                   Q=PROP(9)
!             end if    
!             IF (mod(ISELF,3).EQ.2.AND.d.EQ.1)then 
!                  Q=PROP(9)
!             end if   
!             IF (mod(ISELF,3).EQ.0.AND.LATENT.LT.ISELF) then
!                  Q=PROP(9)   
!             end if
!           end if

!           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
!              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              HLATNT=PROP(1)*TERM2**2*Q

!           ELSE

C-----  Bassani's hardening law
CFIXA
!              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

!              ID=0
!              G=1.
!              DO I=1,NSET
!                 IF (I.EQ.ISET) THEN
!                    GAMMA0=PROP(5)
!                    FAB=PROP(7)
!                 ELSE
!                    GAMMA0=PROP(6)
!                    FAB=PROP(8)
!                 END IF

!                 DO J=1,NSLIP(I)
!                    ID=ID+1
!                    IF (ID.NE.ISELF) THEN
CFIXA
!		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
!		    END IF

!                 END DO
!              END DO

!              HLATNT=F*G*Q

!           END IF

!           RETURN
!           END


C----------------------------------------------------------------------

CFIXA
 !     SUBROUTINE ITERATION (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
 !    2                      NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD, 
 !    3                      DHDGDG)
CFIXB

C-----  This subroutine generates arrays for the Newton-Rhapson 
C     iteration method.

C-----  Users who want to use their own self- and latent-hardening law 
C     may change the function subprograms DHSELF (self hardening) and 
C     DHLATN (latent hardening).  The parameters characterizing these 
C     hardening laws are passed into DHSELF and DHLATN through array 
C     PROP.


C-----  Function subprograms:
C
C       DHSELF -- User-supplied function of the derivative of self-
C                 hardening moduli
C
C       DHLATN -- User-supplied function of the derivative of latent-
C                 hardening moduli

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip system 
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems 
C               (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     ND     -- leading dimension of arrays defined in subroutine UMAT 
C               (INPUT) 
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in 
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set and jth set (i not equal j) 
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip 
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C-----  Arrays for iteration:
C
C       DGAMOD (INPUT)
C
C       DHDGDG (OUTPUT)
C

C-----  Use single precision on cray
C
!      IMPLICIT REAL*8 (A-H,O-Z)
!      EXTERNAL DHSELF, DHLATN
CFIXA
!      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          DGAMOD(NSLPTL), DHDGDG(ND,NSLPTL)
CFIXB

!      CHECK=0.
!      DO I=1,NSET
!         DO J=4,8
!            CHECK=CHECK+ABS(PROP(J,I))
!         END DO
!      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

!      ISELF=0
!      DO I=1,NSET
!         ISET=I
!         DO J=1,NSLIP(I)
!            ISELF=ISELF+1

!            DO KDERIV=1,NSLPTL
!               DHDGDG(ISELF,KDERIV)=0.

!               DO LATENT=1,NSLPTL
!                  IF (LATENT.EQ.ISELF) THEN
CFIXA
!                     DHDG=DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
!     3                           KDERIV)
CFIXB
!                  ELSE
CFIXA
!                     DHDG=DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
!     3                           LATENT,KDERIV)
CFIXB
!                  END IF

!                  DHDGDG(ISELF,KDERIV)=DHDGDG(ISELF,KDERIV)+
!     2                                 DHDG*ABS(DGAMOD(LATENT))
!               END DO

!            END DO
!         END DO
!      END DO

!      RETURN
!      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
!           REAL*8 FUNCTION DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                            NSLIP,PROP,CHECK,ISELF,ISET,
!     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of self-hardening
C     moduli

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
!           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), 
!     2               NSLIP(NSET), PROP(16)
CFIXB

!           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
!              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
!              DHSELF=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3

!           ELSE

C-----  Bassani's hardening law
CFIXA
!              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

!              IF (KDERIV.EQ.ISELF) THEN
!                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
!                 ID=0
!                 G=1.
!                 DO I=1,NSET
!                    IF (I.EQ.ISET) THEN
!                       GAMMA0=PROP(5)
!                       FAB=PROP(7)
!                    ELSE
!                       GAMMA0=PROP(6)
!                       FAB=PROP(8)
!                    END IF

!                    DO J=1,NSLIP(I)
!                       ID=ID+1
CFIXA
!                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
!                    END DO
!                 END DO

!              ELSE
!                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
!                 ILOWER=0
!                 IUPPER=NSLIP(1)
!                 IF (ISET.GT.1) THEN
!                    DO K=2,ISET
!                       ILOWER=ILOWER+NSLIP(K-1)
!                       IUPPER=IUPPER+NSLIP(K)
!                    END DO
!                 END IF

!                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
!                    GAMMA0=PROP(5)
!                    FAB=PROP(7)
!                 ELSE
!                    GAMMA0=PROP(6)
!                    FAB=PROP(8)
!                 END IF

CFIXA
!                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
!                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
!                 G=FAB/GAMMA0*TERM5**2

!              END IF

!              DHSELF=F*G

!           END IF

!           RETURN
!           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
!           REAL*8 FUNCTION DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
!     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT,
!     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of latent-hardening 
C     moduli

C-----  Use single precision on cray
C
!           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
!           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), NSLIP(NSET), 
!     2               PROP(16)
CFIXB

!           Q=PROP(10)
!           d=abs(LATENT-ISELF)
!           if (d.LE.2) then
!             IF (mod(ISELF,3).EQ.1.AND.LATENT.GT.ISELF) then
!                   Q=PROP(9)
!             end if    
!             IF (mod(ISELF,3).EQ.2.AND.d.EQ.1)then 
!                  Q=PROP(9)
!             end if   
!             IF (mod(ISELF,3).EQ.0.AND.LATENT.LT.ISELF) then
!                  Q=PROP(9)   
!             end if
!           end if

!           IF (CHECK.EQ.0.) THEN
C-----  HYPER SECANT hardening law by Asaro, Pierce et al
!              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
!              DHLATN=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3*Q

!           ELSE
C-----  Bassani's hardening law
CFIXA
!              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
!              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
!              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

!              IF (KDERIV.EQ.ISELF) THEN
!                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
!                 ID=0
!                 G=1.
!                 DO I=1,NSET
!                    IF (I.EQ.ISET) THEN
!                       GAMMA0=PROP(5)
!                       FAB=PROP(7)
!                    ELSE
!                       GAMMA0=PROP(6)
!                       FAB=PROP(8)
!                    END IF

!                    DO J=1,NSLIP(I)
!                       ID=ID+1
CFIXA
!                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
!                    END DO
!                 END DO

!              ELSE
!                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
!                 ILOWER=0
!                 IUPPER=NSLIP(1)
!                 IF (ISET.GT.1) THEN
!                    DO K=2,ISET
!                       ILOWER=ILOWER+NSLIP(K-1)
!                       IUPPER=IUPPER+NSLIP(K)
!                    END DO
!                 END IF

!                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
!                    GAMMA0=PROP(5)
!                    FAB=PROP(7)
!                 ELSE
!                    GAMMA0=PROP(6)
!                    FAB=PROP(8)
!                 END IF
CFIXA
!                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
!                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
!                 G=FAB/GAMMA0*TERM5**2

!              END IF

!              DHLATN=F*G*Q

!           END IF

!           RETURN
!           END


C----------------------------------------------------------------------