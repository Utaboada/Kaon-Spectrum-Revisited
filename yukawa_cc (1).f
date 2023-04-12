    



         PROGRAM YUKAWA
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)


         DOUBLE PRECISION eta,ms,mu,md,mc,mb,PAULI,GM
         DOUBLE PRECISION lamda_pi,lamda_sigma,lamda_eta,lamda_k
         double precision INTEG,gch2,m_pi,m_sigma,m_eta,m_k
         DOUBLE PRECISION KM_PI,KM_SIGMA,KM_K,KM_ETA,KM_SO
         DOUBLE PRECISION LAMDAO,MUO,JP,LS_P,bb,adim
         DOUBLE PRECISION ac_conf,muc_con,Dlt_con,as_conf,rg,LL
         INTEGER J,nmax,l,I,k,spin,JJ

         INTEGER INFO, ITYPE, LDA, LDB, LDVL, LDVR, LWORK
         CHARACTER JOBVL,JOBVR
         complex*16 autov
        
         DOUBLE PRECISION  r,rmax,r1,pi,XFACT,jn
         DOUBLE PRECISION mc2,hc,tetap,HB
      
    
         
*     Algunas constantes
         PARAMETER (PI = 3.141592d0)
         PARAMETER (HB= 6.582119569E-22) !MeV·s
         PARAMETER (HC = 197.326963) !MeV·fm
         PARAMETER (bb = 0.518d0 )  ! fm -> ESCALA TIPICA PARA ADIMENSIONALIZAR
         PARAMETER (adim = bb/HC ) !MeV^{-1}
*     Parametros del modelo
*        masas de los quarks ligeros
         PARAMETER (MU=313.d0*adim) !adim
         PARAMETER (MS=555.d0*adim)  !adim
         PARAMETER (MD=313.D0*adim)  !adim
         PARAMETER (MC2 = (MU*MS)/(MU+MS)) !masa reducida (adim)

*        OGE
         PARAMETER (ALPHAO=2.118) !adim
         PARAMETER (LAMDAO=0.113*bb) !adim
         PARAMETER (MUO=36.976*adim) !adim
         PARAMETER (RO=0.181/bb) !adim
         PARAMETER (rg = 0.259/bb) !adim
         PARAMETER (LL=-16./3D0) 

*        Bosones de Goldstone
         PARAMETER (M_PI=0.70*bb) !adim
         PARAMETER (M_SIGMA=3.42*bb) !adim
         PARAMETER (M_ETA=2.51*bb) !adim
         PARAMETER (M_K=2.77*bb) ! adim
         PARAMETER (LAMDA_PI=4.20*bb) !adim
         PARAMETER (LAMDA_ETA=5.20*bb) !adim
         PARAMETER (LAMDA_SIGMA=4.20*bb) !adim
         PARAMETER (LAMDA_K=4.21*bb) !adim               
         PARAMETER (gch2=0.54) !ya va con el 4pi
         PARAMETER (tetap=(-15*PI)/180.D0) !RAD

         
*     Confinamiento  
         PARAMETER (ac_conf=507.4*adim) !adim
         PARAMETER (muc_con=0.576*bb) !adim
         PARAMETER (Dlt_con=184.432*adim) !adim
         PARAMETER (as_conf=0.81d0) !adim
*     Parametros deL GEM
         PARAMETER (nmax = 24)
         PARAMETER (r1 = 0.024d0/bb) !adim
         PARAMETER (rmax = 7.5d0/bb) !adim
*     Parametros para los eigenvalues
         parameter (LDVL=nmax)
         parameter (LDVR=nmax)
         PARAMETER (LWORK=8*nmax) 
         PARAMETER (NFIRST=10)





*      LLAPACK
         DOUBLE PRECISION A(nmax,nmax),BETA(nmax),VL(nmax,nmax)
         DOUBLE PRECISION VR(1:nmax,1:nmax),ALPHAR(nmax),ALPHAI(nmax)
         DOUBLE PRECISION WORK(LWORK),S(1:nmax,1:nmax)
         DOUBLE PRECISION T(1:nmax,1:nmax),V(1:nmax,1:nmax)
         DOUBLE PRECISION XBETA,NORM(nmax)

*      TERMINOS CENTRALES
         DOUBLE PRECISION VC_PI(1:NMAX,1:NMAX),VC_SIGMA(1:NMAX,1:NMAX)
         DOUBLE PRECISION VC_K(1:NMAX,1:NMAX),VC_ETA(1:NMAX,1:NMAX)
         DOUBLE PRECISION VQQ_C(1:NMAX,1:NMAX),VC_OGE(1:NMAX,1:NMAX)
         DOUBLE PRECISION VC_CON(1:NMAX,1:NMAX),VT_OGE(1:NMAX,1:NMAX)

*      TERMINOS TENSORIALES
         DOUBLE PRECISION VT_PI(1:NMAX,1:NMAX)
         DOUBLE PRECISION VT_K(1:NMAX,1:NMAX), VT_ETA(1:NMAX,1:NMAX)
         DOUBLE PRECISION VQQ_T(1:NMAX,1:NMAX)  
         DOUBLE PRECISION V_OGE(1:NMAX,1:NMAX)
         
*      TERMINO ESPIN-ORBITA
         DOUBLE PRECISION Vqq_SO(1:NMAX,1:NMAX)
         DOUBLE PRECISION VSO_CON(1:NMAX,1:NMAX)
         
*      INTERCAMBIO GLUON  VOGE
        
         dimension r(1:nmax), eta(1:nmax)
         dimension autov(1:nmax)
         dimension state(nmax)
         dimension indice(nmax)
         
         JOBVL = 'N'
         JOBVR = 'V'
        
       JJ=2
       L=1
       SPIN=1

      IF ( ( JJ.GE.ABS(L-SPIN) .AND. JJ.LE.L+SPIN )) THEN
         CONTINUE
      ELSE
         STOP
      END IF


        do J = 1, nmax  
           ARG =DBLE(J-1)/DBLE(nmax-1)            
           r(J) = r1*(rmax/r1)**ARG
           eta(J) = 1.d0 /( r(J)**2)
           
*     Normalization coeficient of the gaussian base 
           
           NORM(J) = RNORMA ( L, ETA(J) )

        end do

   

        do i=1,nmax
           do j=1,nmax

              XBETA=eta(i)+eta(j)
*     PRINT*,XBETA
              XFACT=NORM(i)*NORM(j)
              
              S(i,j)=XFACT*INTEG(2*L,XBETA)
              
*     ENERGIA CINETICA
              
              T(i,j)=-1.D0/(2D0*mc2)*XFACT*
     >             (4.d0*eta(j)**2*INTEG(2*l+2,XBETA)-
     >             2.d0*eta(j)*(2*l+3)*INTEG(2*l,XBETA))

c              if(j.eq.1)
c     &             write(*,*) eta(i),NORM(i),T(i,j),S(i,j)
*         POTENCIAL DE YUKAWA


*             TERMINOS CENTRALES
      
       
              if (SPIN.eq.0) then
                 PAULI=-3.D0
              ELSE
                 PAULI=1.D0
              endif
               GM=0.D0

* constantes primeros terminos

*      KM_PI= GCH2*M_PI*(M_PI**2/12*MS*MD)*
*     &    (LAMDA_PI**2/(LAMDA_PI**2-M_PI**2))

*              KM_K= GCH2*M_K*(M_K**2/(12*MS*MD))*
*     &             (LAMDA_k**2/(LAMDA_k**2-M_K**2))

              KM_ETA=GCH2/(12D0*MS*MD)*
     &             (LAMDA_ETA**2/(LAMDA_ETA**2-M_ETA**2))
              
              KM_SIGMA=-GCH2*M_SIGMA*((LAMDA_SIGMA**2)/
     &             (LAMDA_SIGMA**2-M_SIGMA**2))
       
           

*  n=1       


              CALL XINTEG(2*l+1,XBETA,M_ETA,XINTGG)
              CALL XINTEG(2*l+1,XBETA,LAMDA_ETA,XINTGGG)

              CALL XINTEG(2*l+1,XBETA,M_SIGMA,XINTGS)
              CALL XINTEG(2*l+1,XBETA,LAMDA_SIGMA,XINTGS1)
    


              VC_SIGMA(i,j)=KM_SIGMA*((XINTGS/M_SIGMA)-
     &             (XINTGS1/M_SIGMA))

               
              
      

              VC_ETA(i,j)=KM_ETA*(M_ETA**2*XINTGG-
     &            LAMDA_ETA**2*XINTGGG )*PAULI*(-SIN(TETAP))

*      SUMA CONTRIBUCIONES DEL CENTRAL
  
               VQQ_C(i,j)= VC_SIGMA(i,j)+VC_ETA(i,j)

*               PRINT*,VQQ_C(I,J)

 
*   TERMINOS TENSORIALES


              LS_P=0.5D0*(JJ*(JJ+1)-L*(L+1)-SPIN*(SPIN+1))
              LS_N=0
	      
           
 
*              PRINT*,LS_P

              SIJ = 8*((2*L*(L+1)*SPIN*(SPIN+1)-3*LS_P-6*(LS_P**2))/
     &             ((2*L-1)*(2*L+3)))
      
*              PRINT*,SIJ

              CALL XINTEG(2*L,XBETA,M_ETA,XINTGE2)
              CALL XINTEG(2*L-1,XBETA,M_ETA,XINTGE3)
              CALL XINTEG(2*L,XBETA,LAMDA_ETA,XINTGEE2)
              CALL XINTEG(2*L-1,XBETA,LAMDA_ETA,XINTGEE3)
*              PRINT*,XINTGE3

            VT_ETA(I,J)=KM_ETA*(M_ETA**2*XINTGG+3.*M_ETA*XINTGE2 +
     &      3.*XINTGE3 - LAMDA_ETA**2*XINTGGG-
     &      3.*LAMDA_ETA*XINTGEE2 - 3.*XINTGEE3)*SIJ*(-SIN(TETAP))


*              PRINT*,VT_ETA(I,J)
     

             Vqq_T(i,j)= Vt_eta(i,j)




*     INTERCAMBIO GLUON

     

            ALPHAS=ALPHAO/log((MC2**2+MUO**2)/(LAMDAO**2))
*       PRINT*, ALPHAS

            R_MU = (RO*MU)/(2D0*MC2)
            RI_MU=(1D0/R_MU)

            RG_MU=(RG*MU)/(2D0*MC2)
            RGI_MU=(1D0/RG_MU)
*            if(i.eq.1.and.j.eq.1)
*     &           write(*,*) 'alfa ',ALPHAS,R_MU,MC2,RI_MU
             
    
            CALL XINTEG(2*l+1,XBETA,0D0,XINTGG1)
            CALL XINTEG(2*l+1,XBETA,RI_MU,XINTGG2)


            VC_OGE(I,J)=0.25*ALPHAS*LL*(XINTGG1-(PAULI*
     >        XINTGG2)/(6D0*MU*MS*R_MU**2))

            CALL XINTEG(2*L-1,XBETA,0D0,XINTGGT)
            CALL XINTEG(2*L-1,XBETA,RGI_MU,XINTGGT1)
            CALL XINTEG(2*L,XBETA,RGI_MU,XINTGGT2)
            CALL XINTEG(2*L+1,XBETA,RGI_MU,XINTGGT3)
*             PRINT*,XINTGGT

            VT_OGE(I,J)=((-ALPHAS*LL)/(16D0*MU*MS))*(XINTGGT-XINTGGT1
     >           - XINTGGT3*RGI_MU**2/3.D+00 -XINTGGT2*RGI_MU )*SIJ
*           PRINT*,VC_OGE(I,J)
            V_OGE(I,J)=VC_OGE(I,J)+VT_OGE(I,J)

*  TERMINO DE CONFINAMIENTO
 
            CALL XINTEG(2*L+2,XBETA,MUC_CON,XINTGCCON)
       
        
            VC_CON(I,J)=((DLT_CON-AC_CONF)*INTEG(2*l,XBETA)+
     &           AC_CONF*XINTGCCON)*LL



        CALL XINTEG(2*L+1,XBETA,MUC_CON,XINTGCSO)

        VSO_CON(I,J)=-LL*((AC_CONF*MUC_CON)/(4D0*(MU**2)*(MS**2)))*
     >  XINTGCSO*(((MU**2+MS**2)*(1-2*AS_CONF)+4*MU*MS*
     >  (1-AS_CONF))*(LS_P)+(MU**2-MS**2)*(1-AS_CONF)*(LS_N))

c         PRINT*, VSO_CON(I,J)
   
           A(i,j)=(XFACT*(Vqq_C(i,j)+Vqq_T(i,j)+V_OGE(I,J)+VC_CON(I,J))+
     &           T(i,j))/adim

          
         end do
      end do
        

 
c        write(*,*) JOBVL,JOBVR,nr,LDA,LDB,LWORK   
      call DGGEV(JOBVL,JOBVR,NMAX,A,NMAX,S,NMAX,ALPHAR,ALPHAI,
     >     BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)

      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) 'ERROR IN LAPACK'
         WRITE(*,*) 'INFO: ',INFO
      ENDIF

          
      WRITE(*,*) 'FIRST ',NFIRST,' EIGENVALUES FOR L=',l,' : '
      DO i = 1, NFIRST
         indice(i) = -1
         state(i)= 1.D+99
         DO j = 1, nmax
*     este bucle es para que corra el indice del autovalor
            DO K = 1, i-1
               IF ( indice(k) .EQ. j ) GO TO 111
            END DO
            IF(BETA(j) .NE. 0.D+00) TEMP =ALPHAR(j)/BETA(j)+(mu+ms)/adim

          
            IF(TEMP .LT. STATE(i) ) THEN
             
               STATE(I) = TEMP
               indice(i)  = j
           
            END IF
 111        CONTINUE
            
         END DO
         WRITE(*,'(I5,F24.10)') indice(i), STATE(i)

           
         
         
*            IF ((JJ.EQ.0).AND.(SPIN.EQ.0).AND.(L.EQ.0)) THEN
*             WRITE(24,*) INDICE(I),STATE(I)
            
*            ELSE IF ((JJ.EQ.1).AND.(SPIN.EQ.1).AND.(L.EQ.0)) THEN
*              WRITE(25,*) INDICE(I),STATE(I)
           
*            ELSE IF ((L.EQ.1).AND.(JJ.EQ.0).AND.(SPIN.EQ.1)) THEN
*              WRITE(26,*) INDICE(I),STATE(I)
*            ELSE

*              WRITE(27,*) INDICE(I),STATE(I)
*            ENDIF

      END DO
    
       

        
c        cerramos el bucle de los momentos orbitales l=0,1      
c     END DO

      ENDPROGRAM

      SUBROUTINE GammaFun(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function Gamma(x)
C       Input :  x  --- Argument of Gamma(x)
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- Gamma(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END SUBROUTINE



      DOUBLE PRECISION FUNCTION RNORMA( LLL, ETAM )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     
*----------------------------------------------------------------------*
*                                                                      *
*     Author : P. Garcia Ortega                                        *
*                                                                      *
*     Normalization of the meson wave function (in fm^{-(l+1.5)})      *
*                                                                      *
*        Inputs:                                                       *
*     LLL  = Orbital momentum of the meson                             *
*     ETAM = Argument of the exponential (range) [Adim]                *
*                                                                      *
*----------------------------------------------------------------------*   
*
      PARAMETER ( PIPIPI = 3.141592D+00 )
      DOUBLE PRECISION ETAM

      IF( LLL .EQ. 0 ) Sfact =   1.0D+00
      IF( LLL .EQ. 1 ) Sfact =   3.0D+00
      IF( LLL .EQ. 2 ) Sfact =  15.0D+00
      IF( LLL .EQ. 3 ) Sfact = 105.0D+00
      IF( LLL .EQ. 4 ) Sfact = 945.0D+00

      XNUM = 2.D+00**(LLL+2) * (2.D+00 * ETAM)**(LLL+1.5D+00)
      XDEN = SQRT( PIPIPI ) * Sfact
  
      RNORMA = SQRT ( XNUM / XDEN )

      
      RETURN
      END
     
       DOUBLE PRECISION FUNCTION INTEG( ALPHA, BETA )

*     
*----------------------------------------------------------------------*
*                                                                      *
*     Author : P. Garcia Ortega                                        *
*                                                                      *
*     GEM integral                                                     *
*                                                                      *
*----------------------------------------------------------------------*   
*
      INTEGER ALPHA
      DOUBLE PRECISION BETA
      DOUBLE PRECISION ARG, GAMMA

      ARG = DBLE( (ALPHA + 3)/2.D+00 )
      CALL GammaFun(ARG, GAMMA )
      INTEG = GAMMA/(2.D+00*BETA**ARG)
      
      RETURN
      END

 
    
*=== Xinteg ===========================================================*
*
      SUBROUTINE XINTEG (  GMMA, BETA, ARG, XINTG )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

*     
*----------------------------------------------------------------------*
*                                                                      *
*     Author : P. Garcia Ortega                                        *
*                                                                      *
*     Matrix element for terms like r^gamma*e^{-beta*r^2-arg*r}        *
*     for GEM base                                                     *
*                                                                      *
*     Created on 16 March 2022 by    P. Garcia Ortega                  *

*                                                                      *
*        Inputs:                                                       *
*     GMMA  = Power of r                                               *
*     BETA  = Argument of r^2 in the exponential (range) [Adim]        *
*     ARG   = Argument of r in the exponential  [Adim]                 *
*     XINTG =  Int_0^infty r^{gamma}*e^{-beta*r^2-arg*r}dr result      *
*                                                                      *
*----------------------------------------------------------------------*   
*
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      INTEGER GMMA
      DOUBLE PRECISION INTEG,ERF
      DOUBLE PRECISION XINTG, XINTG1,XINTG2,BETA, UUU, RHO, ARG
      
      IF ( ABS(ARG) .LE. 1.E-14 ) THEN
*        If arg is zero, we have the usual Integ(alpha,beta)
*            Aqui pongo GMMA-2 para cuadrar con la
*            definicion de ALPHA en INTEG.
         XINTG = INTEG( GMMA-2, BETA )

      ELSE
*     If arg is not zero, we have a Yukawa-like potential

*        Some definitions (u and rho(u) variables)
         UUU = ARG/(2.D+00*SQRT(BETA))
         RHO = SQRT(PIPIPI)*EXP(UUU**2)*(1.-ERF(UUU))
         IF ( UUU .GT. 5. ) THEN
            RHO = 1./UUU-1./(2.*UUU**3)+3./(4.*UUU**5) +
     &           (-15.)/(8.*UUU**7)+105./(16.*UUU**9)
         END IF
         IF ( GMMA .LT. 0 ) THEN
            XINTG = 0.D+00

         ELSE IF( GMMA .EQ. 0 ) THEN
            XINTG = RHO/(2.D+00*SQRT(BETA))

         ELSE IF ( GMMA .EQ. 1 ) THEN
            XINTG = (1.D+00-UUU*RHO)/(2.D+00*BETA)
         ELSE
            CALL XIMASK ( GMMA-1, BETA, ARG, XINTG1 )
            CALL XIMASK ( GMMA-2, BETA, ARG, XINTG2 )
            XINTG = -UUU*XINTG1/SQRT(BETA)+
     &           (GMMA-1.D+00)/(2.D+00*BETA)*XINTG2
         END IF
         
      END IF

      RETURN
      END SUBROUTINE
*** End of Xinteg subroutine ***


*
*=== Ximask ===========================================================*
*
      SUBROUTINE XIMASK ( GMMA, BETA, ARG, XINTG )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

*     
*----------------------------------------------------------------------*
*     Mask for XINTEG, because Fortran does not allow to call a        *
*     routine inside itself                                            *
*----------------------------------------------------------------------*   
*
      INTEGER GMMA
      DOUBLE PRECISION XINTG, BETA, ARG
      
      CALL XINTEG ( GMMA, BETA, ARG, XINTG )

      RETURN
      END SUBROUTINE
*** End of Xinteg subUTroutine ***
      
        
    

     


      
