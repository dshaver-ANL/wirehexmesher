c-----------------------------------------------------------------------

      subroutine rzero(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END

c-----------------------------------------------------------------------

      subroutine izero(a,n)
      INTEGER A(1)
C
      DO 100 I = 1, N
 100     A(I ) = 0
      return
      END

c-----------------------------------------------------------------------

      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END

c-----------------------------------------------------------------------

      subroutine rone(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      return
      END
c-----------------------------------------------------------------------

      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine chcopy(a,b,n)
      CHARACTER*1 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
C
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END

c-----------------------------------------------------------------------

      subroutine fd_weights_full(xx,x,n,m,c)
c
c     This routine evaluates the derivative based on all points
c     in the stencils.  It is more memory efficient than "fd_weights"
c
c     This set of routines comes from the appendix of
c     A Practical Guide to Pseudospectral Methods, B. Fornberg
c     Cambridge Univ. Press, 1996.   (pff)
c
c     Input parameters:
c       xx -- point at wich the approximations are to be accurate
c       x  -- array of x-ordinates:   x(0:n)
c       n  -- polynomial degree of interpolant (# of points := n+1)
c       m  -- highest order of derivative to be approxxmated at xi
c
c     Output:
c       c  -- set of coefficients c(0:n,0:m).
c             c(j,k) is to be applied at x(j) when
c             the kth derivative is approxxmated by a
c             stencil extending over x(0),x(1),...x(n).
c
c
      real x(0:n),c(0:n,0:m)
c
      c1       = 1.
      c4       = x(0) - xx
c
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
c
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      SUBROUTINE ZWGLL (Z,W,NP)
C
C     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
C     The polynomial degree N=NP-1.
C     Z and W are in single precision, but all the arithmetic
C     operations are done in double precision.
C
      REAL Z(1),W(1)
      ALPHA = 0.
      BETA  = 0.
      CALL ZWGLJ (Z,W,NP,ALPHA,BETA)
      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE ZWGLJ (Z,W,NP,ALPHA,BETA)
C
C     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
C     The polynomial degree N=NP-1.
C     Single precision version.
C
      PARAMETER (NMAX=84)
      PARAMETER (NZD = NMAX)
      REAL*8  ZD(NZD),WD(NZD),ALPHAD,BETAD
      REAL Z(1),W(1),ALPHA,BETA
C
      NPMAX = NZD
      IF (NP.GT.NPMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in ZWGLJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NP=',NP
         STOP
      ENDIF
      ALPHAD = ALPHA
      BETAD  = BETA
      CALL ZWGLJD (ZD,WD,NP,ALPHAD,BETAD)
      DO 100 I=1,NP
         Z(I) = ZD(I)
         W(I) = WD(I)
 100  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE ZWGLJD (Z,W,NP,ALPHA,BETA)
C
C     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
C     The polynomial degree N=NP-1.
C     Double precision version.
C
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  Z(NP),W(NP),ALPHA,BETA
C
      N     = NP-1
      NM1   = N-1
      ONE   = 1.
      TWO   = 2.
C
      IF (NP.LE.1) THEN
       WRITE (6,*) 'ZWGLJD: Minimum number of Gauss-Lobatto points is 2'
       WRITE (6,*) 'ZWGLJD: alpha,beta:',alpha,beta,np
       STOP      
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'ZWGLJD: Alpha and Beta must be greater than -1'
         STOP      
      ENDIF
C
      IF (NM1.GT.0) THEN
         ALPG  = ALPHA+ONE
         BETG  = BETA+ONE
         CALL ZWGJD (Z(2),W(2),NM1,ALPG,BETG)
      ENDIF
      Z(1)  = -ONE
      Z(NP) =  ONE
      DO 100  I=2,NP-1
         W(I) = W(I)/(ONE-Z(I)**2)
 100  CONTINUE
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(1))
      W(1)  = ENDW1 (N,ALPHA,BETA)/(TWO*PD)
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(NP))
      W(NP) = ENDW2 (N,ALPHA,BETA)/(TWO*PD)
C
      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE ZWGJD (Z,W,NP,ALPHA,BETA)
C
C     Generate NP GAUSS JACOBI points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
C     The polynomial degree N=NP-1.
C     Double precision version.
C
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  Z(1),W(1),ALPHA,BETA
C
      N     = NP-1
      DN    = ((N))
      ONE   = 1.
      TWO   = 2.
      APB   = ALPHA+BETA
C
      IF (NP.LE.0) THEN
         WRITE (6,*) 'ZWGJD: Minimum number of Gauss points is 1',np
         STOP      
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'ZWGJD: Alpha and Beta must be greater than -1'
         STOP
      ENDIF
C
      IF (NP.EQ.1) THEN
         Z(1) = (BETA-ALPHA)/(APB+TWO)
         W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO)
     $          * TWO**(APB+ONE)
         RETURN
      ENDIF
C
      CALL JACG (Z,NP,ALPHA,BETA)
C
      NP1   = N+1
      NP2   = N+2
      DNP1  = ((NP1))
      DNP2  = ((NP2))
      FAC1  = DNP1+ALPHA+BETA+ONE
      FAC2  = FAC1+DNP1
      FAC3  = FAC2+ONE
      FNORM = PNORMJ(NP1,ALPHA,BETA)
      RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
      DO 100 I=1,NP
         CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z(I))
         W(I) = -RCOEF/(P*PDM1)
 100  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE JACG (XJAC,NP,ALPHA,BETA)
C
C     Compute NP Gauss points XJAC, which are the zeros of the
C     Jacobi polynomial J(NP) with parameters ALPHA and BETA.
C     ALPHA and BETA determines the specific type of Gauss points.
C     Examples:
C     ALPHA = BETA =  0.0  ->  Legendre points
C     ALPHA = BETA = -0.5  ->  Chebyshev points
C
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  XJAC(1)
      DATA KSTOP /10/
      DATA EPS/1.0e-12/
      N   = NP-1
      one = 1.
      DTH = 4.*ATAN(one)/(2.*((N))+2.)
      DO 40 J=1,NP
         IF (J.EQ.1) THEN
            X = COS((2.*(((J))-1.)+1.)*DTH)
         ELSE
            X1 = COS((2.*(((J))-1.)+1.)*DTH)
            X2 = XLAST
            X  = (X1+X2)/2.
         ENDIF
         DO 30 K=1,KSTOP
            CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X)
            RECSUM = 0.
            JM = J-1
            DO 29 I=1,JM
               RECSUM = RECSUM+1./(X-XJAC(NP-I+1))
 29         CONTINUE
            DELX = -P/(PD-RECSUM*P)
            X    = X+DELX
            IF (ABS(DELX) .LT. EPS) GOTO 31
 30      CONTINUE
 31      CONTINUE
         XJAC(NP-J+1) = X
         XLAST        = X
 40   CONTINUE
      DO 200 I=1,NP
         XMIN = 2.
         DO 100 J=I,NP
            IF (XJAC(J).LT.XMIN) THEN
               XMIN = XJAC(J)
               JMIN = J
            ENDIF
 100     CONTINUE
         IF (JMIN.NE.I) THEN
            SWAP = XJAC(I)
            XJAC(I) = XJAC(JMIN)
            XJAC(JMIN) = SWAP
         ENDIF
 200  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------

      SUBROUTINE JACOBF (POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,
     $                   N,ALP,BET,X)
C
C     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
C     of degree N at X.
C
      IMPLICIT REAL*8  (A-H,O-Z)
      APB  = ALP+BET
      POLY = 1.
      PDER = 0.
      IF (N .EQ. 0) RETURN
      POLYL = POLY
      PDERL = PDER
      POLY  = (ALP-BET+(APB+2.)*X)/2.
      PDER  = (APB+2.)/2.
      IF (N .EQ. 1) RETURN
      DO 20 K=2,N
         DK = ((K))
         A1 = 2.*DK*(DK+APB)*(2.*DK+APB-2.)
         A2 = (2.*DK+APB-1.)*(ALP**2-BET**2)
         B3 = (2.*DK+APB-2.)
         A3 = B3*(B3+1.)*(B3+2.)
         A4 = 2.*(DK+ALP-1.)*(DK+BET-1.)*(2.*DK+APB)
         POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
         PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
         PSAVE  = POLYL
         PDSAVE = PDERL
         POLYL  = POLY
         POLY   = POLYN
         PDERL  = PDER
         PDER   = PDERN
 20   CONTINUE
      POLYM1 = POLYL
      PDERM1 = PDERL
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      RETURN
      END

C--------------------------------------------------------------------

      subroutine zuni(z,np)
c
c     Generate equaly spaced np points on the interval [-1:1]
c
      real z(1)

      dz = 2./(np-1)
      z(1) = -1.
      do i = 2,np-1
         z(i) = z(i-1) + dz
      enddo
      z(np) = 1.

      return
      end

c-----------------------------------------------------------------------

      REAL*8  FUNCTION ENDW1 (N,ALPHA,BETA)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW1 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW1 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2)
     $        * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW1 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW1  = F3
      RETURN
      END

c-----------------------------------------------------------------------

      REAL*8  FUNCTION ENDW2 (N,ALPHA,BETA)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW2 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW2 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2)
     $        * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW2 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW2  = F3
      RETURN
      END

c-----------------------------------------------------------------------

      REAL*8  FUNCTION GAMMAF (X)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  X
      ZERO = 0.0
      HALF = 0.5
      ONE  = 1.0
      TWO  = 2.0
      FOUR = 4.0
      PI   = FOUR*ATAN(ONE)
      GAMMAF = ONE
      IF (X.EQ.-HALF) GAMMAF = -TWO*SQRT(PI)
      IF (X.EQ. HALF) GAMMAF =  SQRT(PI)
      IF (X.EQ. ONE ) GAMMAF =  ONE
      IF (X.EQ. TWO ) GAMMAF =  ONE
      IF (X.EQ. 1.5  ) GAMMAF =  SQRT(PI)/2.
      IF (X.EQ. 2.5) GAMMAF =  1.5*SQRT(PI)/2.
      IF (X.EQ. 3.5) GAMMAF =  0.5*(2.5*(1.5*SQRT(PI)))
      IF (X.EQ. 3. ) GAMMAF =  2.
      IF (X.EQ. 4. ) GAMMAF = 6.
      IF (X.EQ. 5. ) GAMMAF = 24.
      IF (X.EQ. 6. ) GAMMAF = 120.
      RETURN
      END

c-----------------------------------------------------------------------

      REAL*8  FUNCTION PNORMJ (N,ALPHA,BETA)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  ALPHA,BETA
      ONE   = 1.
      TWO   = 2.
      DN    = ((N))
      CONST = ALPHA+BETA+ONE
      IF (N.LE.1) THEN
         PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
         PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
         PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
         RETURN
      ENDIF
      PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
      PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
      PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
      PROD  = PROD*(ONE+BETA)*(TWO+BETA)
      DO 100 I=3,N
         DINDX = ((I))
         FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
         PROD  = PROD*FRAC
 100  CONTINUE
      PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
      RETURN
      END

c-----------------------------------------------------------------------
