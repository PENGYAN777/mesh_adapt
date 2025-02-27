C      ________________________________________________________
C     |                                                        |
C     |  COMPUTE THE DETERMINANT OF A GENERAL FACTORED MATRIX  |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --KFACT'S OUTPUT                         |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         KDET,E--DETERMINANT IS KDET*10.**E (E INTEGER) |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS,ALOG10,DLOG10                |
C     |________________________________________________________|
C
      REAL FUNCTION KDET(E,A)
      REAL A(1),D,F,G
      DOUBLE PRECISION C
      INTEGER E,H,I,J,K,L,M,N
      D = A(1)
      IF ( ABS(D) .EQ. 1236 ) GOTO 10
      WRITE(6,*) 'ERROR: MUST FACTOR WITH KFACT',
     1'BEFORE COMPUTING DETERMINANT'
      STOP
10    E = 0
      IF ( D .LT. 0. ) GOTO 70
      N = A(2)
      IF ( N .EQ. 1 ) GOTO 80
      D = 1.
      F = 2.**64
      G = 1./F
      H = 64
      M = N + 1
      J = 0
      K = 4
      L = 3 - M + M*N
      N = L + M
      DO 40 I = K,L,M
           J = J + 1
           IF ( A(I) .GT. J ) D = -D
           IF ( A(J+N) .GT. J ) D = -D
           D = D*A(I+J)
20         IF ( ABS(D) .LT. F ) GOTO 30
           E = E + H
           D = D*G
           GOTO 20
30         IF ( ABS(D) .GT. G ) GOTO 40
           E = E - H
           D = D*F
           GOTO 30
40    CONTINUE
      D = D*A(L+M)
      IF ( E .NE. 0 ) GOTO 50
      KDET = D
      RETURN
50    IF ( D .EQ. 0. ) GOTO 90
      C = ALOG10(ABS(D)) + E*DLOG10(2.D0)
      E = C
      C = C - E
      IF ( C .LE. 0.D0 ) GOTO 60
      C = C - 1
      E = E + 1
60    F = 10.**C
      IF ( D .LT. 0. ) F = -F
      KDET = F
      RETURN
70    KDET = 0.
      RETURN
80    KDET = A(5)
      RETURN
90    E = 0
      GOTO 70
      END
