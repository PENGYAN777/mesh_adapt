
C      ________________________________________________________
C     |                                                        |
C     |     INVERT A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --WORK ARRAY WITH AT LEAST 2N ELEMENTS   |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      SUBROUTINE KVERT(V,LV,N,W)
      REAL V(LV,1),W(1),S,T
      INTEGER H,I,J,K,L,M,N,O,P,Q
      IF ( N .EQ. 1 ) GOTO 120
      O = N + 1
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      Q = L
      S = ABS(V(L,L))
      DO 20 H = L,N
           DO 20 I = L,N
                T = ABS(V(I,H))
                IF ( T .LE. S ) GOTO 20
                P = I
                Q = H
                S = T
20    CONTINUE
      W(N+L) = P
      W(O-L) = Q
      DO 30 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,Q)
30         V(I,Q) = T
      S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0. ) GOTO 130
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.
      S = 1./S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0. ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K+N)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
C     --------------------
C     |*** PIVOT ROWS ***|
C     --------------------
      DO 110 J = 1,N
           DO 110 I = 2,N
                P = W(I)
                H = O - I
                T = V(P,J)
                V(P,J) = V(H,J)
                V(H,J) = T
110   CONTINUE
      RETURN
120   IF ( V(1,1) .EQ. 0. ) GOTO 130
      V(1,1) = 1./V(1,1)
      RETURN
130   WRITE(6,*) 'MATRIX HAS NO INVERSE'
      STOP
      END
