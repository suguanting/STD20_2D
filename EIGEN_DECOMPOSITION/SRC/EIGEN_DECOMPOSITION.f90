    !摸索特征分解步骤
    SUBROUTINE SPECTRAL_DECOMPOSITION
    IMPLICIT NONE

    INTEGER(KIND=4)::N, LDA
    INTEGER(KIND=4)::INDEX_1,INDEX_2
    REAL(KIND=8),ALLOCATABLE::A(:,:),A_RESTORE(:,:)

    !----解特征值和特征向量-----!

    !DGEBAL
    INTEGER(KIND=4)::ILO, IHI, INFO
    REAL(KIND=8),ALLOCATABLE::SCALE1(:)
    CHARACTER*1::JOB
    !DGEHRD
    INTEGER(KIND=4)::LWORK
    REAL(KIND=8),ALLOCATABLE::WORK(:),TAU(:)
    INTEGER(KIND=4)::SIZEM(2)

    !DHSEQR
    INTEGER(KIND=4)::LDZ,LDH
    REAL(KIND=8),ALLOCATABLE::Z(:,:),H(:,:),WR(:), WI(:)
    CHARACTER*1::COMPZ

    !DTREVC
    INTEGER(KIND=4)::LDT,LDVL,LDVR,MM,M
    CHARACTER*1::SIDE,HOWMNY
    LOGICAL,ALLOCATABLE::SELECT(:)
    REAL(KIND=8),ALLOCATABLE::T(:,:),VL(:,:),VR(:,:),WORK_DTREVC(:)

    !DGEBAK
    INTEGER(KIND=4)::LDV
    REAL(KIND=8),ALLOCATABLE::V(:,:)

    !DGECON
    CHARACTER*1::NORM
    REAL(KIND=8)::ANORM,RCOND
    REAL(KIND=8),ALLOCATABLE::WORK1(:)
    INTEGER(KIND=4),ALLOCATABLE::IWORK1(:)

    !----组成特征分解矩阵-----!
    REAL(KIND=8),ALLOCATABLE::LAMBDA(:,:),RV(:,:)
    INTEGER(KIND=4),ALLOCATABLE::IPIV(:)

    WRITE(*,*) "************************开始特征分解*************************"
    
    OPEN(UNIT=10,FILE="DIM_MATRIX.DAT")
    READ(10,*) LDA,N
    WRITE(*,*) "矩阵为",LDA,"*",N
    CLOSE(10)
    
    ALLOCATE( A(LDA,N),A_RESTORE(LDA,N) )
    
    OPEN(UNIT=10,FILE="MATRIX_4_EIGENDECOMPOSITION.DAT")
    READ(10,*) A!((A(INDEX_1,INDEX_2),INDEX_2=1,N),INDEX_1=1,LDA)
    CLOSE(10)
    
    !OPEN(UNIT=10,FILE="MATRIX_READ.DAT")
    !WRITE(10,*) A
    !CLOSE(10)
    
    WRITE(*,*) "************************成功读入*************************"

    !LDA=3
    !N=3
    !ALLOCATE( A(LDA,N),A_RESTORE(LDA,N) )
    !
    !
    !!A(1,1)=1.0D0
    !!A(1,2)=2.0D0
    !!A(1,3)=0.0D0
    !!A(2,1)=0.0D0
    !!A(2,2)=3.0D0
    !!A(2,3)=0.0D0
    !!A(3,1)=2.0D0
    !!A(3,2)=-4.0D0
    !!A(3,3)=2.0D0
    !
    !A(1,1)=1.0D0
    !A(1,2)=2.0D0
    !A(1,3)=45.0D0
    !A(2,1)=23.0D0
    !A(2,2)=3.0D0
    !A(2,3)=65.0D0
    !A(3,1)=2.0D0
    !A(3,2)=-4.0D0
    !A(3,3)=2.0D0

    A_RESTORE=A

    !WRITE(*,*)A
    !WRITE(*,*) ((A(INDEX_1,INDEX_2),INDEX_2=1,N),INDEX_1=1,LDA)

    !--------------------------------解特征值和特征向量---------------------------------------------!

    !-----------------------DGEBAL-----------------------------!
    JOB='B'
    ALLOCATE( SCALE1(MAX(1, N)) )
    SCALE1=0.0D0
    
    !call dgebal(job, n, a, lda, ilo, ihi, scale, info)
    CALL DGEBAL(JOB, N, A, LDA, ILO, IHI, SCALE1, INFO)


    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DGEBAL*****************"
    END IF

    !----------------------DGEHRD----------------------------!
    LWORK=MAX(1,N)
    ALLOCATE( WORK(LWORK),TAU(MAX(1,N-1)) )
    WORK=0.0D0
    TAU=0.0D0

    !call dgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
    CALL DGEHRD(N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO)

    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DGEHRD*****************"
    END IF

    !READY FOR DHSEQR
    LDH=MAX(1, N)
    ALLOCATE( H(LDH,MAX(1, N)))
    H=0.0D0
    !H=A
    SIZEM=SHAPE(A)
    DO INDEX_2=1,SIZEM(2)
        DO INDEX_1=1,SIZEM(1)
            IF(INDEX_1-INDEX_2<=1) H(INDEX_1,INDEX_2)=A(INDEX_1,INDEX_2)
        END DO
    END DO


    !---------------------DORGHR------------------------------!
    !call dorghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
    !A 被q覆盖
    CALL DORGHR(N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO)
    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DORGHR*****************"
    END IF

    !---------------------DHSEQR-----------------------------!
    JOB='S'
    COMPZ='V'
    LDZ=MAX(1, N)
    ALLOCATE( Z(LDZ,MAX(1, N)),WR(MAX(1, N)), WI(MAX(1, N)))

    !实际是把q给z
    Z=A
    WR=0.0D0
    WI=0.0D0

    !call dhseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info)
    !H 运行结束后包含T了， Z 包含Q*Z
    CALL DHSEQR(JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, LDZ, WORK, LWORK, INFO)

    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DHSEQR*****************"
    END IF

    !---------------------DTREVC--------------------------!
    SIDE='R'
    HOWMNY='B'
    ALLOCATE( SELECT(MAX (1, N)) )
    LDT=MAX(1, N)
    LDVL=1
    LDVR=N
    MM=N
    M=0
    ALLOCATE( T(LDT,MAX (1, N)),VL(LDVL,MAX(1, MM)),VR(LDVR,MAX(1, MM)) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!就是这个傻逼work数组的大小，坑死老子了，这个库函数怎么这么多坑啊!!!!!!!!!!!!!!!!
    ALLOCATE( WORK_DTREVC(3*N) )
    SELECT=0
    !!!!!!!!!!!!!!!!!!关键啊，不说谁TM能注意到H的输出就是T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    T=H!0.0D0
    VL=0.0D0
    !实际是把QZ给vr
    VR=Z

    !call dtrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info)
    CALL DTREVC(SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK_DTREVC, INFO)
    
    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DTREVC*****************"
    END IF

    !------------------DGEBAK-----------------------!
    JOB='B'
    SIDE='R'
    LDV=LDVR
    ALLOCATE( V(LDV,MAX(1, M)) )
    V=VR

    !call dgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
    CALL DGEBAK(JOB, SIDE, N, ILO, IHI, SCALE1, M, V, LDV, INFO)

    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DGEBAK*****************"
    END IF

    !!WRITE(*,*)V
    !WRITE(*,*) ((V(INDEX_1,INDEX_2),INDEX_2=1,MAX(1, M)),INDEX_1=1,LDV)
    !----------验算特征值求解-----------!
    ALLOCATE( LAMBDA(LDA,N) )
    LAMBDA=0.0D0
    SIZEM=SHAPE(LAMBDA)
    DO INDEX_2=1,SIZEM(2)
        DO INDEX_1=1,SIZEM(1)
            IF(INDEX_1==INDEX_2) LAMBDA(INDEX_1,INDEX_2)=WR(INDEX_1)
        END DO
    END DO

    !WRITE(*,*) ((V(INDEX_1,INDEX_2),INDEX_2=1,MAX(1, M)),INDEX_1=1,LDV)
    !WRITE(*,*) 1
    !WRITE(*,*) ((LAMBDA(INDEX_1,INDEX_2),INDEX_2=1,MAX(1, M)),INDEX_1=1,LDV)
    !WRITE(*,*) 1
    !A=MATMUL( V,LAMBDA )
    !WRITE(*,*) ((A(INDEX_1,INDEX_2),INDEX_2=1,N),INDEX_1=1,LDA)
    !WRITE(*,*) 1
    !
    !A=MATMUL( A_RESTORE,V )
    !WRITE(*,*) ((A(INDEX_1,INDEX_2),INDEX_2=1,N),INDEX_1=1,LDA)
    !WRITE(*,*) 1

    A=MATMUL( A_RESTORE,V )-MATMUL( V,LAMBDA )

    IF(MAXVAL(DABS(A))<1.0D-6)THEN
        WRITE(*,*)"************************特征向量和特征值求解完成*************************"

    ELSE IF(MAXVAL(DABS(A))>=1.0D-6)THEN
        WRITE(*,*)"!!!!!!!!!!!特征向量和特征值求解有误!!!!!!!!!!!!!!"
    END IF



    !--------------------------------组成特征分解矩阵---------------------------------------------!
    ALLOCATE( RV(LDA,N) )

    !----------RV-----------!
    RV=V
    ALLOCATE( IPIV(MAX(1,MIN(LDA, N))) )
    CALL DGETRF( LDA ,N, RV, LDA, IPIV, INFO )
    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DGETRF*****************"
    END IF

    !NORM='1'
    !ALLOCATE( WORK1(4*N),IWORK1(N) )
    !CALL DGECON (NORM, N, RV, LDA, ANORM, RCOND, WORK1, IWORK1, INFO)
    !
    !IF(INFO/=0)THEN
    !    WRITE(*,*)"INFO=",INFO
    !    STOP
    !ELSE
    !    WRITE(*,*) "ANORM=",ANORM
    !    WRITE(*,*) "****************DGECON*****************"
    !END IF


    CALL DGETRI( N, RV, LDA, IPIV, WORK, LWORK, INFO )
    IF(INFO/=0)THEN
        WRITE(*,*)"INFO=",INFO
        STOP
    ELSE
        WRITE(*,*) "****************DGETRI*****************"
    END IF


    A=MATMUL( V,LAMBDA )
    A=MATMUL( A,RV )

    A=A-A_RESTORE
    !WRITE(*,*) ((A(INDEX_1,INDEX_2),INDEX_2=1,N),INDEX_1=1,LDA)

    IF(MAXVAL(DABS(A))<1.0D-6)THEN
        WRITE(*,*)"************************特征分解完成*************************"

    ELSE IF(MAXVAL(DABS(A))>=1.0D-6)THEN
        WRITE(*,*)"!!!!!!!!!!!特征分解有误!!!!!!!!!!!!!!"
    END IF

    A=A_RESTORE

    OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTOR.DAT")
    WRITE(10,*) V
    CLOSE(10)

    OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTORINV.DAT")
    WRITE(10,*) RV
    CLOSE(10)

    OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVALUE.DAT")
    WRITE(10,*) WR
    CLOSE(10)
    
    
    DEALLOCATE(A,A_RESTORE      )
    DEALLOCATE(SCALE1 )
    DEALLOCATE(WORK,TAU            )
    DEALLOCATE(Z,H,WR, WI )
    DEALLOCATE(SELECT          )
    DEALLOCATE(T,VL,VR     )
    DEALLOCATE(V               )
    !DEALLOCATE(WORK1                   )
    !DEALLOCATE(IWORK1                 )
    DEALLOCATE(LAMBDA,RV       )
    DEALLOCATE(IPIV                   )


    RETURN
    END SUBROUTINE