    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !****************************************************解析解初始化*********************************************************!
    SUBROUTINE ANALYTICAL_CASE_INITIATION
    USE DECLARATION
    IMPLICIT NONE

    !------U------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            U(I,J)=DSIN( X(I) )*DCOS( YPU(J) )
        END DO
    END DO

    !------施加边界条件------!
    U(1       ,1:JM-1:1)=BCU_AL*U(2       ,1:JM-1:1)+BCU_BL
    U(IM      ,1:JM-1:1)=BCU_AR*U(IM-1    ,1:JM-1:1)+BCU_BR
    !DO J=1,JM-1,1
    !    U(IM,J)=DSIN( X(IM) )*DCOS( YPU(J) )
    !END DO
    U(:       ,0       )=BCU_AB*U(:       ,1       )+BCU_BB
    U(:       ,JM      )=BCU_AT*U(:       ,JM-1    )+BCU_BT

    !------V------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            V(I,J)=-DCOS( XPV(I) )*DSIN( Y(J) )
        END DO
    END DO

    !------施加边界条件------!
    V(1:IM-1:1,1       )=BCV_AB*V(1:IM-1:1,2       )+BCV_BB
    V(1:IM-1:1,JM      )=BCV_AT*V(1:IM-1:1,JM-1    )+BCV_BT
    !DO I=1,IM-1,1
    !    V(I,JM)=-DCOS( XPV(I) )*DSIN( Y(JM) )
    !END DO
    V(0       ,:       )=BCV_AL*V(1       ,:       )+BCV_BL
    V(IM      ,:       )=BCV_AR*V(IM-1    ,:       )+BCV_BR

    !------P------!
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            P(I,J)=0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )!*DEXP(-4.0D0*(-0.5D0*DT)/Re)
        END DO
    END DO

    UN=U
    VN=V
    !!------A-B格式才需要------!
    !UN1=U
    !VN1=V
    UHAT=U
    VHAT=V
    
    
    !!------A-B格式才需要------!
    !!------UN1------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        UN1(I,J)=DSIN( X(I) )*DCOS( YPU(J) )*DEXP(2.0D0*DT/Re)
    !    END DO
    !END DO
    !
    !!------施加边界条件------!
    !UN1(1       ,1:JM-1:1)=BCU_AL*UN1(2       ,1:JM-1:1)+BCU_BL
    !UN1(IM      ,1:JM-1:1)=BCU_AR*UN1(IM-1    ,1:JM-1:1)+BCU_BR
    !!DO J=1,JM-1,1
    !!    UN1(IM,J)=DSIN( X(IM) )*DCOS( YPU(J) )
    !!END DO
    !UN1(:       ,0       )=BCU_AB*UN1(:       ,1       )+BCU_BB
    !UN1(:       ,JM      )=BCU_AT*UN1(:       ,JM-1    )+BCU_BT
    !
    !!------VN1------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        VN1(I,J)=-DCOS( XPV(I) )*DSIN( Y(J) )*DEXP(2.0D0*DT/Re)
    !    END DO
    !END DO
    !
    !!------施加边界条件------!
    !VN1(1:IM-1:1,1       )=BCV_AB*VN1(1:IM-1:1,2       )+BCV_BB
    !VN1(1:IM-1:1,JM      )=BCV_AT*VN1(1:IM-1:1,JM-1    )+BCV_BT
    !!DO I=1,IM-1,1
    !!    VN1(I,JM)=-DCOS( XPV(I) )*DSIN( Y(JM) )
    !!END DO
    !VN1(0       ,:       )=BCV_AL*VN1(1       ,:       )+BCV_BL
    !VN1(IM      ,:       )=BCV_AR*VN1(IM-1    ,:       )+BCV_BR

    RETURN
    END


    !****************************************************生成某时刻解析解流场并计算误差并输出*********************************************************!
    SUBROUTINE ANALYTICAL_CASE_ERROR_OUTPUT
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8),ALLOCATABLE::U_ANA(:,:),V_ANA(:,:),P_ANA(:,:)
    REAL(KIND=8),ALLOCATABLE::U_ERR(:,:),V_ERR(:,:),P_ERR(:,:)

    REAL(KIND=8),ALLOCATABLE::X_OUT(:),Y_OUT(:)
    REAL(KIND=8),ALLOCATABLE::U_OUT(:,:),V_OUT(:,:),UH_OUT(:,:),VH_OUT(:,:),PA_OUT(:,:),PE_OUT(:,:)
    REAL(KIND=8),ALLOCATABLE::PHI_OUT(:,:),DIV_OUT(:,:),DIVH_OUT(:,:)
    REAL(KIND=8),ALLOCATABLE::DIV(:,:),DIVH(:,:)
    REAL(KIND=8)ERRORU_INFINITE,ERRORU_1ST,ERRORU_2ND
    REAL(KIND=8)ERRORV_INFINITE,ERRORV_1ST,ERRORV_2ND
    REAL(KIND=8)ERRORP_INFINITE,ERRORP_1ST,ERRORP_2ND
    REAL(KIND=8)SUM,Q
    
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( U_ANA(IM,0:JM),V_ANA(0:IM,JM),P_ANA(IM-1,JM-1) )
    ALLOCATE( U_ERR(IM,0:JM),V_ERR(0:IM,JM),P_ERR(IM-1,JM-1) )

    !--------------------------------------流场变量-----------------------------------------!
    !------U------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            U_ANA(I,J)=DSIN( X(I) )*DCOS( YPU(J) )*DEXP(-2.0D0*T/Re)
        END DO
    END DO
    !------施加边界条件------!
    U_ANA(1       ,1:JM-1:1)=BCU_AL*U_ANA(2       ,1:JM-1:1)+BCU_BL
    U_ANA(IM      ,1:JM-1:1)=BCU_AR*U_ANA(IM-1    ,1:JM-1:1)+BCU_BR
    U_ANA(:       ,0       )=BCU_AB*U_ANA(:       ,1       )+BCU_BB
    U_ANA(:       ,JM      )=BCU_AT*U_ANA(:       ,JM-1    )+BCU_BT

    U_ERR=U-U_ANA

    !------V------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            V_ANA(I,J)=-DCOS( XPV(I) )*DSIN( Y(J) )*DEXP(-2.0D0*T/Re)
        END DO
    END DO

    !------施加边界条件------!
    V_ANA(1:IM-1:1,1       )=BCV_AB*V_ANA(1:IM-1:1,2       )+BCV_BB
    V_ANA(1:IM-1:1,JM      )=BCV_AT*V_ANA(1:IM-1:1,JM-1    )+BCV_BT
    V_ANA(0       ,:       )=BCV_AL*V_ANA(1       ,:       )+BCV_BL
    V_ANA(IM      ,:       )=BCV_AR*V_ANA(IM-1    ,:       )+BCV_BR

    V_ERR=V-V_ANA

    !------P------!
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            !之前错误的
            !P_ANA(I,J)=0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*(T-0.5D0*DT)/Re)
            !现在正确的
            P_ANA(I,J)=0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*T/Re)
        END DO
    END DO

    P_ERR=P-P_ANA

    ALLOCATE( X_OUT(2*IM+1),Y_OUT(2*JM+1) )
    ALLOCATE( U_OUT(2*IM+1,2*JM+1),V_OUT(2*IM+1,2*JM+1),UH_OUT(2*IM+1,2*JM+1),VH_OUT(2*IM+1,2*JM+1),PA_OUT(2*IM+1,2*JM+1),PE_OUT(2*IM+1,2*JM+1) )
    ALLOCATE( PHI_OUT(2*IM+1,2*JM+1),DIV_OUT(2*IM+1,2*JM+1),DIVH_OUT(2*IM+1,2*JM+1) )
    ALLOCATE( DIV(IM-1,JM-1),DIVH(IM-1,JM-1) )
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

    OPEN(UNIT=10,FILE='ANA&ERR_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","P_ANA","P_ERR","U_ANA","V_ANA","U_ERR","V_ERR"'!
    WRITE(10,*) 'ZONE T="NONAME", I=',2*IM+1,', J=',2*JM+1,', F=POINT'

    !OMEGA=0.0D0
    !OMEGA_OP=0.0D0
    DIV=0.0D0
    DIVH=0.0D0

    X_OUT=0.0D0
    Y_OUT=0.0D0
    U_OUT=0.0D0
    V_OUT=0.0D0
    UH_OUT=0.0D0
    VH_OUT=0.0D0
    PA_OUT=0.0D0
    PE_OUT=0.0D0
    PHI_OUT=0.0D0
    DIV_OUT=0.0D0
    DIVH_OUT=0.0D0


    DO J=1,JM-1,1
        DO I=1,IM-1,1
            DIV(I,J)=( U(I+1,J)-U(I,J) )/( X(I+1)-X(I) )+( V(I,J+1)-V(I,J) )/( Y(J+1)-Y(J) )
        END DO
    END DO
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            DIVH(I,J)=( UHAT(I+1,J)-UHAT(I,J) )/( X(I+1)-X(I) )+( VHAT(I,J+1)-VHAT(I,J) )/( Y(J+1)-Y(J) )
        END DO
    END DO

    DO I=1,IM,1
        X_OUT(2*I-1)=XPV(I-1)
        X_OUT(2*I)=X(I)
    END DO
    X_OUT(2*IM+1)=XPV(IM)
    DO J=1,JM,1
        Y_OUT(2*J-1)=YPU(J-1)
        Y_OUT(2*J)=Y(J)
    END DO
    Y_OUT(2*JM+1)=YPU(JM)

    DO J=1,JM,1
        DO I=1,IM,1
            U_OUT(2*I,2*J-1)=U_ANA(I,J-1)
            V_OUT(2*I-1,2*J)=V_ANA(I-1,J)
            UH_OUT(2*I,2*J-1)=U_ERR(I,J-1)
            VH_OUT(2*I-1,2*J)=V_ERR(I-1,J)
        END DO
    END DO

    DO J=1,JM,1
        V_OUT(2*IM+1,2*J)=V_ANA(IM,J)
        VH_OUT(2*IM+1,2*J)=V_ERR(IM,J)
    END DO
    DO I=1,IM,1
        U_OUT(2*I,2*JM+1)=U_ANA(I,JM)
        UH_OUT(2*I,2*JM+1)=U_ERR(I,JM)
    END DO

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            PA_OUT(2*I+1,2*J+1)=P_ANA(I,J)
            PE_OUT(2*I+1,2*J+1)=P_ERR(I,J)
            DIV_OUT(2*I+1,2*J+1)=DIV(I,J)
            DIVH_OUT(2*I+1,2*J+1)=DIVH(I,J)
        END DO
    END DO

    DO J=0,JM,1
        DO I=0,IM,1
            PHI_OUT(2*I+1,2*J+1)=PHI(I,J)
        END DO
    END DO

    !DO J=1,JM,1
    !    DO I=1,IM,1
    !        OMEGA(I,J)=( V(I,J)-V(I-1,J) )/( XPV(I)-XPV(I-1) )-( U(I,J)-U(I,J-1) )/( YPU(J)-YPU(J-1) )
    !    END DO
    !END DO
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !        OMEGA_OP(I,J)=( OMEGA(I,J)+OMEGA(I+1,J)+OMEGA(I,J+1)+OMEGA(I+1,J+1) )/4.0D0
    !    END DO
    !END DO

    DO J=1,2*JM+1,1
        DO I=1,2*IM+1,1
            WRITE(10,*) X_OUT(I),Y_OUT(J),PA_OUT(I,J),PE_OUT(I,J),U_OUT(I,J),V_OUT(I,J),UH_OUT(I,J),VH_OUT(I,J)
            !WRITE(10,*) X_OUT(I),Y_OUT(J),U_OUT(I,J),V_OUT(I,J),UH_OUT(I,J),VH_OUT(I,J),P_OUT(I,J),PHI_OUT(I,J),DIV_OUT(I,J),DIVH_OUT(I,J)
        END DO
    END DO

    CLOSE(10)
    
    
    !---------------求解误差范数-------------------!
    
    !------U------!
    ERRORU_INFINITE=MAXVAL(DABS(U_ERR))
    SUM=0.0D0
    Q=1.0D0
    DO J=0,JM,1
        DO I=1,IM,1
            SUM=SUM+DABS(U_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORU_1ST=SUM**(1.0D0/Q)
    SUM=0.0D0
    Q=2.0D0
    DO J=0,JM,1
        DO I=1,IM,1
            SUM=SUM+DABS(U_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORU_2ND=SUM**(1.0D0/Q)

    !------V------!
    ERRORV_INFINITE=MAXVAL(DABS(V_ERR))
    SUM=0.0D0
    Q=1.0D0
    DO J=1,JM,1
        DO I=0,IM,1
            SUM=SUM+DABS(V_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORV_1ST=SUM**(1.0D0/Q)
    SUM=0.0D0
    Q=2.0D0
    DO J=1,JM,1
        DO I=0,IM,1
            SUM=SUM+DABS(V_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORV_2ND=SUM**(1.0D0/Q)
    
    !------P------!
    ERRORP_INFINITE=MAXVAL(DABS(P_ERR))
    SUM=0.0D0
    Q=1.0D0
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            SUM=SUM+DABS(P_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORP_1ST=SUM**(1.0D0/Q)
    SUM=0.0D0
    Q=2.0D0
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            SUM=SUM+DABS(P_ERR(I,J))**Q/DBLE(IM-1)**2.0D0
        END DO
    END DO
    ERRORP_2ND=SUM**(1.0D0/Q)


    WRITE(12345677,"( I7,4(1X,E20.12) )")NSTEP,T/PI**2.0D0,ERRORU_INFINITE,ERRORU_1ST,ERRORU_2ND
    WRITE(12345678,"( I7,4(1X,E20.12) )")NSTEP,T/PI**2.0D0,ERRORV_INFINITE,ERRORV_1ST,ERRORV_2ND
    WRITE(12345679,"( I7,4(1X,E20.12) )")NSTEP,T/PI**2.0D0,ERRORP_INFINITE,ERRORP_1ST,ERRORP_2ND


    RETURN
    END
    
    !****************************************************生成某时刻解析解流场并计算误差并输出*********************************************************!
    !SUBROUTINE ANALYTICAL_CASE_ERROR_ANALYSIS
    !USE DECLARATION
    !IMPLICIT NONE
    !
    !REAL(KIND=8),ALLOCATABLE::U_ANA(:,:),V_ANA(:,:),P_ANA(:,:)
    !REAL(KIND=8),ALLOCATABLE::U_ERR(:,:),V_ERR(:,:),P_ERR(:,:)
    !
    !REAL(KIND=8),ALLOCATABLE::X_OUT(:),Y_OUT(:)
    !REAL(KIND=8),ALLOCATABLE::U_OUT(:,:),V_OUT(:,:),UH_OUT(:,:),VH_OUT(:,:),PA_OUT(:,:),PE_OUT(:,:)
    !REAL(KIND=8),ALLOCATABLE::PHI_OUT(:,:),DIV_OUT(:,:),DIVH_OUT(:,:)
    !REAL(KIND=8),ALLOCATABLE::DIV(:,:),DIVH(:,:)
    !CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    !INTEGER REYNOLDS
    !
    !REAL(KIND=8),ALLOCATABLE::RUPN_ANA(:,:),RUCXN_ANA(:,:),RUCYN_ANA(:,:),RUVXN_ANA(:,:),RUVYN_ANA(:,:)
    !REAL(KIND=8),ALLOCATABLE::RUPN_NUM(:,:),RUCXN_NUM(:,:),RUCYN_NUM(:,:),RUVXN_NUM(:,:),RUVYN_NUM(:,:)
    !REAL(KIND=8),ALLOCATABLE::RUPN_ERR(:,:),RUCXN_ERR(:,:),RUCYN_ERR(:,:),RUVXN_ERR(:,:),RUVYN_ERR(:,:)
    !
    !REAL(KIND=8)::UR,UL,UU,UD,VR,VL,VU,VD,D1,D2!各方向合成速度，n时间层
    !
    !ALLOCATE( U_ANA(IM,0:JM),V_ANA(0:IM,JM),P_ANA(IM-1,JM-1) )
    !ALLOCATE( U_ERR(IM,0:JM),V_ERR(0:IM,JM),P_ERR(IM-1,JM-1) )
    !
    !!--------------------------------------流场变量-----------------------------------------!
    !!------U------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        U_ANA(I,J)=DSIN( X(I) )*DCOS( YPU(J) )*DEXP(-2.0D0*T/Re)
    !    END DO
    !END DO
    !!------施加边界条件------!
    !U_ANA(1       ,1:JM-1:1)=BCU_AL*U_ANA(2       ,1:JM-1:1)+BCU_BL
    !U_ANA(IM      ,1:JM-1:1)=BCU_AR*U_ANA(IM-1    ,1:JM-1:1)+BCU_BR
    !!DO J=1,JM-1,1
    !!    U_ANA(IM,J)=DSIN( X(IM) )*DCOS( YPU(J) )
    !!END DO
    !U_ANA(:       ,0       )=BCU_AB*U_ANA(:       ,1       )+BCU_BB
    !U_ANA(:       ,JM      )=BCU_AT*U_ANA(:       ,JM-1    )+BCU_BT
    !
    !U_ERR=U-U_ANA
    !
    !!------V------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        V_ANA(I,J)=-DCOS( XPV(I) )*DSIN( Y(J) )*DEXP(-2.0D0*T/Re)
    !    END DO
    !END DO
    !
    !!------施加边界条件------!
    !V_ANA(1:IM-1:1,1       )=BCV_AB*V_ANA(1:IM-1:1,2       )+BCV_BB
    !V_ANA(1:IM-1:1,JM      )=BCV_AT*V_ANA(1:IM-1:1,JM-1    )+BCV_BT
    !!DO I=1,IM-1,1
    !!    V_ANA(I,JM)=-DCOS( XPV(I) )*DSIN( Y(JM) )
    !!END DO
    !V_ANA(0       ,:       )=BCV_AL*V_ANA(1       ,:       )+BCV_BL
    !V_ANA(IM      ,:       )=BCV_AR*V_ANA(IM-1    ,:       )+BCV_BR
    !
    !V_ERR=V-V_ANA
    !
    !!------P------!
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !        P_ANA(I,J)=0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*(T-0.5D0*DT)/Re)!+
    !    END DO
    !END DO
    !
    !P_ERR=P-P_ANA
    !
    !ALLOCATE( X_OUT(2*IM+1),Y_OUT(2*JM+1) )
    !ALLOCATE( U_OUT(2*IM+1,2*JM+1),V_OUT(2*IM+1,2*JM+1),UH_OUT(2*IM+1,2*JM+1),VH_OUT(2*IM+1,2*JM+1),PA_OUT(2*IM+1,2*JM+1),PE_OUT(2*IM+1,2*JM+1) )
    !ALLOCATE( PHI_OUT(2*IM+1,2*JM+1),DIV_OUT(2*IM+1,2*JM+1),DIVH_OUT(2*IM+1,2*JM+1) )
    !ALLOCATE( DIV(IM-1,JM-1),DIVH(IM-1,JM-1) )
    !WRITE(CHAR_STEP,'(I6.6)') NSTEP
    !REYNOLDS=IDNINT(Re)
    !WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS
    !
    !OPEN(UNIT=10,FILE='ANA&ERR_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
    !WRITE(10,*) 'TITLE="NONAME"'
    !WRITE(10,*) 'VARIABLES="X","Y","P_ANA","P_ERR","U_ANA","V_ANA","U_ERR","V_ERR"'!
    !WRITE(10,*) 'ZONE T="NONAME", I=',2*IM+1,', J=',2*JM+1,', F=POINT'
    !
    !!OMEGA=0.0D0
    !!OMEGA_OP=0.0D0
    !DIV=0.0D0
    !DIVH=0.0D0
    !
    !X_OUT=0.0D0
    !Y_OUT=0.0D0
    !U_OUT=0.0D0
    !V_OUT=0.0D0
    !UH_OUT=0.0D0
    !VH_OUT=0.0D0
    !PA_OUT=0.0D0
    !PE_OUT=0.0D0
    !PHI_OUT=0.0D0
    !DIV_OUT=0.0D0
    !DIVH_OUT=0.0D0
    !
    !
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !        DIV(I,J)=( U(I+1,J)-U(I,J) )/( X(I+1)-X(I) )+( V(I,J+1)-V(I,J) )/( Y(J+1)-Y(J) )
    !    END DO
    !END DO
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !        DIVH(I,J)=( UHAT(I+1,J)-UHAT(I,J) )/( X(I+1)-X(I) )+( VHAT(I,J+1)-VHAT(I,J) )/( Y(J+1)-Y(J) )
    !    END DO
    !END DO
    !
    !DO I=1,IM,1
    !    X_OUT(2*I-1)=XPV(I-1)
    !    X_OUT(2*I)=X(I)
    !END DO
    !X_OUT(2*IM+1)=XPV(IM)
    !DO J=1,JM,1
    !    Y_OUT(2*J-1)=YPU(J-1)
    !    Y_OUT(2*J)=Y(J)
    !END DO
    !Y_OUT(2*JM+1)=YPU(JM)
    !
    !DO J=1,JM,1
    !    DO I=1,IM,1
    !        U_OUT(2*I,2*J-1)=U_ANA(I,J-1)
    !        V_OUT(2*I-1,2*J)=V_ANA(I-1,J)
    !        UH_OUT(2*I,2*J-1)=U_ERR(I,J-1)
    !        VH_OUT(2*I-1,2*J)=V_ERR(I-1,J)
    !    END DO
    !END DO
    !
    !DO J=1,JM,1
    !    V_OUT(2*IM+1,2*J)=V_ANA(IM,J)
    !    VH_OUT(2*IM+1,2*J)=V_ERR(IM,J)
    !END DO
    !DO I=1,IM,1
    !    U_OUT(2*I,2*JM+1)=U_ANA(I,JM)
    !    UH_OUT(2*I,2*JM+1)=U_ERR(I,JM)
    !END DO
    !
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !        PA_OUT(2*I+1,2*J+1)=P_ANA(I,J)
    !        PE_OUT(2*I+1,2*J+1)=P_ERR(I,J)
    !        DIV_OUT(2*I+1,2*J+1)=DIV(I,J)
    !        DIVH_OUT(2*I+1,2*J+1)=DIVH(I,J)
    !    END DO
    !END DO
    !
    !DO J=0,JM,1
    !    DO I=0,IM,1
    !        PHI_OUT(2*I+1,2*J+1)=PHI(I,J)
    !    END DO
    !END DO
    !
    !!DO J=1,JM,1
    !!    DO I=1,IM,1
    !!        OMEGA(I,J)=( V(I,J)-V(I-1,J) )/( XPV(I)-XPV(I-1) )-( U(I,J)-U(I,J-1) )/( YPU(J)-YPU(J-1) )
    !!    END DO
    !!END DO
    !!DO J=1,JM-1,1
    !!    DO I=1,IM-1,1
    !!        OMEGA_OP(I,J)=( OMEGA(I,J)+OMEGA(I+1,J)+OMEGA(I,J+1)+OMEGA(I+1,J+1) )/4.0D0
    !!    END DO
    !!END DO
    !
    !DO J=1,2*JM+1,1
    !    DO I=1,2*IM+1,1
    !        WRITE(10,*) X_OUT(I),Y_OUT(J),PA_OUT(I,J),PE_OUT(I,J),U_OUT(I,J),V_OUT(I,J),UH_OUT(I,J),VH_OUT(I,J)
    !        !WRITE(10,*) X_OUT(I),Y_OUT(J),U_OUT(I,J),V_OUT(I,J),UH_OUT(I,J),VH_OUT(I,J),P_OUT(I,J),PHI_OUT(I,J),DIV_OUT(I,J),DIVH_OUT(I,J)
    !    END DO
    !END DO
    !
    !
    !CLOSE(10)
    !
    !
    !WRITE(12345678,"( I7,3(1X,E20.12) )") NSTEP,MAXVAL(DABS(U_ERR)),MAXVAL(DABS(V_ERR)),MAXVAL(DABS(P_ERR))
    !
    !!--------------------------------------U动量方程组分-----------------------------------------!
    !ALLOCATE( RUPN_ANA(2*IM+1,2*JM+1),RUCXN_ANA(2*IM+1,2*JM+1),RUCYN_ANA(2*IM+1,2*JM+1),RUVXN_ANA(2*IM+1,2*JM+1),RUVYN_ANA(2*IM+1,2*JM+1) )
    !ALLOCATE( RUPN_NUM(2*IM+1,2*JM+1),RUCXN_NUM(2*IM+1,2*JM+1),RUCYN_NUM(2*IM+1,2*JM+1),RUVXN_NUM(2*IM+1,2*JM+1),RUVYN_NUM(2*IM+1,2*JM+1) )
    !ALLOCATE( RUPN_ERR(2*IM+1,2*JM+1),RUCXN_ERR(2*IM+1,2*JM+1),RUCYN_ERR(2*IM+1,2*JM+1),RUVXN_ERR(2*IM+1,2*JM+1),RUVYN_ERR(2*IM+1,2*JM+1) )
    !
    !RUPN_ANA=0.0D0
    !RUCXN_ANA=0.0D0
    !RUCYN_ANA=0.0D0
    !RUVXN_ANA=0.0D0
    !RUVYN_ANA=0.0D0
    !RUPN_NUM=0.0D0
    !RUCXN_NUM=0.0D0
    !RUCYN_NUM=0.0D0
    !RUVXN_NUM=0.0D0
    !RUVYN_NUM=0.0D0
    !RUPN_ERR=0.0D0
    !RUCXN_ERR=0.0D0
    !RUCYN_ERR=0.0D0
    !RUVXN_ERR=0.0D0
    !RUVYN_ERR=0.0D0
    !
    !!--------------------------------------纯解析解-----------------------------------------!
    !!DO J=1,JM-1,1
    !!    DO I=2,IM-1,1
    !!        
    !!    RUPN_ACC(2*I,2*J+1)=
    !!    RUCXN_ACC(2*I,2*J+1)=( D2*D2*UL *UL +(D1*D1-D2*D2)*U_ACC (I,J)*U_ACC (I,J)-D1*D1*UR *UR  )/( D1*D2*(D2-D1) )
    !!    RUCYN_ACC(2*I,2*J+1)=( UU *VU -UD *VD  )/( Y(J+1)-Y(J) )
    !!    RUVXN_ACC(2*I,2*J+1)=2.0D0*( ( XPV(I)-X(I) )*UL-( XPV(I)-XPV(I-1) )*U_ACC(I,J)+( X(I)-XPV(I-1) )*UR )/( ( X(I)-XPV(I-1) )*( XPV(I)-XPV(I-1) )*( XPV(I)-X(I) ) )
    !!    RUVYN_ACC(2*I,2*J+1)=(4.0D0*UD-8.0D0*U_ACC(I,J)+4.0D0*UU)/( Y(J+1)-Y(J) )**2.0D0
    !!
    !!    END IF
    !!
    !!    END DO
    !!END DO
    !
    !!--------------------------------------解析解差分值-----------------------------------------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !
    !    UR=0.5D0*( U_ANA(I+1,J)+U_ANA(I,J) )
    !    UL=0.5D0*( U_ANA(I-1,J)+U_ANA(I,J) )
    !    CALL LINEAR_INTERPOLATION( U_ANA(I,J+1),UU,U_ANA(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !    CALL LINEAR_INTERPOLATION( U_ANA(I,J  ),UD,U_ANA(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !    CALL LINEAR_INTERPOLATION( V_ANA(I,J+1),VU,V_ANA(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !    CALL LINEAR_INTERPOLATION( V_ANA(I,J  ),VD,V_ANA(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !    D1=XPV(I-1)-X(I)
    !    D2=XPV(I)-X(I)
    !    RUPN_ANA(2*I,2*J+1)=( P_ANA(I,J)-P_ANA(I-1,J) )/( XPV(I)-XPV(I-1) )
    !    RUCXN_ANA(2*I,2*J+1)=( D2*D2*UL *UL +(D1*D1-D2*D2)*U_ANA (I,J)*U_ANA (I,J)-D1*D1*UR *UR  )/( D1*D2*(D2-D1) )
    !    RUCYN_ANA(2*I,2*J+1)=( UU *VU -UD *VD  )/( Y(J+1)-Y(J) )
    !    RUVXN_ANA(2*I,2*J+1)=2.0D0*( ( XPV(I)-X(I) )*UL-( XPV(I)-XPV(I-1) )*U_ANA(I,J)+( X(I)-XPV(I-1) )*UR )/( ( X(I)-XPV(I-1) )*( XPV(I)-XPV(I-1) )*( XPV(I)-X(I) ) )
    !    RUVYN_ANA(2*I,2*J+1)=(4.0D0*UD-8.0D0*U_ANA(I,J)+4.0D0*UU)/( Y(J+1)-Y(J) )**2.0D0
    !
    !    IF(I==8 .AND. J==32)THEN
    !        WRITE(12345676,"( I7 )") NSTEP
    !        WRITE(12345676,"( 3(1X,E20.12) )") UL,U_ANA(I,J),UR
    !        WRITE(12345676,"( 3(1X,E20.12) )") XPV(I-1),X(I),XPV(I)
    !        WRITE(12345676,"( 2(1X,E20.12) )") RUCXN_ANA(2*I,2*J+1),RUVXN_ANA(2*I,2*J+1)
    !    END IF
    !    IF(I==16 .AND. J==32)THEN
    !        WRITE(12345675,"( I7 )") NSTEP
    !        WRITE(12345675,"( 3(1X,E20.12) )") UL,U_ANA(I,J),UR
    !        WRITE(12345675,"( 3(1X,E20.12) )") XPV(I-1),X(I),XPV(I)
    !        WRITE(12345675,"( 2(1X,E20.12) )") RUCXN_ANA(2*I,2*J+1),RUVXN_ANA(2*I,2*J+1)
    !    END IF
    !
    !    END DO
    !END DO
    !
    !
    !!--------------------------------------数值解差分值-----------------------------------------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !
    !    UR=0.5D0*( U(I+1,J)+U(I,J) )
    !    UL=0.5D0*( U(I-1,J)+U(I,J) )
    !    CALL LINEAR_INTERPOLATION( U(I,J+1),UU,U(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !    CALL LINEAR_INTERPOLATION( U(I,J  ),UD,U(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !    CALL LINEAR_INTERPOLATION( V(I,J+1),VU,V(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !    CALL LINEAR_INTERPOLATION( V(I,J  ),VD,V(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !    D1=XPV(I-1)-X(I)
    !    D2=XPV(I)-X(I)
    !    RUPN_NUM(2*I,2*J+1)=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
    !    RUCXN_NUM(2*I,2*J+1)=( D2*D2*UL *UL +(D1*D1-D2*D2)*U (I,J)*U (I,J)-D1*D1*UR *UR  )/( D1*D2*(D2-D1) )
    !    RUCYN_NUM(2*I,2*J+1)=( UU *VU -UD *VD  )/( Y(J+1)-Y(J) )
    !    RUVXN_NUM(2*I,2*J+1)=2.0D0*( ( XPV(I)-X(I) )*UL-( XPV(I)-XPV(I-1) )*U(I,J)+( X(I)-XPV(I-1) )*UR )/( ( X(I)-XPV(I-1) )*( XPV(I)-XPV(I-1) )*( XPV(I)-X(I) ) )
    !    RUVYN_NUM(2*I,2*J+1)=(4.0D0*UD-8.0D0*U(I,J)+4.0D0*UU)/( Y(J+1)-Y(J) )**2.0D0
    !
    !    IF(I==8 .AND. J==32)THEN
    !        WRITE(12345676,*) ' '
    !        WRITE(12345676,"( 3(1X,E20.12) )") UL,U(I,J),UR
    !        WRITE(12345676,"( 3(1X,E20.12) )") XPV(I-1),X(I),XPV(I)
    !        WRITE(12345676,"( 2(1X,E20.12) )") RUCXN_NUM(2*I,2*J+1),RUVXN_NUM(2*I,2*J+1)
    !    END IF
    !    IF(I==16 .AND. J==32)THEN
    !        WRITE(12345675,*) ' '
    !        WRITE(12345675,"( 3(1X,E20.12) )") UL,U(I,J),UR
    !        WRITE(12345675,"( 3(1X,E20.12) )") XPV(I-1),X(I),XPV(I)
    !        WRITE(12345675,"( 2(1X,E20.12) )") RUCXN_NUM(2*I,2*J+1),RUVXN_NUM(2*I,2*J+1)
    !    END IF
    !
    !    END DO
    !END DO
    !
    !!--------------------------------------误差值-----------------------------------------!
    !
    !RUPN_ERR= RUPN_NUM- RUPN_ANA
    !RUCXN_ERR=RUCXN_NUM-RUCXN_ANA
    !RUCYN_ERR=RUCYN_NUM-RUCYN_ANA
    !RUVXN_ERR=RUVXN_NUM-RUVXN_ANA
    !RUVYN_ERR=RUVYN_NUM-RUVYN_ANA
    !
    !OPEN(UNIT=10,FILE='ERR_RU下一时间步_RE'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
    !WRITE(10,*) 'TITLE="NONAME"'
    !WRITE(10,'(A,$)') 'VARIABLES="X","Y","RUPN_ERR","RUCXN_ERR","RUCYN_ERR","RUVXN_ERR","RUVYN_ERR","ZERO"'!
    !WRITE(10,*) ' '
    !WRITE(10,*) 'ZONE T="NONAME", I=',2*IM+1,', J=',2*JM+1,', F=POINT'
    !
    !DO J=1,2*JM+1,1
    !    DO I=1,2*IM+1,1
    !        WRITE(10,*) X_OUT(I),Y_OUT(J),RUPN_ERR(I,J),RUCXN_ERR(I,J),RUCYN_ERR(I,J),RUVXN_ERR(I,J)/Re,RUVYN_ERR(I,J)/Re,0.0D0
    !    END DO
    !END DO
    !
    !CLOSE(10)
    !
    !WRITE(12345677,"( I7,5(1X,E20.12) )") NSTEP,MAXVAL(DABS(RUPN_ERR)),MAXVAL(DABS(RUCXN_ERR)),MAXVAL(DABS(RUCYN_ERR)),MAXVAL(DABS(RUVXN_ERR))/Re,MAXVAL(DABS(RUVYN_ERR))/Re
    !
    !I=8
    !J=32
    !WRITE(12345676,*) ' '
    !WRITE(12345676,"( 2(1X,E20.12) )") RUCXN_ERR(2*I,2*J+1),RUVXN_ERR(2*I,2*J+1)
    !
    !I=16
    !J=32
    !WRITE(12345675,*) ' '
    !WRITE(12345675,"( 2(1X,E20.12) )") RUCXN_ERR(2*I,2*J+1),RUVXN_ERR(2*I,2*J+1)
    !
    !
    !RETURN
    !END