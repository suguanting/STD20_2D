    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***********************************************求解输出二维流场信息********************************************************!
    SUBROUTINE OUTPUT_PLT_1_STAGGERED
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::OMEGA(:,:)
    REAL(KIND=8),ALLOCATABLE::OMEGA_OP(:,:),DIV(:,:),DIVH(:,:),U_OP(:,:),V_OP(:,:)!,U_OPHAT(:,:),V_OPHAT(:,:),U_OPHATB(:,:),V_OPHATB(:,:)
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( OMEGA(IM,JM) )
    ALLOCATE( OMEGA_OP(IM-1,JM-1),DIV(IM-1,JM-1),DIVH(IM-1,JM-1),U_OP(IM-1,JM-1),V_OP(IM-1,JM-1 ) )!,U_OPHAT(IM-1,JM-1),V_OPHAT(IM-1,JM-1 ),U_OPHATB(IM-1,JM-1),V_OPHATB(IM-1,JM-1 )
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

    OPEN(UNIT=10,FILE='2DXYRe'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","U","V","P","PHI","OMEGA","DIV","DIVH"'!
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    OMEGA=0.0D0
    OMEGA_OP=0.0D0
    DIV=0.0D0
    DIVH=0.0D0
    U_OP=0.0D0
    V_OP=0.0D0

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
    DO J=1,JM,1
        DO I=1,IM,1
            OMEGA(I,J)=( V(I,J)-V(I-1,J) )/( XPV(I)-XPV(I-1) )-( U(I,J)-U(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            OMEGA_OP(I,J)=( OMEGA(I,J)+OMEGA(I+1,J)+OMEGA(I,J+1)+OMEGA(I+1,J+1) )/4.0D0
        END DO
    END DO
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            U_OP(I,J)=( U(I,J)+U(I+1,J) )/2.0D0
            V_OP(I,J)=( V(I,J)+V(I,J+1) )/2.0D0
        END DO
    END DO



    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),U_OP(I,J),V_OP(I,J),P(I,J),PHI(I,J),OMEGA_OP(I,J),DIV(I,J),DIVH(I,J)!
        END DO
    END DO

    !输出续算所需速度边界信息
    WRITE(10,*) "U左边界"
    DO J=1,JM-1,1
        WRITE(10,*) U(1,J)
    END DO
    WRITE(10,*) "V下边界"
    DO I=1,IM-1,1
        WRITE(10,*) V(I,1)
    END DO  

    CLOSE(10)

    RETURN
    END

    !***********************************************求解输出二维流场（包含相对速度）信息********************************************************!
    SUBROUTINE OUTPUT_PLT_1_STAGGERED_RELATIVE
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::OMEGA(:,:)
    REAL(KIND=8),ALLOCATABLE::OMEGA_OP(:,:),DIV(:,:),DIVH(:,:),U_OP(:,:),V_OP(:,:)!,U_OPHAT(:,:),V_OPHAT(:,:),U_OPHATB(:,:),V_OPHATB(:,:)
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( OMEGA(IM,JM) )
    ALLOCATE( OMEGA_OP(IM-1,JM-1),DIV(IM-1,JM-1),DIVH(IM-1,JM-1),U_OP(IM-1,JM-1),V_OP(IM-1,JM-1 ) )!,U_OPHAT(IM-1,JM-1),V_OPHAT(IM-1,JM-1 ),U_OPHATB(IM-1,JM-1),V_OPHATB(IM-1,JM-1 )
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

    OPEN(UNIT=10,FILE='2DRELAXYRe'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","U","V","UR","VR","P","PHI","OMEGA","DIV","DIVH"'!
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    OMEGA=0.0D0
    OMEGA_OP=0.0D0
    DIV=0.0D0
    U_OP=0.0D0
    V_OP=0.0D0

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
    DO J=1,JM,1
        DO I=1,IM,1
            OMEGA(I,J)=( V(I,J)-V(I-1,J) )/( XPV(I)-XPV(I-1) )-( U(I,J)-U(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            OMEGA_OP(I,J)=( OMEGA(I,J)+OMEGA(I+1,J)+OMEGA(I,J+1)+OMEGA(I+1,J+1) )/4.0D0
        END DO
    END DO
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            U_OP(I,J)=( U(I,J)+U(I+1,J) )/2.0D0
            V_OP(I,J)=( V(I,J)+V(I,J+1) )/2.0D0
        END DO
    END DO



    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),U_OP(I,J),V_OP(I,J),U_OP(I,J)+1.0D0,V_OP(I,J)+0.0D0,P(I,J),PHI(I,J),OMEGA_OP(I,J),DIV(I,J),DIVH(I,J)
        END DO
    END DO

    !输出续算所需速度边界信息
    WRITE(10,*) "U左边界"
    DO J=1,JM-1,1
        WRITE(10,*) U(1,J)
    END DO
    WRITE(10,*) "V下边界"
    DO I=1,IM-1,1
        WRITE(10,*) V(I,1)
    END DO  

    CLOSE(10)

    RETURN
    END


    !***************************************************屏幕输出************************************************************!
    SUBROUTINE OUTPUT_SCREEN
    USE DECLARATION
    IMPLICIT NONE

    WRITE(*,"( I6,(1X,F10.6),(1X,F10.6),(1X,F12.8) )") NSTEP,T,VELOMAX,ERRORVELOMAX
    WRITE(*,*) "             该时间步计算结束                "
    WRITE(*,*) "============================================="
    !IF( ERRORVELOMAX<=100.0*CRITERIA)THEN
    !    WRITE(*,*) '升阻力系数:'
    !    WRITE(*,"( I6,(1X,F14.10),(1X,F14.10))") NSTEP,CX1,CY1
    !END IF

    RETURN
    END

    !***************************************************输出设置说明文件************************************************************!
    SUBROUTINE OUTPUT_CONFIGURATION_LOG
    USE DECLARATION
    USE CAL_QUADRIC_DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    OPEN(UNIT=10,FILE='设置说明.TXT')
    IF(CONTINUOUS_MESH==1)THEN
        WRITE(10,*) "连续网格"
        WRITE(10,*) "BETAL:",BL
        WRITE(10,*) "BETAR:",BR
        WRITE(10,*) "BETAB:",BB
        WRITE(10,*) "BETAT:",BT
    ELSE IF(CONTINUOUS_MESH==0)THEN
        WRITE(10,*) "分级网格"
    END IF

    IF(CASE_TYPE==1)THEN
        WRITE(10,*) "初场为均匀流场"
    ELSE IF(CASE_TYPE==2)THEN
        WRITE(10,*) "初场为真实流场"
    END IF

    IF(IB_LOCOMOTION==0)THEN
        WRITE(10,*) "间歇性扑翼：静置初场"
    ELSE IF(IB_LOCOMOTION==-3)THEN
        WRITE(10,*) "顶盖驱动流"
    ELSE IF(IB_LOCOMOTION==-2)THEN
        WRITE(10,*) "解析解：T&G涡"
    ELSE IF(IB_LOCOMOTION==-1)THEN
        WRITE(10,*) "静置圆柱绕流"
    ELSE IF(IB_LOCOMOTION==1 .OR. IB_LOCOMOTION==2)THEN
        WRITE(10,*) "1-周期运动（扑翼），2-间歇性飞行（扑翼）"
    ELSE IF(IB_LOCOMOTION==3)THEN
        WRITE(10,*) "圆柱突然启动"
    ELSE IF(IB_LOCOMOTION==4)THEN
        WRITE(10,*) "HEAVING&PLUNGING"
    ELSE IF(IB_LOCOMOTION==5 .OR. IB_LOCOMOTION==7)THEN
        WRITE(10,*) "CAVITY_OSCILLATING/PURE_ROTATING"
    ELSE IF(IB_LOCOMOTION==6)THEN
        WRITE(10,*) "X_OSCILLATING静水振荡圆柱"
    END IF

    IF(IB_LOCOMOTION>=4)THEN
        WRITE(10,*) "KH: ",KH
        WRITE(10,*) "K: ",K
        WRITE(10,*) "H:",H
    END IF

    IF(VISCOUS_TERM_METHOD==1)THEN
        WRITE(10,*) "粘性项离散方式正常"
    ELSE IF(VISCOUS_TERM_METHOD==2)THEN
        WRITE(10,*) "粘性项根据上一时间步决定其离散方式"
    END IF
    
    IF(IB_SHAPE==1)THEN
        WRITE(10,*) "圆"
    ELSE IF(IB_SHAPE==2)THEN
        WRITE(10,*) "椭圆"
    END IF

    IF(BOUNDARY_EXISTENCE_1==1)THEN
        WRITE(10,*) "边界1出现"
        END IF
    IF(BOUNDARY_EXISTENCE_2==1)THEN
        WRITE(10,*) "边界2出现"
    END IF
    
    WRITE(10,*) "Re: ",Re
    WRITE(10,*) "DX3: ",DX3
    WRITE(10,*) "NCYCLE:",NCYCLE
    WRITE(10,*) "DT:",DT

    WRITE(10,*) "IM:",IM
    WRITE(10,*) "JM:",JM


    CLOSE(10)

    RETURN
    END

    !***********************************************输出数值探针数据********************************************************!
    SUBROUTINE OUTPUT_PROBE_STAGGERD_CONTINUOUS
    USE DECLARATION
    IMPLICIT NONE
    INTEGER IP,JP,IU,JU,IV,JV
    REAL(KIND=8)::VALUE_U,VALUE_V,VALUE_P

    !探针1
    IP=PROBE_IPV1
    JP=PROBE_JPU1
    IU=PROBE_IU1
    JU=PROBE_JPU1
    IV=PROBE_IPV1
    JV=PROBE_JV1
    CALL BILINEAR_INTERPOLATION(IU,JU,IV,JV,IP,JP,PROBE_X1,PROBE_Y1,VALUE_U,VALUE_V,VALUE_P)
    WRITE(60,"( I6,3(1X,E15.7E3) )") NSTEP,VALUE_U,VALUE_V,VALUE_P

    !探针2
    IP=PROBE_IPV2
    JP=PROBE_JPU2
    IU=PROBE_IU2
    JU=PROBE_JPU2
    IV=PROBE_IPV2
    JV=PROBE_JV2
    CALL BILINEAR_INTERPOLATION(IU,JU,IV,JV,IP,JP,PROBE_X2,PROBE_Y2,VALUE_U,VALUE_V,VALUE_P)
    WRITE(70,"( I6,3(1X,E15.7E3) )") NSTEP,VALUE_U,VALUE_V,VALUE_P

    !探针3
    IP=PROBE_IPV3
    JP=PROBE_JPU3
    IU=PROBE_IU3
    JU=PROBE_JPU3
    IV=PROBE_IPV3
    JV=PROBE_JV3
    CALL BILINEAR_INTERPOLATION(IU,JU,IV,JV,IP,JP,PROBE_X3,PROBE_Y3,VALUE_U,VALUE_V,VALUE_P)
    WRITE(80,"( I6,3(1X,E15.7E3) )") NSTEP,VALUE_U,VALUE_V,VALUE_P

    !探针4
    IP=PROBE_IPV4
    JP=PROBE_JPU4
    IU=PROBE_IU4
    JU=PROBE_JPU4
    IV=PROBE_IPV4
    JV=PROBE_JV4
    CALL BILINEAR_INTERPOLATION(IU,JU,IV,JV,IP,JP,PROBE_X4,PROBE_Y4,VALUE_U,VALUE_V,VALUE_P)
    WRITE(90,"( I6,3(1X,E15.7E3) )") NSTEP,VALUE_U,VALUE_V,VALUE_P


    RETURN
    END


    !***********************************************输出IB参数分布情况********************************************************!
    SUBROUTINE OUTPUT_IB_STAGGERED
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

    OPEN(UNIT=10,FILE='IB_TYPE_C_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","TUX","TUY","TVX","TVY","IB_CUX","IB_CUY","IB_CVX","IB_CVY"'
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),TYPEUX(I,J),TYPEUY(I,J),TYPEVX(I,J),TYPEVY(I,J),IB_CUX(I,J),IB_CUY(I,J),IB_CVX(I,J),IB_CVY(I,J)
        END DO
    END DO

    CLOSE(10)

    OPEN(UNIT=10,FILE='IB_ITSCT_IPSVL_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","S_UX","S_UY","S_VX","S_VY","V_UX","V_UY","V_VX","V_VY"'
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),IB_ITSCT_UX(I,J),IB_ITSCT_UY(I,J),IB_ITSCT_VX(I,J),IB_ITSCT_VY(I,J),IB_IPSVL_UX(I,J),IB_IPSVL_UY(I,J),IB_IPSVL_VX(I,J),IB_IPSVL_VY(I,J)
        END DO
    END DO

    CLOSE(10)

    OPEN(UNIT=10,FILE='IB_AB_R_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","AUX","BUY","AVX","BVY","RUX","RUY","RVX","RVY"'
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),IB_AUX(I,J),IB_BUY(I,J),IB_AVX(I,J),IB_BVY(I,J),IB_RUX(I,J),IB_RUY(I,J),IB_RVX(I,J),IB_RVY(I,J)
        END DO
    END DO

    CLOSE(10)

    OPEN(UNIT=10,FILE='IBN1_IPSVL_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    WRITE(10,*) 'TITLE="NONAME"'
    WRITE(10,*) 'VARIABLES="X","Y","BN1_VN_UX","BN1_VN_UY","BN1_VN_VX","BN1_VN_VY"'
    WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            WRITE(10,*) XPV(I),YPU(J),IBN1_IPSVL_UXN(I,J),IBN1_IPSVL_VXN(I,J),IBN1_IPSVL_UYN(I,J),IBN1_IPSVL_VYN(I,J)
        END DO
    END DO

    CLOSE(10)

    RETURN
    END