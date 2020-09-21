    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************求解椭圆升力推力和效率******************************************************!
    SUBROUTINE CAL_CLCT_ELLIPTIC(ELLIPSE_GEOMETRICAL,ELLIPSE_KINETIC,BOUNDARY_ID)
    USE QUADRIC_PARAMETER
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::ELLIPSE_GEOMETRICAL(36),ELLIPSE_KINETIC(9)
    INTEGER::BOUNDARY_ID
    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    CHARACTER(LEN=2)CHAR_ID
    INTEGER::REYNOLDS

    !运算区域
    REAL(KIND=8)::CEN_ELP(2)
    REAL(KIND=8)::CEN_DEVIATION(2)
    REAL(KIND=8)::LAXIS,SAXIS
    !坐标转换矩阵
    REAL(KIND=8)::MAT_ABS2REL(2,2),MAT_REL2ABS(2,2),MAT_ABS2TRU(2,2),MAT_ABS2FLP(2,2)

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !网格密度
    REAL(KIND=8)::RDA,RDN!圆心角，法向
    !坐标值
    REAL(KIND=8),ALLOCATABLE::XR(:,:),YR(:,:)
    REAL(KIND=8),ALLOCATABLE::XA(:,:),YA(:,:)
    REAL(KIND=8),ALLOCATABLE::NXR(:),NYR(:)
    REAL(KIND=8),ALLOCATABLE::NXA(:),NYA(:)
    REAL(KIND=8),ALLOCATABLE::ARC(:)
    REAL(KIND=8)::XATEMP,YATEMP
    !网格数
    INTEGER::RIM
    INTEGER::RJM
    INTEGER::RI,RJ,RJIM
    !---------插值涉及的绝对网格相关---------!
    !角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !真实世界坐标系（仍是相对系）下翅膀上一点的速度
    REAL(KIND=8)::UT,VT
    !速度压力场
    REAL(KIND=8),ALLOCATABLE::UR(:,:),VR(:,:),PR(:,:)
    REAL(KIND=8),ALLOCATABLE::UA(:,:),VA(:,:)
    REAL(KIND=8),ALLOCATABLE::VELON(:,:),VELOT(:,:)

    !---------气动参数相关(被二分之一rou*u方无量纲化---------!
    REAL(KIND=8),ALLOCATABLE::CU(:),CP(:)
    REAL(KIND=8),ALLOCATABLE::CUX(:),CPX(:)
    REAL(KIND=8),ALLOCATABLE::CUY(:),CPY(:)
    REAL(KIND=8)::CUX_TOTAL,CPX_TOTAL
    REAL(KIND=8)::CUY_TOTAL,CPY_TOTAL
    REAL(KIND=8)::CU_TOTAL,CP_TOTAL
    !------力系数------!
    REAL(KIND=8),ALLOCATABLE::CX(:),CY(:)
    REAL(KIND=8),ALLOCATABLE::CXT(:),CYT(:)
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!计算域绝对坐标系下的力
    REAL(KIND=8)::CXT_TOTAL,CYT_TOTAL!真实世界坐标系下的力
    REAL(KIND=8)::CXF_TOTAL,CYF_TOTAL!扑翼坐标系下的力
    !------功系数------!
    REAL(KIND=8),ALLOCATABLE::CP_PROP(:),CP_ELEV(:),CP_EFFE(:),CP_AERO(:)
    REAL(KIND=8),ALLOCATABLE::EFFICIENCY(:)
    REAL(KIND=8)::CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL,ELEV_EFFICIENCY_TOTAL,PROP_EFFICIENCY_TOTAL,EFFICIENCY_TOTAL


    !------初始化------!
    MAT_ABS2REL(1,1)=ELLIPSE_GEOMETRICAL(1)
    MAT_ABS2REL(1,2)=ELLIPSE_GEOMETRICAL(2)
    MAT_ABS2REL(2,1)=ELLIPSE_GEOMETRICAL(4)
    MAT_ABS2REL(2,2)=ELLIPSE_GEOMETRICAL(5)
    MAT_REL2ABS=TRANSPOSE(MAT_ABS2REL)

    MAT_ABS2TRU(1,1)=DCOS(ABSX_TRUX_ANGLE)
    MAT_ABS2TRU(1,2)=DSIN(ABSX_TRUX_ANGLE)
    MAT_ABS2TRU(2,1)=-DSIN(ABSX_TRUX_ANGLE)
    MAT_ABS2TRU(2,2)=DCOS(ABSX_TRUX_ANGLE)

    MAT_ABS2FLP(1,1)=DCOS(ABSX_UPSTROKE_ANGLE)
    MAT_ABS2FLP(1,2)=DSIN(ABSX_UPSTROKE_ANGLE)
    MAT_ABS2FLP(2,1)=-DSIN(ABSX_UPSTROKE_ANGLE)
    MAT_ABS2FLP(2,2)=DCOS(ABSX_UPSTROKE_ANGLE)

    COX2=ELLIPSE_GEOMETRICAL(7)
    COY2=ELLIPSE_GEOMETRICAL(8)
    COXY=ELLIPSE_GEOMETRICAL(9)
    COX =ELLIPSE_GEOMETRICAL(10)
    COY =ELLIPSE_GEOMETRICAL(11)
    COM =ELLIPSE_GEOMETRICAL(13)

    LAXIS=ELLIPSE_GEOMETRICAL(33)
    SAXIS=ELLIPSE_GEOMETRICAL(34)
    CEN_DEVIATION(1)=ELLIPSE_GEOMETRICAL(35)
    CEN_DEVIATION(2)=ELLIPSE_GEOMETRICAL(36)

    CEN_ELP(1)=ELLIPSE_KINETIC(4)
    CEN_ELP(2)=ELLIPSE_KINETIC(5)

    RDA=DX3/LAXIS

    RIM=IDNINT( 2.0D0*PI/RDA )
    RJM=2
    RJIM=1!与翅膀重合的网格层数

    RDA=2.0D0*PI/DBLE(RIM)
    RDN=DSQRT(2.0D0)*DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响

    ALLOCATE( XR(RIM,RJM),YR(RIM,RJM) )
    ALLOCATE( XA(RIM,RJM),YA(RIM,RJM) )
    ALLOCATE( UR(RIM,RJM),VR(RIM,RJM),PR(RIM,RJM) )
    ALLOCATE( UA(RIM,RJM),VA(RIM,RJM) )
    ALLOCATE( VELON(RIM,RJM),VELOT(RIM,RJM) )
    ALLOCATE( NXR(RIM),NYR(RIM) )
    ALLOCATE( NXA(RIM),NYA(RIM) )
    ALLOCATE( ARC(RIM) )

    !------绘制贴体网格------!
    !翅表面层及法向量
    RJ=RJIM
    DO RI=1,RIM,1
        !从X轴正方向逆时针旋转
        XR(RI,RJ)=LAXIS*DCOS( RDA *(DBLE(RI)-0.5D0) )+CEN_DEVIATION(1)
        YR(RI,RJ)=SAXIS*DSIN( RDA *(DBLE(RI)-0.5D0) )+CEN_DEVIATION(2)
        CALL CRDNT_TRANSFORM(MAT_REL2ABS,XR(RI,RJ),YR(RI,RJ),XATEMP,YATEMP)
        XA(RI,RJ)=XATEMP+CEN_ELP(1)
        YA(RI,RJ)=YATEMP+CEN_ELP(2)
        CALL NORMALVECTOR_2D(ELLIPSE_GEOMETRICAL,XA(RI,RJ),YA(RI,RJ),NXA(RI),NYA(RI))
        CALL CRDNT_TRANSFORM(MAT_ABS2REL,NXA(RI),NYA(RI),NXR(RI),NYR(RI))
    END DO
    !弧长
    DO RI=1,RIM,1
        XATEMP=LAXIS*DABS( DCOS( RDA * DBLE(RI) )-DCOS( RDA * DBLE(RI-1) ) )
        YATEMP=SAXIS*DABS( DSIN( RDA * DBLE(RI) )-DSIN( RDA * DBLE(RI-1) ) )
        ARC(RI)=DSQRT(XATEMP**2.0D0+YATEMP**2.0D0)
    END DO
    !翅外层
    RJ=RJM
    DO RI=1,RIM,1
        XA(RI,RJ)=XA(RI,RJIM)+RDN*NXA(RI)
        YA(RI,RJ)=YA(RI,RJIM)+RDN*NYA(RI)
        XR(RI,RJ)=XR(RI,RJIM)+RDN*NXR(RI)
        YR(RI,RJ)=YR(RI,RJIM)+RDN*NYR(RI)
    END DO

    !------流场及气动参数初始化------!
    UR=0.0D0
    VR=0.0D0
    PR=0.0D0
    UA=0.0D0
    VA=0.0D0

    ALLOCATE( CU(RIM),CP(RIM) )
    ALLOCATE( CUX(RIM),CPX(RIM) )
    ALLOCATE( CUY(RIM),CPY(RIM) )
    ALLOCATE( CX(RIM),CY(RIM) )
    ALLOCATE( CXT(RIM),CYT(RIM) )
    ALLOCATE( CP_PROP(RIM),CP_ELEV(RIM),CP_EFFE(RIM),CP_AERO(RIM) )
    ALLOCATE( EFFICIENCY(RIM) )

    CU=0.0D0
    CP=0.0D0
    CUX=0.0D0
    CPX=0.0D0
    CUY=0.0D0
    CPY=0.0D0
    CUX_TOTAL=0.0D0
    CPX_TOTAL=0.0D0
    CUY_TOTAL=0.0D0
    CPY_TOTAL=0.0D0
    CU_TOTAL =0.0D0
    CP_TOTAL =0.0D0

    CX=0.0D0
    CY=0.0D0
    CXT=0.0D0
    CYT=0.0D0
    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0
    CXT_TOTAL=0.0D0
    CYT_TOTAL=0.0D0
    CXF_TOTAL=0.0D0
    CYF_TOTAL=0.0D0

    CP_PROP=0.0D0
    CP_ELEV=0.0D0
    CP_EFFE=0.0D0
    CP_AERO=0.0D0
    EFFICIENCY=0.0D0
    CP_PROP_TOTAL=0.0D0
    CP_ELEV_TOTAL=0.0D0
    CP_EFFE_TOTAL=0.0D0
    CP_AERO_TOTAL=0.0D0
    ELEV_EFFICIENCY_TOTAL=0.0D0
    PROP_EFFICIENCY_TOTAL=0.0D0
    EFFICIENCY_TOTAL=0.0D0

    !------外层点直接插值并求解切法向速度------!
    DO RJ=1,RJM,1
        IF(RJ/=RJIM)THEN
            DO RI=1,RIM,1
                IF( XA(RI,RJ)>LEIN-CRITERIA .AND. XA(RI,RJ)<RIIN+CRITERIA )THEN
                    IAU=IL+FLOOR( (XA(RI,RJ)-LEIN) / DX3  )
                    IAV=IL+FLOOR( (XA(RI,RJ)-LEIN) / DX3 - 0.5D0  )
                    IAP=IL+FLOOR( (XA(RI,RJ)-LEIN) / DX3 - 0.5D0  )
                ELSE
                    WRITE(*,*)"XA ERROR"
                    STOP
                END IF
                IF( YA(RI,RJ)>BOIN-CRITERIA .AND. YA(RI,RJ)<TOIN+CRITERIA )THEN
                    JAU=JB+FLOOR( (YA(RI,RJ)-BOIN) / DX3 - 0.5D0  )
                    JAV=JB+FLOOR( (YA(RI,RJ)-BOIN) / DX3  )
                    JAP=JB+FLOOR( (YA(RI,RJ)-BOIN) / DX3 - 0.5D0  )
                ELSE
                    WRITE(*,*)"YA ERROR"
                    STOP
                END IF

                CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA(RI,RJ),YA(RI,RJ),UA(RI,RJ),VA(RI,RJ),PR(RI,RJ))

                IF(NYR(RI)>=CRITERIA)THEN
                    VELOT(RI,RJ)= UA(RI,RJ)*NYA(RI)-VA(RI,RJ)*NXA(RI)!(NYA,-NXA)
                    VELON(RI,RJ)= UA(RI,RJ)*NXA(RI)+VA(RI,RJ)*NYA(RI)!(NXA, NYA)
                ELSE
                    VELOT(RI,RJ)=-UA(RI,RJ)*NYA(RI)+VA(RI,RJ)*NXA(RI)!(-NYA,NXA)
                    VELON(RI,RJ)= UA(RI,RJ)*NXA(RI)+VA(RI,RJ)*NYA(RI)!( NXA,NYA)
                END IF

            END DO
        END IF
    END DO

    !------翅膀层赋予动边界运动速度并求解切法向速度------!
    RJ=RJIM
    DO RI=1,RIM,1
        CALL VELOCITY_LB(ELLIPSE_KINETIC,XA(RI,RJ),YA(RI,RJ),UA(RI,RJ),VA(RI,RJ))
        IF(NYR(RI)>=CRITERIA)THEN
            VELOT(RI,RJ)=UA(RI,RJ)*NYA(RI)-VA(RI,RJ)*NXA(RI)!(NYA,-NXA)
            VELON(RI,RJ)=UA(RI,RJ)*NXA(RI)+VA(RI,RJ)*NYA(RI)!(NXA, NYA)
        ELSE
            VELOT(RI,RJ)=-UA(RI,RJ)*NYA(RI)+VA(RI,RJ)*NXA(RI)!(-NYA,NXA)
            VELON(RI,RJ)= UA(RI,RJ)*NXA(RI)+VA(RI,RJ)*NYA(RI)!( NXA,NYA)
        END IF
    END DO

    !------求解相对系下速度------!
    DO RJ=1,RJM,1
        DO RI=1,RIM,1
            CALL CRDNT_TRANSFORM(MAT_ABS2REL,UA(RI,RJ),VA(RI,RJ),UR(RI,RJ),VR(RI,RJ))
        END DO
    END DO


    !METHOD==1
    !------力和效率的计算------!
    RJ=RJIM
    DO RI=1,RIM,1
        CU(RI)=2.0D0/Re*( VELOT(RI,RJ+1)-VELOT(RI,RJ) )/RDN
        CP(RI)=2.0D0*PR(RI,RJ+1)
        !CU沿后掠切向量方向(NYA,-NXA),CP沿法向量负方向(-NXA,-NYA)
        IF(NYR(RI)>=CRITERIA)THEN
            CUX(RI)= CU(RI)*NYA(RI)
            CUY(RI)=-CU(RI)*NXA(RI)
        ELSE!CU沿后掠切向量方向(-NYA,NXA),CP沿法向量负方向(-NXA,-NYA)
            CUX(RI)=-CU(RI)*NYA(RI)
            CUY(RI)= CU(RI)*NXA(RI)
        END IF
        CPX(RI)=-CP(RI)*NXA(RI)
        CPY(RI)=-CP(RI)*NYA(RI)
        !求解压差力和粘性力
        CUX_TOTAL=CUX_TOTAL+CUX(RI)*ARC(RI)/(2.0D0*LAXIS)
        CPX_TOTAL=CPX_TOTAL+CPX(RI)*ARC(RI)/(2.0D0*LAXIS)

        CUY_TOTAL=CUY_TOTAL+CUY(RI)*ARC(RI)/(2.0D0*LAXIS)
        CPY_TOTAL=CPY_TOTAL+CPY(RI)*ARC(RI)/(2.0D0*LAXIS)
        CU_TOTAL=CU_TOTAL+CU(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_TOTAL=CP_TOTAL+CP(RI)*ARC(RI)/(2.0D0*LAXIS)

        !求解计算域绝对坐标系下力分量
        CX(RI)=CUX(RI)+CPX(RI)
        CY(RI)=CUY(RI)+CPY(RI)
        CXC_TOTAL=CXC_TOTAL+CX(RI)*ARC(RI)/(2.0D0*LAXIS)
        CYC_TOTAL=CYC_TOTAL+CY(RI)*ARC(RI)/(2.0D0*LAXIS)

        !求解真实坐标系下力分量与总量
        CALL CRDNT_TRANSFORM(MAT_ABS2TRU,CX(RI),CY(RI),CXT(RI),CYT(RI))
        CXT_TOTAL=CXT_TOTAL+CXT(RI)*ARC(RI)/(2.0D0*LAXIS)
        CYT_TOTAL=CYT_TOTAL+CYT(RI)*ARC(RI)/(2.0D0*LAXIS)

        !求解真实地面坐标系下有用功分量
        CP_PROP(RI)=CXT(RI)*VELO_RATIO*DCOS(TRUX_FLIGHT_ANGLE)
        CP_ELEV(RI)=CYT(RI)*VELO_RATIO*DSIN(TRUX_FLIGHT_ANGLE)
        CP_EFFE(RI)=CP_PROP(RI)+CP_ELEV(RI)
        !求解真实地面坐标系下气动功分量
        CALL CRDNT_TRANSFORM(MAT_ABS2TRU,UA(RI,RJ),VA(RI,RJ),UT,VT)
        CP_AERO(RI)=-CXT(RI)*(VELO_RATIO*DCOS(TRUX_FLIGHT_ANGLE)+UT)-CYT(RI)*(VELO_RATIO*DSIN(TRUX_FLIGHT_ANGLE)+VT)
        !效率分量
        EFFICIENCY(RI)=CP_EFFE(RI)/CP_AERO(RI)
        !求解真实地面坐标系下气动功总量
        CP_PROP_TOTAL=CP_PROP_TOTAL+CP_PROP(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_ELEV_TOTAL=CP_ELEV_TOTAL+CP_ELEV(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_EFFE_TOTAL=CP_EFFE_TOTAL+CP_EFFE(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_AERO_TOTAL=CP_AERO_TOTAL+CP_AERO(RI)*ARC(RI)/(2.0D0*LAXIS)

    END DO

    !求解扑翼坐标系下力总量
    CALL CRDNT_TRANSFORM(MAT_ABS2FLP,CXC_TOTAL,CYC_TOTAL,CXF_TOTAL,CYF_TOTAL)
    !求解真实坐标系下总效率
    PROP_EFFICIENCY_TOTAL=CP_PROP_TOTAL/CP_AERO_TOTAL
    ELEV_EFFICIENCY_TOTAL=CP_ELEV_TOTAL/CP_AERO_TOTAL
    EFFICIENCY_TOTAL=CP_EFFE_TOTAL/CP_AERO_TOTAL

    !*************************输出系数********************************!
    IF(BOUNDARY_ID==1)THEN
        WRITE(300,"( I6,6(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL,CXT_TOTAL,CYT_TOTAL,CXF_TOTAL,CYF_TOTAL
        WRITE(301,"( I6,6(1X,F9.5)                       )")NSTEP,CU_TOTAL,CP_TOTAL,CUX_TOTAL,CUY_TOTAL,CPX_TOTAL,CPY_TOTAL
        WRITE(302,"( I6,4(1X,F9.5)                       )")NSTEP,CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL
        WRITE(303,"( I6,4(1X,F9.5)                       )")NSTEP,ELEV_EFFICIENCY_TOTAL,PROP_EFFICIENCY_TOTAL,EFFICIENCY_TOTAL
    ELSE IF(BOUNDARY_ID==2)THEN
        WRITE(305,"( I6,6(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL,CXT_TOTAL,CYT_TOTAL,CXF_TOTAL,CYF_TOTAL
        WRITE(306,"( I6,6(1X,F9.5)                       )")NSTEP,CU_TOTAL,CP_TOTAL,CUX_TOTAL,CUY_TOTAL,CPX_TOTAL,CPY_TOTAL
        WRITE(307,"( I6,4(1X,F9.5)                       )")NSTEP,CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL
        WRITE(308,"( I6,4(1X,F9.5)                       )")NSTEP,PROP_EFFICIENCY_TOTAL,ELEV_EFFICIENCY_TOTAL,EFFICIENCY_TOTAL
    END IF

    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS
    WRITE(CHAR_ID,'(I2.2)') BOUNDARY_ID
    !*************************输出相对流场********************************!
    IF( MOD(NSTEP,200)==0 )THEN
        OPEN(UNIT=10,FILE='RELA_SECTION'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
        WRITE(10,*) 'TITLE="NONAME"'
        WRITE(10,*) 'VARIABLES="XR","YR","PR","UR","VR","CUX","CUY","CPX","CPY","CX","CY"'
        WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
        RJ=RJIM
        DO RI=1,RIM,1
            WRITE(10,*) XR(RI,RJ),YR(RI,RJ),PR(RI,RJ),UR(RI,RJ),VR(RI,RJ),CUX(RI),CUY(RI),CPX(RI),CPY(RI),CX(RI),CY(RI)
        END DO

        DO RJ=2,RJM,1
            DO RI=1,RIM,1
                WRITE(10,*) XR(RI,RJ),YR(RI,RJ),PR(RI,RJ),UR(RI,RJ),VR(RI,RJ),0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
            END DO
        END DO
        CLOSE(10)
    END IF
    !*************************输出绝对流场********************************!
    IF( MOD(NSTEP,200)==0 )THEN
        OPEN(UNIT=10,FILE='ABSL_SECTION'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
        WRITE(10,*) 'TITLE="NONAME"'
        WRITE(10,*) 'VARIABLES="XA","YA","PR","UA","VA","CUX","CUY","CPX","CPY","CX","CY"'
        WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
        RJ=RJIM
        DO RI=1,RIM,1
            WRITE(10,*) XA(RI,RJ),YA(RI,RJ),PR(RI,RJ),UA(RI,RJ),VA(RI,RJ),CUX(RI),CUY(RI),CPX(RI),CPY(RI),CX(RI),CY(RI)
        END DO

        DO RJ=2,RJM,1
            DO RI=1,RIM,1
                WRITE(10,*) XA(RI,RJ),YA(RI,RJ),PR(RI,RJ),UA(RI,RJ),VA(RI,RJ),0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
            END DO
        END DO
        CLOSE(10)
    END IF
    !!*************************输出压力系数********************************!
    !!IF( MOD(NSTEP,200)==0 )THEN
    !    OPEN(UNIT=10,FILE='P_COEFFICIENT'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.DAT')
    !    RJ=RJIM
    !    DO RI=1,RIM,1
    !        WRITE(10,*) DBLE(RI)/DBLE(RIM),PR(RI,RJ+1)
    !    END DO
    !!END IF

    !------校验法向量等------!
    !IF( MOD(NSTEP,5)==0 )THEN
    !    OPEN(UNIT=10,FILE='SECTIONNXNY'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    !    WRITE(10,*) 'TITLE="NONAME"'
    !    WRITE(10,*) 'VARIABLES="XR","YR","NX","NY","ARC"'
    !    WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
    !    RJ=RJIM
    !    DO RI=1,RIM,1
    !        WRITE(10,*) XR(RI,RJ),YR(RI,RJ),NXR(RI),NYR(RI),ARC(RI)
    !    END DO
    !
    !    DO RJ=2,RJM,1
    !        DO RI=1,RIM,1
    !            WRITE(10,*) XR(RI,RJ),YR(RI,RJ),0.0D0,0.0D0,0.0D0
    !        END DO
    !    END DO
    !    CLOSE(10)
    !END IF
    !------校验切法向速度等------!
    !IF( MOD(NSTEP,5)==0 )THEN
    !    OPEN(UNIT=10,FILE='SECTIONVELO'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
    !    WRITE(10,*) 'TITLE="NONAME"'
    !    WRITE(10,*) 'VARIABLES="XR","YR","VELOT","VELON"'
    !    WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
    !    RJ=RJIM
    !    DO RI=1,RIM,1
    !        WRITE(10,*) XR(RI,RJ),YR(RI,RJ),VELOT(RI,RJ),VELON(RI,RJ)
    !    END DO
    !
    !    DO RJ=2,RJM,1
    !        DO RI=1,RIM,1
    !            WRITE(10,*) XR(RI,RJ),YR(RI,RJ),VELOT(RI,RJ),VELON(RI,RJ)
    !        END DO
    !    END DO
    !    CLOSE(10)
    !END IF

    RETURN
    END SUBROUTINE

    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************坐标转换******************************************************!
    SUBROUTINE CRDNT_TRANSFORM(MAT,CRDNT11,CRDNT12,CRDNT21,CRDNT22)
    IMPLICIT NONE

    REAL(KIND=8)::MAT(2,2)
    REAL(KIND=8)::CRDNT11,CRDNT12,CRDNT21,CRDNT22
    REAL(KIND=8)::COORDINATE1(2),COORDINATE2(2)

    COORDINATE1(1)=CRDNT11
    COORDINATE1(2)=CRDNT12

    COORDINATE2=MATMUL( MAT,COORDINATE1 )

    CRDNT21=COORDINATE2(1)
    CRDNT22=COORDINATE2(2)

    RETURN
    END SUBROUTINE

    !***************************************************三线性插值（若U,V,W,P中某输入值大于五百则******************************************************!
    SUBROUTINE BILINEAR_INTERPOLATION(INDEX_XU,INDEX_YU,INDEX_XV,INDEX_YV,INDEX_XP,INDEX_YP,VALUE_X,VALUE_Y,VALUE_U,VALUE_V,VALUE_P)
    USE DECLARATION
    IMPLICIT NONE

    INTEGER::INDEX_XU,INDEX_YU,INDEX_XV,INDEX_YV,INDEX_XP,INDEX_YP
    REAL(KIND=8)::VALUE_X,VALUE_Y,VALUE_U,VALUE_V,VALUE_P
    REAL(KIND=8)::DX,DY,VX,VY

    !------U------!
    DX=X  (INDEX_XU+1)-X  (INDEX_XU)
    DY=YPU(INDEX_YU+1)-YPU(INDEX_YU)

    VX=VALUE_X-X  (INDEX_XU)
    VY=VALUE_Y-YPU(INDEX_YU)

    VALUE_U=&
        U(INDEX_XU  ,INDEX_YU  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
        U(INDEX_XU+1,INDEX_YU  )*       VX/DX *(1.0D0-VY/DY)+&
        U(INDEX_XU  ,INDEX_YU+1)*(1.0D0-VX/DX)*       VY/DY +&
        U(INDEX_XU+1,INDEX_YU+1)*       VX/DX *       VY/DY

    !------V------!
    DX=XPV(INDEX_XV+1)-XPV(INDEX_XV)
    DY=Y  (INDEX_YV+1)-Y  (INDEX_YV)

    VX=VALUE_X-XPV(INDEX_XV)
    VY=VALUE_Y-Y  (INDEX_YV)

    VALUE_V=&
        V(INDEX_XV  ,INDEX_YV  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
        V(INDEX_XV+1,INDEX_YV  )*       VX/DX *(1.0D0-VY/DY)+&
        V(INDEX_XV  ,INDEX_YV+1)*(1.0D0-VX/DX)*       VY/DY +&
        V(INDEX_XV+1,INDEX_YV+1)*       VX/DX *       VY/DY

    !------P------!
    DX=XPV(INDEX_XP+1)-XPV(INDEX_XP)
    DY=YPU(INDEX_YP+1)-YPU(INDEX_YP)

    VX=VALUE_X-XPV(INDEX_XP)
    VY=VALUE_Y-YPU(INDEX_YP)

    VALUE_P=&
        P(INDEX_XP  ,INDEX_YP  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
        P(INDEX_XP+1,INDEX_YP  )*       VX/DX *(1.0D0-VY/DY)+&
        P(INDEX_XP  ,INDEX_YP+1)*(1.0D0-VX/DX)*       VY/DY +&
        P(INDEX_XP+1,INDEX_YP+1)*       VX/DX *       VY/DY



    RETURN
    END SUBROUTINE


    !***************************************************求解升力推力和效率******************************************************!
    SUBROUTINE CAL_CLCT_CIRCULAR
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !运算区域
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CYLINDER_KINETIC(9)
    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !网格密度
    REAL(KIND=8)::RDC,RDR!周向，径向
    !坐标值
    REAL(KIND=8),ALLOCATABLE::XR(:,:),YR(:,:)
    REAL(KIND=8),ALLOCATABLE::NXR(:),NYR(:)
    REAL(KIND=8)::XA,YA
    !网格数
    INTEGER::RIM
    INTEGER::RJM
    INTEGER::RI,RJ,RJIM
    !---------插值涉及的绝对网格相关---------!
    !角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !速度压力场
    REAL(KIND=8),ALLOCATABLE::UR(:,:),VR(:,:),PR(:,:)
    !REAL(KIND=8)::UA,VA
    REAL(KIND=8),ALLOCATABLE::UA(:,:),VA(:,:)

    !---------力相关(被二分之一rou*u方无量纲化---------!
    REAL(KIND=8),ALLOCATABLE::CU(:),CP(:)
    REAL(KIND=8),ALLOCATABLE::CUX(:),CPX(:)
    REAL(KIND=8),ALLOCATABLE::CUY(:),CPY(:)
    REAL(KIND=8),ALLOCATABLE::CX(:),CY(:)
    !REAL(KIND=8),ALLOCATABLE::CP_PROP(:),CP_ELEV(:),CP_EFFE(:),CP_AERO(:)
    !REAL(KIND=8),ALLOCATABLE::EFFICIENCY(:)
    REAL(KIND=8)::CXS_TOTAL,CYS_TOTAL!,CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL,EFFICIENCY_TOTAL
    REAL(KIND=8)::CUX_TOTAL,CPX_TOTAL
    REAL(KIND=8)::CUY_TOTAL,CPY_TOTAL
    !REAL(KIND=8)::CU_TOTAL,CP_TOTAL
    !INTEGER,ALLOCATABLE::POINT_COUNT_CHORD(:)
    !INTEGER::POINT_COUNT_TOTAL
    !INTEGER::METHOD!粘性力计算方法，1为外层网格插值，2为外两层网格插值

    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------边界1-----------------------------------------!
    IF( BOUNDARY_EXISTENCE_1==1 )THEN
        CYLINDER_KINETIC=QUADRIC_KINETIC_1
        CEN_CRC(1)=QUADRIC_KINETIC_1(4)
        CEN_CRC(2)=QUADRIC_KINETIC_1(5)
        COX2=QUADRIC_GEOMETRICAL_1(7)
        COY2=QUADRIC_GEOMETRICAL_1(8)
        COXY=QUADRIC_GEOMETRICAL_1(9)
        COX=QUADRIC_GEOMETRICAL_1(10)
        COY=QUADRIC_GEOMETRICAL_1(11)
        COM=QUADRIC_GEOMETRICAL_1(13)
        !认为COXY==0，COX2==COY2
        RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

        RJM=3
        RIM=IDNINT( 2.0D0*RADIUS/DX3 )
        RJIM=1!与翅膀重合的网格层数

        RDC=2.0D0*RADIUS*PI/DBLE(RIM)
        RDR=DSQRT(2.0D0)*DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响

        !METHOD=2

        WRITE(*,*) "升阻力系数格点数1：",RIM,RJM

        !NXT=QUADRIC_GEOMETRICAL_1(10)
        !NYT=QUADRIC_GEOMETRICAL_1(11)
        !NZT=QUADRIC_GEOMETRICAL_1(12)
        !NXB=-QUADRIC_GEOMETRICAL_1(10)
        !NYB=-QUADRIC_GEOMETRICAL_1(11)
        !NZB=-QUADRIC_GEOMETRICAL_1(12)

        ALLOCATE( XR(RIM,RJM),YR(RIM,RJM) )
        ALLOCATE( UR(RIM,RJM),VR(RIM,RJM),PR(RIM,RJM) )
        ALLOCATE( UA(RIM,RJM),VA(RIM,RJM) )
        ALLOCATE( NXR(RIM),NYR(RIM) )
        ALLOCATE( CU(RIM),CP(RIM) )
        ALLOCATE( CUX(RIM),CPX(RIM) )
        ALLOCATE( CUY(RIM),CPY(RIM) )
        ALLOCATE( CX(RIM),CY(RIM) )
        !ALLOCATE( CP_PROP(RIM),CP_ELEV(RIM),CP_EFFE(RIM),CP_AERO(RIM) )
        !ALLOCATE( EFFICIENCY(RIM) )

        !------绘制贴体网格------!
        DO RJ=1,RJM,1
            DO RI=1,RIM,1
                !从X轴逆时针旋转
                XR(RI,RJ)=CEN_CRC(1)+(RADIUS+DBLE(RJ-1)*RDR)*DCOS(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                YR(RI,RJ)=CEN_CRC(2)+(RADIUS+DBLE(RJ-1)*RDR)*DSIN(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                NXR(RI)=DCOS(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                NYR(RI)=DSIN(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
            END DO
        END DO

        !------初始化------!
        UR=0.0D0
        VR=0.0D0
        PR=0.0D0
        UA=0.0D0
        VA=0.0D0

        CU=0.0D0
        CP=0.0D0
        CX=0.0D0
        CY=0.0D0
        !CP_PROP=0.0D0
        !CP_ELEV=0.0D0
        !CP_EFFE=0.0D0
        !CP_AERO=0.0D0
        !EFFICIENCY=0.0D0
        CXS_TOTAL=0.0D0
        CYS_TOTAL=0.0D0
        CUX_TOTAL=0.0D0
        CPX_TOTAL=0.0D0
        CUY_TOTAL=0.0D0
        CPY_TOTAL=0.0D0
        !CU_TOTAL =0.0D0
        !CP_TOTAL =0.0D0
        !CP_PROP_TOTAL=0.0D0
        !CP_ELEV_TOTAL=0.0D0
        !CP_EFFE_TOTAL=0.0D0
        !CP_AERO_TOTAL=0.0D0
        !EFFICIENCY_TOTAL=0.0D0


        !------外层点直接插值------!
        DO RJ=1,RJM,1
            IF(RJ/=RJIM)THEN
                DO RI=1,RIM,1
                    XA=XR(RI,RJ)
                    YA=YR(RI,RJ)

                    IF( XA>LEIN-CRITERIA .AND. XA<RIIN+CRITERIA )THEN
                        IAU=IL+FLOOR( (XA-LEIN) / DX3  )
                        IAV=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
                        IAP=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
                    ELSE
                        WRITE(*,*)"XA ERROR"
                        STOP
                    END IF
                    IF( YA>BOIN-CRITERIA .AND. YA<TOIN+CRITERIA )THEN
                        JAU=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
                        JAV=JB+FLOOR( (YA-BOIN) / DX3  )
                        JAP=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
                    ELSE
                        WRITE(*,*)"YA ERROR"
                        STOP
                    END IF

                    CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA,YA,UA(RI,RJ),VA(RI,RJ),PR(RI,RJ))

                    IF(NYR(RI)>=CRITERIA)THEN
                        UR(RI,RJ)=UA(RI,RJ)*NYR(RI)-VA(RI,RJ)*NXR(RI)!(NYR,-NXR)
                        VR(RI,RJ)=UA(RI,RJ)*NXR(RI)+VA(RI,RJ)*NYR(RI)!(NXR, NYR)
                    ELSE
                        UR(RI,RJ)=-UA(RI,RJ)*NYR(RI)+VA(RI,RJ)*NXR(RI)!(-NYR,NXR)
                        VR(RI,RJ)= UA(RI,RJ)*NXR(RI)+VA(RI,RJ)*NYR(RI)!( NXR,NYR)
                    END IF

                END DO
            END IF
        END DO

        !------翅膀层赋予动边界运动速度------!
        RJ=RJIM
        DO RI=1,RIM,1
            XA=XR(RI,RJ)
            YA=YR(RI,RJ)

            CALL VELOCITY_LB(CYLINDER_KINETIC,XA,YA,UA(RI,RJ),VA(RI,RJ))
            IF(NYR(RI)>=CRITERIA)THEN
                UR(RI,RJ)=UA(RI,RJ)*NYR(RI)-VA(RI,RJ)*NXR(RI)!(NYR,-NXR)
                VR(RI,RJ)=UA(RI,RJ)*NXR(RI)+VA(RI,RJ)*NYR(RI)!(NXR, NYR)
            ELSE
                UR(RI,RJ)=-UA(RI,RJ)*NYR(RI)+VA(RI,RJ)*NXR(RI)!(-NYR,NXR)
                VR(RI,RJ)= UA(RI,RJ)*NXR(RI)+VA(RI,RJ)*NYR(RI)!( NXR,NYR)
            END IF
        END DO

        !METHOD==2
        !------力和效率的计算------!
        RJ=RJIM
        DO RI=1,RIM,1
            !IF(XR(RI)<=RXTRL .AND. XR(RI)>=RXLED)THEN
            CU(RI)=2.0D0/Re*( UR(RI,RJ+2)-UR(RI,RJ+1) )/RDR
            CP(RI)=2.0D0*PR(RI,RJ+1)
            !CU沿后掠切向量方向(NYR,-NXR),CP沿法向量负方向(-NXR,-NYR)
            IF(NYR(RI)>=CRITERIA)THEN
                CUX(RI)= CU(RI)*NYR(RI)
                CUY(RI)=-CU(RI)*NXR(RI)
            ELSE!CU沿后掠切向量方向(-NYR,NXR),CP沿法向量负方向(-NXR,-NYR)
                CUX(RI)=-CU(RI)*NYR(RI)
                CUY(RI)= CU(RI)*NXR(RI)
            END IF
            CPX(RI)=-CP(RI)*NXR(RI)
            CPY(RI)=-CP(RI)*NYR(RI)

            CX(RI)=CUX(RI)+CPX(RI)
            CY(RI)=CUY(RI)+CPY(RI)
            !CX(RI)= CU(RI)*NYR(RI)-CP(RI)*NXR(RI)
            !CY(RI)=-CU(RI)*NXR(RI)-CP(RI)*NYR(RI)
            !CX(RI)=-CU(RI)*NYR(RI)-CP(RI)*NXR(RI)
            !CY(RI)= CU(RI)*NXR(RI)-CP(RI)*NYR(RI)

            !CALL CRDNT_TRANSFORM(MAT_REL2ABS,CRX(RI),CRY(RI),CX(RI),CY(RI))
            !CALL CRDNT_TRANSFORM(MAT_ADVANCE,CX(RI),CY(RI),CXS(RI),CYS(RI))
            !CU_TOTAL=CU_TOTAL+CU(RI)*PI/DBLE(RIM)
            !CP_TOTAL=CP_TOTAL+CP(RI)*PI/DBLE(RIM)

            CUX_TOTAL=CUX_TOTAL+CUX(RI)*PI/DBLE(RIM)
            CPX_TOTAL=CPX_TOTAL+CPX(RI)*PI/DBLE(RIM)

            CUY_TOTAL=CUY_TOTAL+CUY(RI)*PI/DBLE(RIM)
            CPY_TOTAL=CPY_TOTAL+CPY(RI)*PI/DBLE(RIM)

            CXS_TOTAL=CXS_TOTAL+CX(RI)*PI/DBLE(RIM)
            CYS_TOTAL=CYS_TOTAL+CY(RI)*PI/DBLE(RIM)

            !EFFICIENCY(RI)=-CXS(RI)*VELO_RATIO/( -CRX(RI,RK)*UR(RI,RJ,RK)-CRY(RI,RK)*VR(RI,RJ+2,RK)-FRZ(RI,RK)*WR(RI,RJ+2,RK) )

            !CP_PROP(RI)=-CXS(RI)*VELO_RATIO*DSIN(TRUX_FLIGHT_ANGLE)
            !CP_ELEV(RI)= CYS(RI)*VELO_RATIO*DCOS(TRUX_FLIGHT_ANGLE)
            !CP_EFFE(RI)=CP_PROP(RI)+CP_ELEV(RI)
            !CP_AERO(RI)=-CRX(RI)*UR(RI,RJ)-CRY(RI)*VR(RI,RJ)
            !EFFICIENCY(RI)=CP_EFFE(RI)/CP_AERO(RI)
            !
            !CP_PROP_TOTAL=CP_PROP_TOTAL+CP_PROP(RI)/DBLE(RIM)
            !CP_ELEV_TOTAL=CP_ELEV_TOTAL+CP_ELEV(RI)/DBLE(RIM)
            !CP_EFFE_TOTAL=CP_EFFE_TOTAL+CP_EFFE(RI)/DBLE(RIM)
            !CP_AERO_TOTAL=CP_AERO_TOTAL+CP_AERO(RI)/DBLE(RIM)

            !END IF

        END DO

        !EFFICIENCY_TOTAL=CP_EFFE_TOTAL/CP_AERO_TOTAL

        !*************************输出系数********************************!
        WRITE(300,"( I6,6(1X,F9.5)                       )")NSTEP,CXS_TOTAL,CYS_TOTAL!,CXT_TOTAL,CYT_TOTAL,CXF_TOTAL,CYF_TOTAL!,CU_TOTAL,CP_TOTAL
        !WRITE(301,"( I6,4(1X,F9.5)                       )")NSTEP,CUX_TOTAL,CUY_TOTAL,CPX_TOTAL,CPY_TOTAL
        !WRITE(302,"( I6,4(1X,F9.5)                       )")NSTEP,CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL
        !WRITE(303,"( I6,4(1X,F9.5)                       )")NSTEP,ELEV_EFFICIENCY_TOTAL,PROP_EFFICIENCY_TOTAL,EFFICIENCY_TOTAL

        WRITE(CHAR_STEP,'(I6.6)') NSTEP
        REYNOLDS=IDNINT(Re)
        WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

        !!OPEN(UNIT=10,FILE='CLCTE_N'//TRIM(CHAR_STEP)//'.DAT')
        !!DO RI=1,RIM,1
        !!    WRITE(10,"( (F7.4),2(1X,F9.5),4(1X,F11.5),(1X,F9.5) )") XR(RI)-R_LE,CXS(RI),CYS(RI),CP_PROP(RI),CP_ELEV(RI),CP_EFFE(RI),CP_AERO(RI),EFFICIENCY(RI)
        !!END DO
        !!CLOSE(10)
        !WRITE(310,"( I6,120(1X,F9.5) )") NSTEP,-CXS
        !WRITE(320,"( I6,120(1X,F9.5) )") NSTEP,CYS
        !WRITE(330,"( I6,120(1X,F9.5) )") NSTEP,CP_PROP
        !WRITE(340,"( I6,120(1X,F9.5) )") NSTEP,CP_ELEV
        !WRITE(350,"( I6,120(1X,F9.5) )") NSTEP,CP_EFFE
        !WRITE(360,"( I6,120(1X,F9.5) )") NSTEP,CP_AERO
        !WRITE(370,"( I6,120(1X,F9.5) )") NSTEP,EFFICIENCY

        !*************************输出流场********************************!
        IF( MOD(NSTEP,1000)==0 )THEN
            OPEN(UNIT=10,FILE='SECTION2_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
            WRITE(10,*) 'TITLE="NONAME"'
            WRITE(10,*) 'VARIABLES="XR","YR","PR","UA","VA","CUX","CUY","CPX","CPY","CX","CY"'
            WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
            RJ=RJIM
            DO RI=1,RIM,1
                WRITE(10,*) XR(RI,RJ),YR(RI,RJ),PR(RI,RJ),UA(RI,RJ),VA(RI,RJ),CUX(RI),CUY(RI),CPX(RI),CPY(RI),CX(RI),CY(RI)
                !!!!!!!!注意此处UA UR的区别
            END DO

            DO RJ=2,RJM,1
                DO RI=1,RIM,1
                    WRITE(10,*) XR(RI,RJ),YR(RI,RJ),PR(RI,RJ),UA(RI,RJ),VA(RI,RJ),0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
                END DO
            END DO
            CLOSE(10)
        END IF

    END IF

    RETURN
    END SUBROUTINE

    !*********************求解升力推力（表面力积分法4）******************!
    SUBROUTINE CAL_CLCT_SRFC_INTGRTN_SOME(BOUNDARY_GEOMETRICAL,&
        BOUNDARY_KINETIC,BOUNDARY_ID)
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::BOUNDARY_GEOMETRICAL(36),BOUNDARY_KINETIC(9)
    INTEGER::BOUNDARY_ID
    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    CHARACTER(LEN=2)CHAR_ID
    INTEGER::REYNOLDS

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM,COZ2,COXZ,COYZ,COZ
    !网格密度
    REAL(KIND=8),ALLOCATABLE::DS(:)
    !网格数
    INTEGER::RSM
    !插值涉及的绝对网格角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !坐标值
    REAL(KIND=8),ALLOCATABLE::X_SRFC(:),Y_SRFC(:),Z_SRFC(:)
    REAL(KIND=8)::XA,YA
    !法向量
    REAL(KIND=8),ALLOCATABLE::N1X(:),N1Y(:),N1Z(:)
    !偏导系数
    REAL(KIND=8),ALLOCATABLE::AXX(:),AXY(:),AXZ(:)
    REAL(KIND=8),ALLOCATABLE::AYX(:),AYY(:),AYZ(:)
    REAL(KIND=8),ALLOCATABLE::AZX(:),AZY(:),AZZ(:)
    !固壁处物理量
    REAL(KIND=8),ALLOCATABLE::U_SRFC(:),V_SRFC(:),W_SRFC(:),P_SRFC(:)
    !包括半径
    REAL(KIND=8)::SEARCH_RADIUS_U,SEARCH_RADIUS_V
    !包括点数
    INTEGER::NUM_ELGBL_U,NUM_ELGBL_V,COUNT_ELGBL

    !一些中间变量
    REAL(KIND=8)::ELLIPTIC_H,ELLIPTIC_PERIMETER!椭圆周长
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_DA(:)!椭圆各段弧对应圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AS(:)!椭圆各段弧起始处的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AM(:)!椭圆各段弧中间的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_RADIUS(:)!椭圆弧长对应半径
    REAL(KIND=8)::AVERAGE_RADIUS!椭圆平均半径
    LOGICAL::ELIGIBILITY_U(-2:2,-2:2),ELIGIBILITY_V(-2:2,-2:2)!可视流体点
    REAL(KIND=8)::DISTANCE_U(-2:2,-2:2),DISTANCE_V(-2:2,-2:2)!各点距离
    INTEGER,ALLOCATABLE::INDEX_U(:,:),INDEX_V(:,:)!求解偏导的几个点的脚标
    REAL(KIND=8),ALLOCATABLE::TEMP_1(:,:),TEMP_2(:,:)!偏导系数求解临时矩阵
    REAL(KIND=8)::TEMP_SOLUTION(3,1)!求解结果临时矩阵
    LOGICAL :: OK_FLAG
    REAL(KIND=8)::DEVIATION_X,DEVIATION_Y!坐标偏移
    REAL(KIND=8)::MAX_X,MAX_Y!坐标偏移
    REAL(KIND=8)::MIN_X,MIN_Y!坐标偏移
    REAL(KIND=8)::XATEMP,YATEMP

    !对象特殊性
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::CEN_ELP(2)
    REAL(KIND=8)::CEN_DEVIATION(2)
    REAL(KIND=8)::LAXIS,SAXIS
    !坐标转换矩阵
    REAL(KIND=8)::MAT_ABS2REL(2,2),MAT_REL2ABS(2,2)

    !力系数
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL,CZC_TOTAL!计算域绝对坐标系下的力
    REAL(KIND=8)::CXP,CYP,CZP!压力
    REAL(KIND=8)::CXV,CYV,CZV!粘性

    !------初始化------!
    COX2=BOUNDARY_GEOMETRICAL(7)
    COY2=BOUNDARY_GEOMETRICAL(8)
    COXY=BOUNDARY_GEOMETRICAL(9)
    COX =BOUNDARY_GEOMETRICAL(10)
    COY =BOUNDARY_GEOMETRICAL(11)
    COM =BOUNDARY_GEOMETRICAL(13)
    COZ2=0.0D0
    COXZ=0.0D0
    COYZ=0.0D0
    COZ =0.0D0

    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0

    !------离散界面------!
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!
    IF(IB_SHAPE==1)THEN!1圆
        RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)&
            **0.5D0
        RSM=IDNINT( 2.0D0*RADIUS*PI/DX3 )

        ALLOCATE( DS(RSM) )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        DS=2.0D0*RADIUS*PI/DBLE(RSM)

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        CEN_CRC(1)=BOUNDARY_KINETIC(4)
        CEN_CRC(2)=BOUNDARY_KINETIC(5)

        DO I=1,RSM,1
            X_SRFC(I)=CEN_CRC(1)+RADIUS*DCOS(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
            Y_SRFC(I)=CEN_CRC(2)+RADIUS*DSIN(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
        END DO

    ELSE IF(IB_SHAPE==2)THEN!2椭圆
        LAXIS=BOUNDARY_GEOMETRICAL(33)
        SAXIS=BOUNDARY_GEOMETRICAL(34)

        CEN_DEVIATION(1)=BOUNDARY_GEOMETRICAL(35)
        CEN_DEVIATION(2)=BOUNDARY_GEOMETRICAL(36)

        CEN_ELP(1)=BOUNDARY_KINETIC(4)
        CEN_ELP(2)=BOUNDARY_KINETIC(5)

        MAT_ABS2REL(1,1)=BOUNDARY_GEOMETRICAL(1)
        MAT_ABS2REL(1,2)=BOUNDARY_GEOMETRICAL(2)
        MAT_ABS2REL(2,1)=BOUNDARY_GEOMETRICAL(4)
        MAT_ABS2REL(2,2)=BOUNDARY_GEOMETRICAL(5)
        MAT_REL2ABS=TRANSPOSE(MAT_ABS2REL)

        !估求椭圆周长
        ELLIPTIC_H=((LAXIS-SAXIS)/(LAXIS+SAXIS))**2.0D0
        ELLIPTIC_H=3.0D0*ELLIPTIC_H/(10.0D0+(4.0D0-3.0D0*ELLIPTIC_H)**0.5D0)
        ELLIPTIC_PERIMETER=PI*(LAXIS+SAXIS)*(1.0D0+ELLIPTIC_H)
        !确定网格数
        RSM=IDNINT( ELLIPTIC_PERIMETER/DX3 )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        ALLOCATE( DS(RSM) )

        ALLOCATE( ELLIPTIC_DA(RSM) )
        ALLOCATE( ELLIPTIC_AS(RSM),ELLIPTIC_AM(RSM),ELLIPTIC_RADIUS(RSM) )

        DS=0.0D0
        ELLIPTIC_DA=0.0D0
        ELLIPTIC_AM=0.0D0
        ELLIPTIC_AS=0.0D0

        !确定平均半径周长/2PI
        AVERAGE_RADIUS=ELLIPTIC_PERIMETER/PI/2.0D0

        !初始化
        !ELLIPTIC_A自x轴正方向为零，顺时针方向为正方向
        DO I=1,RSM,1
            ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)
            ELLIPTIC_AS(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.5D0)
            ELLIPTIC_AM(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.0D0)
            ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
        END DO

        !迭代确定各段弧
        DO N=1,10,1

            !修正圆心角
            DO I=1,RSM,1
                ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)*AVERAGE_RADIUS/ELLIPTIC_RADIUS(I)
            END DO
            ELLIPTIC_DA=ELLIPTIC_DA*2.0D0*PI/SUM(ELLIPTIC_DA)
            !重新获得ELLIPTIC_AS和ELLIPTIC_AM
            ELLIPTIC_AS(1)=-ELLIPTIC_DA(1)/2.0D0
            ELLIPTIC_AM(1)=0.0D0
            DO I=2,RSM,1
                ELLIPTIC_AS(I)=ELLIPTIC_AS(I-1)+ELLIPTIC_DA(I-1)
                ELLIPTIC_AM(I)=ELLIPTIC_AS(I)+ELLIPTIC_DA(I)/2.0D0
            END DO
            !先求起始点坐标
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AS(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AS(I)))**2.0D0 )**0.5D0
            END DO
            DO I=1,RSM,1
                X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS( ELLIPTIC_AS(I) )
                Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN( ELLIPTIC_AS(I) )
            END DO
            !求解弧长
            DO I=1,RSM-1,1
                DS(I)=DSQRT(&
                    (X_SRFC(I+1)-X_SRFC(I))**2.0D0+&
                    (Y_SRFC(I+1)-Y_SRFC(I))**2.0D0)
            END DO
            DS(RSM)=DSQRT(&
                (X_SRFC(1)-X_SRFC(RSM))**2.0D0+&
                (Y_SRFC(1)-Y_SRFC(RSM))**2.0D0)
            !检查弧长差异,可以的话就跳出
            IF(MAXVAL(DS)/MINVAL(DS)<1.1D0)EXIT

            !不行的话再重新修正
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
            END DO

        END DO
        !检查弧长差异,太过了就报错
        IF(MAXVAL(DS)>=2.0D0*DX3)THEN
            WRITE(*,*)'too large max(ds)/dx=',MAXVAL(DS)/DX3
            STOP
        END IF
        !求各弧中心点坐标
        DO I=1,RSM,1
            X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(1)
            Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(2)
            CALL CRDNT_TRANSFORM(MAT_REL2ABS,X_SRFC(I),Y_SRFC(I),&
                XATEMP,YATEMP)
            X_SRFC(I)=XATEMP+CEN_ELP(1)
            Y_SRFC(I)=YATEMP+CEN_ELP(2)
        END DO

    END IF
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!

    !------求解内法向量------!
    DO I=1,RSM,1
        CALL NORMALVECTOR_2D(BOUNDARY_GEOMETRICAL,X_SRFC(I),Y_SRFC(I),&
            N1X(I),N1Y(I))
        N1X(I)=-N1X(I)
        N1Y(I)=-N1Y(I)
        N1Z(I)=0.0D0
    END DO
    !------求解压力和偏导系数------!
    DO N=1,RSM,1

        !------求解点所处位置及压力------!
        XA=X_SRFC(N)
        YA=Y_SRFC(N)

        IF( XA>LEIN-CRITERIA .AND. XA<RIIN+CRITERIA )THEN
            IAU=IL+FLOOR( (XA-LEIN) / DX3  )
            IAV=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
            IAP=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"XA OUT OF INNER REGION"
            STOP
        END IF
        IF( YA>BOIN-CRITERIA .AND. YA<TOIN+CRITERIA )THEN
            JAU=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
            JAV=JB+FLOOR( (YA-BOIN) / DX3  )
            JAP=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"YA OUT OF INNER REGION"
            STOP
        END IF

        CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA,YA,&
            U_SRFC(N),V_SRFC(N),P_SRFC(N))

        !------求解点刚体速度------!
        CALL VELOCITY_LB(BOUNDARY_KINETIC,XA,YA,U_SRFC(N),V_SRFC(N))

        !------获取附近符合要求的点------!
        CALL DETERMINE_IF_ELIGIBLE_U(IAU,JAU,ELIGIBILITY_U)
        CALL DETERMINE_IF_ELIGIBLE_V(IAV,JAV,ELIGIBILITY_V)

        !------求解距离------!
        DISTANCE_U=500.0D0
        DISTANCE_V=500.0D0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    DISTANCE_U(I,J)=((X(IAU+I)-XA)**2.0D0&
                        +(YPU(JAU+J)-YA)**2.0D0)**0.5D0
                END IF
                IF(ELIGIBILITY_V(I,J))THEN
                    DISTANCE_V(I,J)=((XPV(IAV+I)-XA)**2.0D0&
                        +(Y(JAV+J)-YA)**2.0D0)**0.5D0
                END IF
            END DO
        END DO

        !------确定搜索半径最小值------!
        SEARCH_RADIUS_U=0.0D0
        SEARCH_RADIUS_V=0.0D0
        DO J=0,1,1
            DO I=0,1,1
                IF(ELIGIBILITY_U(I,J))THEN
                    SEARCH_RADIUS_U=MAX(SEARCH_RADIUS_U,DISTANCE_U(I,J))
                END IF
                IF(ELIGIBILITY_V(I,J))THEN
                    SEARCH_RADIUS_V=MAX(SEARCH_RADIUS_V,DISTANCE_V(I,J))
                END IF
            END DO
        END DO

        !------确定合理搜索半径及搜索半径之内的点------!
        !还可以保证坐标偏移量
        !U
        DO WHILE (.TRUE.)

            NUM_ELGBL_U=0
            MAX_X=XA
            MAX_Y=YA
            MIN_X=XA
            MIN_Y=YA
            SEARCH_RADIUS_U=SEARCH_RADIUS_U+0.1D0*DX3
            !确定点数
            DO J=-2,+2,1
                DO I=-2,+2,1

                    IF(DISTANCE_U(I,J)<=SEARCH_RADIUS_U)THEN
                        NUM_ELGBL_U=NUM_ELGBL_U+1
                        ELIGIBILITY_U(I,J)=.TRUE.
                        MAX_X=MAX(MAX_X,X  (I+IAU))
                        MAX_Y=MAX(MAX_Y,YPU(J+JAU))
                        MIN_X=MIN(MIN_X,X  (I+IAU))
                        MIN_Y=MIN(MIN_Y,YPU(J+JAU))
                    ELSE
                        ELIGIBILITY_U(I,J)=.FALSE.
                    END IF

                END DO
            END DO
            IF ( NUM_ELGBL_U>=2 ) THEN
                IF( MAX_X-MIN_X>=1000.0D0*CRITERIA .AND.&
                    MAX_Y-MIN_Y>=1000.0D0*CRITERIA ) EXIT
            END IF

            IF ( SEARCH_RADIUS_U>2.1D0*DX3 ) THEN
                WRITE(*,*)N,'u error: not enough eligible points'
                STOP
            END IF

        END DO
        !V
        DO WHILE (.TRUE.)

            NUM_ELGBL_V=0
            MAX_X=XA
            MAX_Y=YA
            MIN_X=XA
            MIN_Y=YA
            SEARCH_RADIUS_V=SEARCH_RADIUS_V+0.1D0*DX3
            !确定点数
            DO J=-2,+2,1
                DO I=-2,+2,1

                    IF(DISTANCE_V(I,J)<=SEARCH_RADIUS_V)THEN
                        NUM_ELGBL_V=NUM_ELGBL_V+1
                        ELIGIBILITY_V(I,J)=.TRUE.
                        MAX_X=MAX(MAX_X,XPV(I+IAV))
                        MAX_Y=MAX(MAX_Y,Y  (J+JAV))
                        MIN_X=MIN(MIN_X,XPV(I+IAV))
                        MIN_Y=MIN(MIN_Y,Y  (J+JAV))
                    ELSE
                        ELIGIBILITY_V(I,J)=.FALSE.
                    END IF

                END DO
            END DO
            IF ( NUM_ELGBL_V>=2 ) THEN
                IF( MAX_X-MIN_X>=1000.0D0*CRITERIA .AND.&
                    MAX_Y-MIN_Y>=1000.0D0*CRITERIA ) EXIT
            END IF

            IF ( SEARCH_RADIUS_V>2.1D0*DX3 ) THEN
                WRITE(*,*)N,'v error: not enough eligible points'
                STOP
            END IF

        END DO
        !分配矩阵
        IF(ALLOCATED(INDEX_U)) DEALLOCATE(INDEX_U,INDEX_V)
        ALLOCATE( INDEX_U(NUM_ELGBL_U,2),INDEX_V(NUM_ELGBL_V,2) )
        !赋值矩阵
        COUNT_ELGBL=0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    COUNT_ELGBL=COUNT_ELGBL+1
                    INDEX_U(COUNT_ELGBL,1)=I
                    INDEX_U(COUNT_ELGBL,2)=J
                END IF
            END DO
        END DO
        INDEX_U(:,1)=INDEX_U(:,1)+IAU
        INDEX_U(:,2)=INDEX_U(:,2)+JAU

        COUNT_ELGBL=0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_V(I,J))THEN
                    COUNT_ELGBL=COUNT_ELGBL+1
                    INDEX_V(COUNT_ELGBL,1)=I
                    INDEX_V(COUNT_ELGBL,2)=J
                END IF
            END DO
        END DO
        INDEX_V(:,1)=INDEX_V(:,1)+IAV
        INDEX_V(:,2)=INDEX_V(:,2)+JAV

        !------构造系数矩阵并求解------!
        !U
        IF(ALLOCATED(TEMP_1)) DEALLOCATE(TEMP_1,TEMP_2)
        ALLOCATE( TEMP_1(NUM_ELGBL_U+1,3),TEMP_2(NUM_ELGBL_U+1,1) )

        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(1,3)=YA
        DO COUNT_ELGBL=1,NUM_ELGBL_U,1
            TEMP_1(1+COUNT_ELGBL,2)=X(INDEX_U(COUNT_ELGBL,1))
            TEMP_1(1+COUNT_ELGBL,3)=YPU(INDEX_U(COUNT_ELGBL,2))
        END DO
        TEMP_2(1,1)=U_SRFC(N)
        DO COUNT_ELGBL=1,NUM_ELGBL_U,1
            TEMP_2(1+COUNT_ELGBL,1)=&
                U(INDEX_U(COUNT_ELGBL,1),INDEX_U(COUNT_ELGBL,2))
        END DO
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dy'
        END IF
        !求解
        TEMP_SOLUTION=0.0D0
        CALL SOLVE_LINEAR_SYSTEM(TEMP_1,TEMP_SOLUTION,TEMP_2,&
            NUM_ELGBL_U+1,3,N)
        AXX(N)=TEMP_SOLUTION(2,1)
        AXY(N)=TEMP_SOLUTION(3,1)

        !V
        IF(ALLOCATED(TEMP_1)) DEALLOCATE(TEMP_1,TEMP_2)
        ALLOCATE( TEMP_1(NUM_ELGBL_V+1,3),TEMP_2(NUM_ELGBL_V+1,1) )

        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(1,3)=YA
        DO COUNT_ELGBL=1,NUM_ELGBL_V,1
            TEMP_1(1+COUNT_ELGBL,2)=XPV(INDEX_V(COUNT_ELGBL,1))
            TEMP_1(1+COUNT_ELGBL,3)=Y(INDEX_V(COUNT_ELGBL,2))
        END DO
        TEMP_2(1,1)=V_SRFC(N)
        DO COUNT_ELGBL=1,NUM_ELGBL_V,1
            TEMP_2(1+COUNT_ELGBL,1)=&
                V(INDEX_V(COUNT_ELGBL,1),INDEX_V(COUNT_ELGBL,2))
        END DO
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dy'
        END IF
        !求解
        TEMP_SOLUTION=0.0D0
        CALL SOLVE_LINEAR_SYSTEM(TEMP_1,TEMP_SOLUTION,TEMP_2,&
            NUM_ELGBL_V+1,3,N)
        AYX(N)=TEMP_SOLUTION(2,1)
        AYY(N)=TEMP_SOLUTION(3,1)

    END DO

    !------求解气动力------!
    DO N=1,RSM,1

        CXP=-P_SRFC(N)*N1X(N)
        CYP=-P_SRFC(N)*N1Y(N)
        CZP=0.0D0!-P_SRFC*N1Z(N)
        CXV=(2.0D0*AXX(N)*N1X(N)+(AXY(N)+AYX(N))*N1Y(N))/Re
        CYV=((AXY(N)+AYX(N))*N1X(N)+2.0D0*AYY(N)*N1Y(N))/Re
        CZV=0.0D0

        CXC_TOTAL=CXC_TOTAL+2.0D0*DS(N)*(CXP+CXV)
        CYC_TOTAL=CYC_TOTAL+2.0D0*DS(N)*(CYP+CYV)

    END DO

    !------物体所受气动力为反作用力------!
    CXC_TOTAL=-CXC_TOTAL
    CYC_TOTAL=-CYC_TOTAL

    !------输出系数------!
    WRITE(360+BOUNDARY_ID*10,"( I6,2(1X,F9.5))")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(361+BOUNDARY_ID*10,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    !------输出分布------!
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS
    WRITE(CHAR_ID,'(I2.2)') BOUNDARY_ID
    IF( MOD(NSTEP,2000)==0 )THEN
        OPEN(UNIT=10,FILE='SRFC_SOME_OBJECT_'//TRIM(CHAR_ID)&
            //'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
        WRITE(10,*) 'TITLE="NONAME"'
        !WRITE(10,*) 'VARIABLES="X","Y","P","U","V","AXX","AXY","AXZ",'
        !WRITE(10,*) '"AYX","AYY","AYZ","AZX","AZY","AZZ"'
        WRITE(10,*) 'VARIABLES="X","Y","P","U","V","NX","NY","AXX","AXY",'
        WRITE(10,*) '"AYX","AYY"'
        WRITE(10,*) 'ZONE T="NONAME", I=',RSM,', J=',1,', F=POINT'
        DO N=1,RSM,1
            !WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
            !    P_SRFC(N),U_SRFC(N),V_SRFC(N),&
            !    AXX(N),AXY(N),AXZ(N),&
            !    AYX(N),AYY(N),AYZ(N),&
            !    AZX(N),AZY(N),AZZ(N)
            WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
                P_SRFC(N),U_SRFC(N),V_SRFC(N),&
                N1X(N),N1Y(N),&
                AXX(N),AXY(N),&
                AYX(N),AYY(N)
        END DO
        CLOSE(10)
    END IF

    RETURN
    END SUBROUTINE

    !*********************求解升力推力（表面力积分法3）******************!
    SUBROUTINE CAL_CLCT_SRFC_INTGRTN_ALL(BOUNDARY_GEOMETRICAL,&
        BOUNDARY_KINETIC,BOUNDARY_ID)
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::BOUNDARY_GEOMETRICAL(36),BOUNDARY_KINETIC(9)
    INTEGER::BOUNDARY_ID
    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    CHARACTER(LEN=2)CHAR_ID
    INTEGER::REYNOLDS

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM,COZ2,COXZ,COYZ,COZ
    !网格密度
    REAL(KIND=8),ALLOCATABLE::DS(:)
    !网格数
    INTEGER::RSM
    !插值涉及的绝对网格角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !坐标值
    REAL(KIND=8),ALLOCATABLE::X_SRFC(:),Y_SRFC(:),Z_SRFC(:)
    REAL(KIND=8)::XA,YA
    !法向量
    REAL(KIND=8),ALLOCATABLE::N1X(:),N1Y(:),N1Z(:)
    !偏导系数
    REAL(KIND=8),ALLOCATABLE::AXX(:),AXY(:),AXZ(:)
    REAL(KIND=8),ALLOCATABLE::AYX(:),AYY(:),AYZ(:)
    REAL(KIND=8),ALLOCATABLE::AZX(:),AZY(:),AZZ(:)
    !固壁处物理量
    REAL(KIND=8),ALLOCATABLE::U_SRFC(:),V_SRFC(:),W_SRFC(:),P_SRFC(:)
    !包括半径
    REAL(KIND=8)::SEARCH_RADIUS
    !包括点数
    INTEGER::NUM_ELGBL_U,NUM_ELGBL_V,COUNT_ELGBL

    !一些中间变量
    REAL(KIND=8)::ELLIPTIC_H,ELLIPTIC_PERIMETER!椭圆周长
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_DA(:)!椭圆各段弧对应圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AS(:)!椭圆各段弧起始处的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AM(:)!椭圆各段弧中间的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_RADIUS(:)!椭圆弧长对应半径
    REAL(KIND=8)::AVERAGE_RADIUS!椭圆平均半径
    LOGICAL::ELIGIBILITY_U(-2:2,-2:2),ELIGIBILITY_V(-2:2,-2:2)!可视流体点
    REAL(KIND=8)::DISTANCE_U(-2:2,-2:2),DISTANCE_V(-2:2,-2:2)!各点距离
    INTEGER,ALLOCATABLE::INDEX_U(:,:),INDEX_V(:,:)!求解偏导的几个点的脚标
    REAL(KIND=8),ALLOCATABLE::TEMP_1(:,:),TEMP_2(:,:)!偏导系数求解临时矩阵
    REAL(KIND=8)::TEMP_SOLUTION(3,1)!求解结果临时矩阵
    LOGICAL :: OK_FLAG
    REAL(KIND=8)::DEVIATION_X,DEVIATION_Y!坐标偏移
    REAL(KIND=8)::XATEMP,YATEMP

    !对象特殊性
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::CEN_ELP(2)
    REAL(KIND=8)::CEN_DEVIATION(2)
    REAL(KIND=8)::LAXIS,SAXIS
    !坐标转换矩阵
    REAL(KIND=8)::MAT_ABS2REL(2,2),MAT_REL2ABS(2,2)

    !力系数
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL,CZC_TOTAL!计算域绝对坐标系下的力
    REAL(KIND=8)::CXP,CYP,CZP!压力
    REAL(KIND=8)::CXV,CYV,CZV!粘性

    !------初始化------!
    COX2=BOUNDARY_GEOMETRICAL(7)
    COY2=BOUNDARY_GEOMETRICAL(8)
    COXY=BOUNDARY_GEOMETRICAL(9)
    COX =BOUNDARY_GEOMETRICAL(10)
    COY =BOUNDARY_GEOMETRICAL(11)
    COM =BOUNDARY_GEOMETRICAL(13)
    COZ2=0.0D0
    COXZ=0.0D0
    COYZ=0.0D0
    COZ =0.0D0

    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0

    SEARCH_RADIUS=(2.0D0)*DX3!(2.0D0)**0.5D0*DX3

    !------离散界面------!
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!
    IF(IB_SHAPE==1)THEN!1圆
        RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)&
            **0.5D0
        RSM=IDNINT( 2.0D0*RADIUS*PI/DX3 )

        ALLOCATE( DS(RSM) )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        DS=2.0D0*RADIUS*PI/DBLE(RSM)

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        CEN_CRC(1)=BOUNDARY_KINETIC(4)
        CEN_CRC(2)=BOUNDARY_KINETIC(5)

        DO I=1,RSM,1
            X_SRFC(I)=CEN_CRC(1)+RADIUS*DCOS(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
            Y_SRFC(I)=CEN_CRC(2)+RADIUS*DSIN(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
        END DO

    ELSE IF(IB_SHAPE==2)THEN!2椭圆
        LAXIS=BOUNDARY_GEOMETRICAL(33)
        SAXIS=BOUNDARY_GEOMETRICAL(34)

        CEN_DEVIATION(1)=BOUNDARY_GEOMETRICAL(35)
        CEN_DEVIATION(2)=BOUNDARY_GEOMETRICAL(36)

        CEN_ELP(1)=BOUNDARY_KINETIC(4)
        CEN_ELP(2)=BOUNDARY_KINETIC(5)

        MAT_ABS2REL(1,1)=BOUNDARY_GEOMETRICAL(1)
        MAT_ABS2REL(1,2)=BOUNDARY_GEOMETRICAL(2)
        MAT_ABS2REL(2,1)=BOUNDARY_GEOMETRICAL(4)
        MAT_ABS2REL(2,2)=BOUNDARY_GEOMETRICAL(5)
        MAT_REL2ABS=TRANSPOSE(MAT_ABS2REL)

        !估求椭圆周长
        ELLIPTIC_H=((LAXIS-SAXIS)/(LAXIS+SAXIS))**2.0D0
        ELLIPTIC_H=3.0D0*ELLIPTIC_H/(10.0D0+(4.0D0-3.0D0*ELLIPTIC_H)**0.5D0)
        ELLIPTIC_PERIMETER=PI*(LAXIS+SAXIS)*(1.0D0+ELLIPTIC_H)
        !确定网格数
        RSM=IDNINT( ELLIPTIC_PERIMETER/DX3 )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        ALLOCATE( DS(RSM) )

        ALLOCATE( ELLIPTIC_DA(RSM) )
        ALLOCATE( ELLIPTIC_AS(RSM),ELLIPTIC_AM(RSM),ELLIPTIC_RADIUS(RSM) )

        DS=0.0D0
        ELLIPTIC_DA=0.0D0
        ELLIPTIC_AM=0.0D0
        ELLIPTIC_AS=0.0D0

        !确定平均半径周长/2PI
        AVERAGE_RADIUS=ELLIPTIC_PERIMETER/PI/2.0D0

        !初始化
        !ELLIPTIC_A自x轴正方向为零，顺时针方向为正方向
        DO I=1,RSM,1
            ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)
            ELLIPTIC_AS(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.5D0)
            ELLIPTIC_AM(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.0D0)
            ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
        END DO

        !迭代确定各段弧
        DO N=1,10,1

            !修正圆心角
            DO I=1,RSM,1
                ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)*AVERAGE_RADIUS/ELLIPTIC_RADIUS(I)
            END DO
            ELLIPTIC_DA=ELLIPTIC_DA*2.0D0*PI/SUM(ELLIPTIC_DA)
            !重新获得ELLIPTIC_AS和ELLIPTIC_AM
            ELLIPTIC_AS(1)=-ELLIPTIC_DA(1)/2.0D0
            ELLIPTIC_AM(1)=0.0D0
            DO I=2,RSM,1
                ELLIPTIC_AS(I)=ELLIPTIC_AS(I-1)+ELLIPTIC_DA(I-1)
                ELLIPTIC_AM(I)=ELLIPTIC_AS(I)+ELLIPTIC_DA(I)/2.0D0
            END DO
            !先求起始点坐标
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AS(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AS(I)))**2.0D0 )**0.5D0
            END DO
            DO I=1,RSM,1
                X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS( ELLIPTIC_AS(I) )
                Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN( ELLIPTIC_AS(I) )
            END DO
            !求解弧长
            DO I=1,RSM-1,1
                DS(I)=DSQRT(&
                    (X_SRFC(I+1)-X_SRFC(I))**2.0D0+&
                    (Y_SRFC(I+1)-Y_SRFC(I))**2.0D0)
            END DO
            DS(RSM)=DSQRT(&
                (X_SRFC(1)-X_SRFC(RSM))**2.0D0+&
                (Y_SRFC(1)-Y_SRFC(RSM))**2.0D0)
            !检查弧长差异,可以的话就跳出
            IF(MAXVAL(DS)/MINVAL(DS)<1.1D0)EXIT

            !不行的话再重新修正
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
            END DO

        END DO
        !检查弧长差异,太过了就报错
        IF(MAXVAL(DS)>=2.0D0*DX3)THEN
            WRITE(*,*)'too large max(ds)/dx=',MAXVAL(DS)/DX3
            STOP
        END IF
        !求各弧中心点坐标
        DO I=1,RSM,1
            X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(1)
            Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(2)
            CALL CRDNT_TRANSFORM(MAT_REL2ABS,X_SRFC(I),Y_SRFC(I),&
                XATEMP,YATEMP)
            X_SRFC(I)=XATEMP+CEN_ELP(1)
            Y_SRFC(I)=YATEMP+CEN_ELP(2)
        END DO

    END IF
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!

    !------求解内法向量------!
    DO I=1,RSM,1
        CALL NORMALVECTOR_2D(BOUNDARY_GEOMETRICAL,X_SRFC(I),Y_SRFC(I),&
            N1X(I),N1Y(I))
        N1X(I)=-N1X(I)
        N1Y(I)=-N1Y(I)
        N1Z(I)=0.0D0
    END DO
    !------求解压力和偏导系数------!
    DO N=1,RSM,1

        !------求解点所处位置及压力------!
        XA=X_SRFC(N)
        YA=Y_SRFC(N)

        IF( XA>LEIN-CRITERIA .AND. XA<RIIN+CRITERIA )THEN
            IAU=IL+FLOOR( (XA-LEIN) / DX3  )
            IAV=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
            IAP=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"XA OUT OF INNER REGION"
            STOP
        END IF
        IF( YA>BOIN-CRITERIA .AND. YA<TOIN+CRITERIA )THEN
            JAU=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
            JAV=JB+FLOOR( (YA-BOIN) / DX3  )
            JAP=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"YA OUT OF INNER REGION"
            STOP
        END IF

        CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA,YA,&
            U_SRFC(N),V_SRFC(N),P_SRFC(N))

        !------求解点刚体速度------!
        CALL VELOCITY_LB(BOUNDARY_KINETIC,XA,YA,U_SRFC(N),V_SRFC(N))

        !------获取附近符合要求的点------!
        CALL DETERMINE_IF_ELIGIBLE_U(IAU,JAU,ELIGIBILITY_U)
        CALL DETERMINE_IF_ELIGIBLE_V(IAV,JAV,ELIGIBILITY_V)

        !------求解距离------!
        DISTANCE_U=500.0D0
        DISTANCE_V=500.0D0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    DISTANCE_U(I,J)=(X(IAU+I)-XA)**2.0D0&
                        +(YPU(JAU+J)-YA)**2.0D0
                END IF
                IF(ELIGIBILITY_V(I,J))THEN
                    DISTANCE_V(I,J)=(XPV(IAV+I)-XA)**2.0D0&
                        +(Y(JAV+J)-YA)**2.0D0
                END IF
            END DO
        END DO

        !------确定2*DX3之内的点------!
        NUM_ELGBL_U=0
        NUM_ELGBL_V=0
        !确定点数
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    IF(DISTANCE_U(I,J)**0.5D0<=SEARCH_RADIUS)THEN
                        NUM_ELGBL_U=NUM_ELGBL_U+1
                    ELSE
                        ELIGIBILITY_U(I,J)=.FALSE.
                    END IF
                END IF
                IF(ELIGIBILITY_V(I,J))THEN
                    IF(DISTANCE_V(I,J)**0.5D0<=SEARCH_RADIUS)THEN
                        NUM_ELGBL_V=NUM_ELGBL_V+1
                    ELSE
                        ELIGIBILITY_V(I,J)=.FALSE.
                    END IF
                END IF
            END DO
        END DO
        IF ( NUM_ELGBL_U<2 ) THEN
            WRITE(*,*)N,'u error: not enough eligible points'
            STOP
        END IF
        IF ( NUM_ELGBL_V<2 ) THEN
            WRITE(*,*)N,'v error: not enough eligible points'
            STOP
        END IF
        !分配矩阵
        IF(ALLOCATED(INDEX_U)) DEALLOCATE(INDEX_U,INDEX_V)
        ALLOCATE( INDEX_U(NUM_ELGBL_U,2),INDEX_V(NUM_ELGBL_V,2) )
        !赋值矩阵
        COUNT_ELGBL=0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    COUNT_ELGBL=COUNT_ELGBL+1
                    INDEX_U(COUNT_ELGBL,1)=I
                    INDEX_U(COUNT_ELGBL,2)=J
                END IF
            END DO
        END DO
        INDEX_U(:,1)=INDEX_U(:,1)+IAU
        INDEX_U(:,2)=INDEX_U(:,2)+JAU

        COUNT_ELGBL=0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_V(I,J))THEN
                    COUNT_ELGBL=COUNT_ELGBL+1
                    INDEX_V(COUNT_ELGBL,1)=I
                    INDEX_V(COUNT_ELGBL,2)=J
                END IF
            END DO
        END DO
        INDEX_V(:,1)=INDEX_V(:,1)+IAV
        INDEX_V(:,2)=INDEX_V(:,2)+JAV

        !------构造系数矩阵并求解------!
        !U
        IF(ALLOCATED(TEMP_1)) DEALLOCATE(TEMP_1,TEMP_2)
        ALLOCATE( TEMP_1(NUM_ELGBL_U+1,3),TEMP_2(NUM_ELGBL_U+1,1) )

        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(1,3)=YA
        DO COUNT_ELGBL=1,NUM_ELGBL_U,1
            TEMP_1(1+COUNT_ELGBL,2)=X(INDEX_U(COUNT_ELGBL,1))
            TEMP_1(1+COUNT_ELGBL,3)=YPU(INDEX_U(COUNT_ELGBL,2))
        END DO
        TEMP_2(1,1)=U_SRFC(N)
        DO COUNT_ELGBL=1,NUM_ELGBL_U,1
            TEMP_2(1+COUNT_ELGBL,1)=&
                U(INDEX_U(COUNT_ELGBL,1),INDEX_U(COUNT_ELGBL,2))
        END DO
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dy'
        END IF
        !求解
        TEMP_SOLUTION=0.0D0
        CALL SOLVE_LINEAR_SYSTEM(TEMP_1,TEMP_SOLUTION,TEMP_2,&
            NUM_ELGBL_U+1,3,N)
        AXX(N)=TEMP_SOLUTION(2,1)
        AXY(N)=TEMP_SOLUTION(3,1)

        !V
        IF(ALLOCATED(TEMP_1)) DEALLOCATE(TEMP_1,TEMP_2)
        ALLOCATE( TEMP_1(NUM_ELGBL_V+1,3),TEMP_2(NUM_ELGBL_V+1,1) )

        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(1,3)=YA
        DO COUNT_ELGBL=1,NUM_ELGBL_V,1
            TEMP_1(1+COUNT_ELGBL,2)=XPV(INDEX_V(COUNT_ELGBL,1))
            TEMP_1(1+COUNT_ELGBL,3)=Y(INDEX_V(COUNT_ELGBL,2))
        END DO
        TEMP_2(1,1)=V_SRFC(N)
        DO COUNT_ELGBL=1,NUM_ELGBL_V,1
            TEMP_2(1+COUNT_ELGBL,1)=&
                V(INDEX_V(COUNT_ELGBL,1),INDEX_V(COUNT_ELGBL,2))
        END DO
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dy'
        END IF
        !求解
        TEMP_SOLUTION=0.0D0
        CALL SOLVE_LINEAR_SYSTEM(TEMP_1,TEMP_SOLUTION,TEMP_2,&
            NUM_ELGBL_V+1,3,N)
        AYX(N)=TEMP_SOLUTION(2,1)
        AYY(N)=TEMP_SOLUTION(3,1)

    END DO

    !------求解气动力------!
    DO N=1,RSM,1

        CXP=-P_SRFC(N)*N1X(N)
        CYP=-P_SRFC(N)*N1Y(N)
        CZP=0.0D0!-P_SRFC*N1Z(N)
        CXV=(2.0D0*AXX(N)*N1X(N)+(AXY(N)+AYX(N))*N1Y(N))/Re
        CYV=((AXY(N)+AYX(N))*N1X(N)+2.0D0*AYY(N)*N1Y(N))/Re
        CZV=0.0D0

        CXC_TOTAL=CXC_TOTAL+2.0D0*DS(N)*(CXP+CXV)
        CYC_TOTAL=CYC_TOTAL+2.0D0*DS(N)*(CYP+CYV)

    END DO

    !------物体所受气动力为反作用力------!
    CXC_TOTAL=-CXC_TOTAL
    CYC_TOTAL=-CYC_TOTAL

    !------输出系数------!
    WRITE(340+BOUNDARY_ID*10,"( I6,2(1X,F9.5))")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(341+BOUNDARY_ID*10,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    !------输出分布------!
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS
    WRITE(CHAR_ID,'(I2.2)') BOUNDARY_ID
    IF( MOD(NSTEP,2000)==0 )THEN
        OPEN(UNIT=10,FILE='SRFC_ALL_OBJECT_'//TRIM(CHAR_ID)&
            //'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
        WRITE(10,*) 'TITLE="NONAME"'
        !WRITE(10,*) 'VARIABLES="X","Y","P","U","V","AXX","AXY","AXZ",'
        !WRITE(10,*) '"AYX","AYY","AYZ","AZX","AZY","AZZ"'
        WRITE(10,*) 'VARIABLES="X","Y","P","U","V","NX","NY","AXX","AXY",'
        WRITE(10,*) '"AYX","AYY"'
        WRITE(10,*) 'ZONE T="NONAME", I=',RSM,', J=',1,', F=POINT'
        DO N=1,RSM,1
            !WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
            !    P_SRFC(N),U_SRFC(N),V_SRFC(N),&
            !    AXX(N),AXY(N),AXZ(N),&
            !    AYX(N),AYY(N),AYZ(N),&
            !    AZX(N),AZY(N),AZZ(N)
            WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
                P_SRFC(N),U_SRFC(N),V_SRFC(N),&
                N1X(N),N1Y(N),&
                AXX(N),AXY(N),&
                AYX(N),AYY(N)
        END DO
        CLOSE(10)
    END IF

    RETURN
    END SUBROUTINE

    !*********************求解升力推力（表面力积分法2）******************!
    SUBROUTINE CAL_CLCT_SRFC_INTGRTN_EXACT(BOUNDARY_GEOMETRICAL,&
        BOUNDARY_KINETIC,BOUNDARY_ID)
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::BOUNDARY_GEOMETRICAL(36),BOUNDARY_KINETIC(9)
    INTEGER::BOUNDARY_ID
    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    CHARACTER(LEN=2)CHAR_ID
    INTEGER::REYNOLDS

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM,COZ2,COXZ,COYZ,COZ
    !网格密度
    REAL(KIND=8),ALLOCATABLE::DS(:)
    !网格数
    INTEGER::RSM
    !插值涉及的绝对网格角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !坐标值
    REAL(KIND=8),ALLOCATABLE::X_SRFC(:),Y_SRFC(:),Z_SRFC(:)
    REAL(KIND=8)::XA,YA
    !法向量
    REAL(KIND=8),ALLOCATABLE::N1X(:),N1Y(:),N1Z(:)
    !偏导系数
    REAL(KIND=8),ALLOCATABLE::AXX(:),AXY(:),AXZ(:)
    REAL(KIND=8),ALLOCATABLE::AYX(:),AYY(:),AYZ(:)
    REAL(KIND=8),ALLOCATABLE::AZX(:),AZY(:),AZZ(:)
    !固壁处物理量
    REAL(KIND=8),ALLOCATABLE::U_SRFC(:),V_SRFC(:),W_SRFC(:),P_SRFC(:)

    !一些中间变量
    REAL(KIND=8)::ELLIPTIC_H,ELLIPTIC_PERIMETER!椭圆周长
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_DA(:)!椭圆各段弧对应圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AS(:)!椭圆各段弧起始处的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_AM(:)!椭圆各段弧中间的圆心角
    REAL(KIND=8),ALLOCATABLE::ELLIPTIC_RADIUS(:)!椭圆弧长对应半径
    REAL(KIND=8)::AVERAGE_RADIUS!椭圆平均半径
    LOGICAL::ELIGIBILITY_U(-2:2,-2:2),ELIGIBILITY_V(-2:2,-2:2)!可视流体点
    REAL(KIND=8)::DISTANCE_U(-2:2,-2:2),DISTANCE_V(-2:2,-2:2)!各点距离
    INTEGER::INDEX_U(2,2),INDEX_V(2,2)!最终求解偏导的两个点的脚标
    REAL(KIND=8)::TEMP_1(3,3),TEMP_2(3,3),TEMP_3(3,1),TEMP_4(3,1)!偏导系数
    LOGICAL :: OK_FLAG
    REAL(KIND=8)::DEVIATION_X,DEVIATION_Y!坐标偏移
    REAL(KIND=8)::XATEMP,YATEMP

    !对象特殊性
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::CEN_ELP(2)
    REAL(KIND=8)::CEN_DEVIATION(2)
    REAL(KIND=8)::LAXIS,SAXIS
    !坐标转换矩阵
    REAL(KIND=8)::MAT_ABS2REL(2,2),MAT_REL2ABS(2,2)

    !力系数
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL,CZC_TOTAL!计算域绝对坐标系下的力
    REAL(KIND=8)::CXP,CYP,CZP!压力
    REAL(KIND=8)::CXV,CYV,CZV!粘性

    !------初始化------!
    COX2=BOUNDARY_GEOMETRICAL(7)
    COY2=BOUNDARY_GEOMETRICAL(8)
    COXY=BOUNDARY_GEOMETRICAL(9)
    COX =BOUNDARY_GEOMETRICAL(10)
    COY =BOUNDARY_GEOMETRICAL(11)
    COM =BOUNDARY_GEOMETRICAL(13)
    COZ2=0.0D0
    COXZ=0.0D0
    COYZ=0.0D0
    COZ =0.0D0

    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0

    !------离散界面------!
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!
    IF(IB_SHAPE==1)THEN!1圆
        RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)&
            **0.5D0
        RSM=IDNINT( 2.0D0*RADIUS*PI/DX3 )

        ALLOCATE( DS(RSM) )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        DS=2.0D0*RADIUS*PI/DBLE(RSM)

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        CEN_CRC(1)=BOUNDARY_KINETIC(4)
        CEN_CRC(2)=BOUNDARY_KINETIC(5)

        DO I=1,RSM,1
            X_SRFC(I)=CEN_CRC(1)+RADIUS*DCOS(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
            Y_SRFC(I)=CEN_CRC(2)+RADIUS*DSIN(-2.0D0*PI*DBLE(I-1)/DBLE(RSM))
        END DO

    ELSE IF(IB_SHAPE==2)THEN!2椭圆
        LAXIS=BOUNDARY_GEOMETRICAL(33)
        SAXIS=BOUNDARY_GEOMETRICAL(34)

        CEN_DEVIATION(1)=BOUNDARY_GEOMETRICAL(35)
        CEN_DEVIATION(2)=BOUNDARY_GEOMETRICAL(36)

        CEN_ELP(1)=BOUNDARY_KINETIC(4)
        CEN_ELP(2)=BOUNDARY_KINETIC(5)

        MAT_ABS2REL(1,1)=BOUNDARY_GEOMETRICAL(1)
        MAT_ABS2REL(1,2)=BOUNDARY_GEOMETRICAL(2)
        MAT_ABS2REL(2,1)=BOUNDARY_GEOMETRICAL(4)
        MAT_ABS2REL(2,2)=BOUNDARY_GEOMETRICAL(5)
        MAT_REL2ABS=TRANSPOSE(MAT_ABS2REL)

        !估求椭圆周长
        ELLIPTIC_H=((LAXIS-SAXIS)/(LAXIS+SAXIS))**2.0D0
        ELLIPTIC_H=3.0D0*ELLIPTIC_H/(10.0D0+(4.0D0-3.0D0*ELLIPTIC_H)**0.5D0)
        ELLIPTIC_PERIMETER=PI*(LAXIS+SAXIS)*(1.0D0+ELLIPTIC_H)
        !确定网格数
        RSM=IDNINT( ELLIPTIC_PERIMETER/DX3 )

        ALLOCATE( X_SRFC(RSM),Y_SRFC(RSM),Z_SRFC(RSM) )
        ALLOCATE( N1X(RSM),N1Y(RSM),N1Z(RSM) )

        ALLOCATE( AXX(RSM),AXY(RSM),AXZ(RSM) )
        ALLOCATE( AYX(RSM),AYY(RSM),AYZ(RSM) )
        ALLOCATE( AZX(RSM),AZY(RSM),AZZ(RSM) )
        ALLOCATE( U_SRFC(RSM),V_SRFC(RSM),W_SRFC(RSM),P_SRFC(RSM) )

        X_SRFC=0.0D0
        Y_SRFC=0.0D0
        Z_SRFC=0.0D0
        N1X=0.0D0
        N1Y=0.0D0
        N1Z=0.0D0

        AXX=0.0D0
        AXY=0.0D0
        AXZ=0.0D0
        AYX=0.0D0
        AYY=0.0D0
        AYZ=0.0D0
        AZX=0.0D0
        AZY=0.0D0
        AZZ=0.0D0
        U_SRFC=0.0D0
        V_SRFC=0.0D0
        W_SRFC=0.0D0
        P_SRFC=0.0D0

        ALLOCATE( DS(RSM) )

        ALLOCATE( ELLIPTIC_DA(RSM) )
        ALLOCATE( ELLIPTIC_AS(RSM),ELLIPTIC_AM(RSM),ELLIPTIC_RADIUS(RSM) )

        DS=0.0D0
        ELLIPTIC_DA=0.0D0
        ELLIPTIC_AM=0.0D0
        ELLIPTIC_AS=0.0D0

        !确定平均半径周长/2PI
        AVERAGE_RADIUS=ELLIPTIC_PERIMETER/PI/2.0D0

        !初始化
        !ELLIPTIC_A自x轴正方向为零，顺时针方向为正方向
        DO I=1,RSM,1
            ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)
            ELLIPTIC_AS(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.5D0)
            ELLIPTIC_AM(I)=2.0D0*PI/DBLE(RSM)*(DBLE(I)-1.0D0)
            ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
        END DO

        !迭代确定各段弧
        DO N=1,10,1

            !修正圆心角
            DO I=1,RSM,1
                ELLIPTIC_DA(I)=2.0D0*PI/DBLE(RSM)*AVERAGE_RADIUS/ELLIPTIC_RADIUS(I)
            END DO
            ELLIPTIC_DA=ELLIPTIC_DA*2.0D0*PI/SUM(ELLIPTIC_DA)
            !重新获得ELLIPTIC_AS和ELLIPTIC_AM
            ELLIPTIC_AS(1)=-ELLIPTIC_DA(1)/2.0D0
            ELLIPTIC_AM(1)=0.0D0
            DO I=2,RSM,1
                ELLIPTIC_AS(I)=ELLIPTIC_AS(I-1)+ELLIPTIC_DA(I-1)
                ELLIPTIC_AM(I)=ELLIPTIC_AS(I)+ELLIPTIC_DA(I)/2.0D0
            END DO
            !先求起始点坐标
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AS(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AS(I)))**2.0D0 )**0.5D0
            END DO
            DO I=1,RSM,1
                X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS( ELLIPTIC_AS(I) )
                Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN( ELLIPTIC_AS(I) )
            END DO
            !求解弧长
            DO I=1,RSM-1,1
                DS(I)=DSQRT(&
                    (X_SRFC(I+1)-X_SRFC(I))**2.0D0+&
                    (Y_SRFC(I+1)-Y_SRFC(I))**2.0D0)
            END DO
            DS(RSM)=DSQRT(&
                (X_SRFC(1)-X_SRFC(RSM))**2.0D0+&
                (Y_SRFC(1)-Y_SRFC(RSM))**2.0D0)
            !检查弧长差异,可以的话就跳出
            IF(MAXVAL(DS)/MINVAL(DS)<1.1D0)EXIT

            !不行的话再重新修正
            DO I=1,RSM,1
                ELLIPTIC_RADIUS(I)=LAXIS*SAXIS/&
                    ( (SAXIS*DCOS(ELLIPTIC_AM(I)))**2.0D0  &
                    + (LAXIS*DSIN(ELLIPTIC_AM(I)))**2.0D0 )**0.5D0
            END DO

        END DO
        !检查弧长差异,太过了就报错
        IF(MAXVAL(DS)>=2.0D0*DX3)THEN
            WRITE(*,*)'too large max(ds)/dx=',MAXVAL(DS)/DX3
            STOP
        END IF
        !求各弧中心点坐标
        DO I=1,RSM,1
            X_SRFC(I)=ELLIPTIC_RADIUS(I)*DCOS(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(1)
            Y_SRFC(I)=ELLIPTIC_RADIUS(I)*DSIN(ELLIPTIC_AM(I))&
                +CEN_DEVIATION(2)
            CALL CRDNT_TRANSFORM(MAT_REL2ABS,X_SRFC(I),Y_SRFC(I),&
                XATEMP,YATEMP)
            X_SRFC(I)=XATEMP+CEN_ELP(1)
            Y_SRFC(I)=YATEMP+CEN_ELP(2)
        END DO

    END IF
    !---!!!!!!!!!!!!!!!!!!---这之间带有对象特殊性---!!!!!!!!!!!!!!!!!!---!

    !------求解内法向量------!
    DO I=1,RSM,1
        CALL NORMALVECTOR_2D(BOUNDARY_GEOMETRICAL,X_SRFC(I),Y_SRFC(I),&
            N1X(I),N1Y(I))
        N1X(I)=-N1X(I)
        N1Y(I)=-N1Y(I)
        N1Z(I)=0.0D0
    END DO
    !------求解压力和偏导系数------!
    DO N=1,RSM,1

        !------求解点所处位置及压力------!
        XA=X_SRFC(N)
        YA=Y_SRFC(N)

        IF( XA>LEIN-CRITERIA .AND. XA<RIIN+CRITERIA )THEN
            IAU=IL+FLOOR( (XA-LEIN) / DX3  )
            IAV=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
            IAP=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"XA OUT OF INNER REGION"
            STOP
        END IF
        IF( YA>BOIN-CRITERIA .AND. YA<TOIN+CRITERIA )THEN
            JAU=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
            JAV=JB+FLOOR( (YA-BOIN) / DX3  )
            JAP=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
        ELSE
            WRITE(*,*)"YA OUT OF INNER REGION"
            STOP
        END IF

        CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA,YA,&
            U_SRFC(N),V_SRFC(N),P_SRFC(N))

        !------求解点刚体速度------!
        CALL VELOCITY_LB(BOUNDARY_KINETIC,XA,YA,U_SRFC(N),V_SRFC(N))

        !------获取附近符合要求的点------!
        CALL DETERMINE_IF_ELIGIBLE_U(IAU,JAU,ELIGIBILITY_U)
        CALL DETERMINE_IF_ELIGIBLE_V(IAV,JAV,ELIGIBILITY_V)

        !------求解距离------!
        DISTANCE_U=500.0D0
        DISTANCE_V=500.0D0
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(ELIGIBILITY_U(I,J))THEN
                    DISTANCE_U(I,J)=(X(IAU+I)-XA)**2.0D0&
                        +(YPU(JAU+J)-YA)**2.0D0
                END IF
                IF(ELIGIBILITY_V(I,J))THEN
                    DISTANCE_V(I,J)=(XPV(IAV+I)-XA)**2.0D0&
                        +(Y(JAV+J)-YA)**2.0D0
                END IF
            END DO
        END DO

        !------确定最近的两个点------!
        !刨去重合点
        DO J=-2,+2,1
            DO I=-2,+2,1
                IF(DISTANCE_U(I,J)<CRITERIA)THEN
                    DISTANCE_U(I,J)=500.0D0
                END IF
                IF(DISTANCE_V(I,J)<CRITERIA)THEN
                    DISTANCE_V(I,J)=500.0D0
                END IF
            END DO
        END DO

        DO WHILE (.TRUE.)

            IF ( MINVAL(DISTANCE_U)>=500.0D0 ) THEN
                WRITE(*,*)N,'u error: no eligible points'
                EXIT
            END IF

            CALL FIND_2SMALLEST_LOC(DISTANCE_U,INDEX_U)

            !检查坐标偏移量
            TEMP_1(1,2)=XA
            TEMP_1(2,2)=X(INDEX_U(1,1)+IAU)
            TEMP_1(3,2)=X(INDEX_U(2,1)+IAU)
            TEMP_1(1,3)=YA
            TEMP_1(2,3)=YPU(INDEX_U(1,2)+JAU)
            TEMP_1(3,3)=YPU(INDEX_U(2,2)+JAU)

            DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
            DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
            IF (DEVIATION_X<1000.0D0*CRITERIA .OR.&
                DEVIATION_Y<1000.0D0*CRITERIA ) THEN
                DISTANCE_U(INDEX_U(2,1),INDEX_U(2,2))=500.0D0
            ELSE
                EXIT
            END IF

        END DO

        INDEX_U(:,1)=INDEX_U(:,1)+IAU
        INDEX_U(:,2)=INDEX_U(:,2)+JAU

        DO WHILE (.TRUE.)

            IF ( MINVAL(DISTANCE_V)>=500.0D0 ) THEN
                WRITE(*,*)N,'v error: no eligible points'
                EXIT
            END IF

            CALL FIND_2SMALLEST_LOC(DISTANCE_V,INDEX_V)

            !检查坐标偏移量
            TEMP_1(1,2)=XA
            TEMP_1(2,2)=XPV(INDEX_V(1,1)+IAV)
            TEMP_1(3,2)=XPV(INDEX_V(2,1)+IAV)
            TEMP_1(1,3)=YA
            TEMP_1(2,3)=Y(INDEX_V(1,2)+JAV)
            TEMP_1(3,3)=Y(INDEX_V(2,2)+JAV)

            DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
            DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
            IF (DEVIATION_X<1000.0D0*CRITERIA .OR.&
                DEVIATION_Y<1000.0D0*CRITERIA ) THEN
                DISTANCE_V(INDEX_V(2,1),INDEX_V(2,2))=500.0D0
            ELSE
                EXIT
            END IF

        END DO

        INDEX_V(:,1)=INDEX_V(:,1)+IAV
        INDEX_V(:,2)=INDEX_V(:,2)+JAV

        !------求解偏导系数------!
        !U
        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(2,2)=X(INDEX_U(1,1))
        TEMP_1(3,2)=X(INDEX_U(2,1))
        TEMP_1(1,3)=YA
        TEMP_1(2,3)=YPU(INDEX_U(1,2))
        TEMP_1(3,3)=YPU(INDEX_U(2,2))
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'u error: negligible dy'
        END IF

        CALL M33INV (TEMP_1, TEMP_2, OK_FLAG)
        IF (.NOT. OK_FLAG) THEN
            WRITE(*,*) N,'u error: singular matrix'
        END IF

        TEMP_3(1,1)=U_SRFC(N)
        TEMP_3(2,1)=U(INDEX_U(1,1),INDEX_U(1,2))
        TEMP_3(3,1)=U(INDEX_U(2,1),INDEX_U(2,2))

        TEMP_4=MATMUL(TEMP_2,TEMP_3)

        AXX(N)=TEMP_4(2,1)
        AXY(N)=TEMP_4(3,1)

        !V
        TEMP_1(:,1)=1.0D0
        TEMP_1(1,2)=XA
        TEMP_1(2,2)=XPV(INDEX_V(1,1))
        TEMP_1(3,2)=XPV(INDEX_V(2,1))
        TEMP_1(1,3)=YA
        TEMP_1(2,3)=Y(INDEX_V(1,2))
        TEMP_1(3,3)=Y(INDEX_V(2,2))
        !检查坐标偏移量
        DEVIATION_X=MAXVAL(TEMP_1(:,2))-MINVAL(TEMP_1(:,2))
        IF (DEVIATION_X<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dx'
        END IF
        DEVIATION_Y=MAXVAL(TEMP_1(:,3))-MINVAL(TEMP_1(:,3))
        IF (DEVIATION_Y<1000.0D0*CRITERIA) THEN
            WRITE(*,*) N,'v error: negligible dy'
        END IF

        CALL M33INV (TEMP_1, TEMP_2, OK_FLAG)
        IF (.NOT. OK_FLAG) THEN
            WRITE(*,*) N,'v error: singular matrix'
        END IF

        TEMP_3(1,1)=V_SRFC(N)
        TEMP_3(2,1)=V(INDEX_V(1,1),INDEX_V(1,2))
        TEMP_3(3,1)=V(INDEX_V(2,1),INDEX_V(2,2))

        TEMP_4=MATMUL(TEMP_2,TEMP_3)

        AYX(N)=TEMP_4(2,1)
        AYY(N)=TEMP_4(3,1)

    END DO

    !------求解气动力------!
    DO N=1,RSM,1

        CXP=-P_SRFC(N)*N1X(N)
        CYP=-P_SRFC(N)*N1Y(N)
        CZP=0.0D0!-P_SRFC*N1Z(N)
        CXV=(2.0D0*AXX(N)*N1X(N)+(AXY(N)+AYX(N))*N1Y(N))/Re
        CYV=((AXY(N)+AYX(N))*N1X(N)+2.0D0*AYY(N)*N1Y(N))/Re
        CZV=0.0D0

        CXC_TOTAL=CXC_TOTAL+2.0D0*DS(N)*(CXP+CXV)
        CYC_TOTAL=CYC_TOTAL+2.0D0*DS(N)*(CYP+CYV)

    END DO

    !------物体所受气动力为反作用力------!
    CXC_TOTAL=-CXC_TOTAL
    CYC_TOTAL=-CYC_TOTAL

    !------输出系数------!
    WRITE(320+BOUNDARY_ID*10,"( I6,2(1X,F9.5))")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(321+BOUNDARY_ID*10,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    !------输出分布------!
    WRITE(CHAR_STEP,'(I6.6)') NSTEP
    REYNOLDS=IDNINT(Re)
    WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS
    WRITE(CHAR_ID,'(I2.2)') BOUNDARY_ID
    IF( MOD(NSTEP,2000)==0 )THEN
        OPEN(UNIT=10,FILE='SRFC_EXACT_OBJECT_'//TRIM(CHAR_ID)&
            //'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
        WRITE(10,*) 'TITLE="NONAME"'
        !WRITE(10,*) 'VARIABLES="X","Y","P","U","V","AXX","AXY","AXZ",'
        !WRITE(10,*) '"AYX","AYY","AYZ","AZX","AZY","AZZ"'
        WRITE(10,*) 'VARIABLES="X","Y","P","U","V","NX","NY","AXX","AXY",'
        WRITE(10,*) '"AYX","AYY"'
        WRITE(10,*) 'ZONE T="NONAME", I=',RSM,', J=',1,', F=POINT'
        DO N=1,RSM,1
            !WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
            !    P_SRFC(N),U_SRFC(N),V_SRFC(N),&
            !    AXX(N),AXY(N),AXZ(N),&
            !    AYX(N),AYY(N),AYZ(N),&
            !    AZX(N),AZY(N),AZZ(N)
            WRITE(10,*) X_SRFC(N),Y_SRFC(N),&
                P_SRFC(N),U_SRFC(N),V_SRFC(N),&
                N1X(N),N1Y(N),&
                AXX(N),AXY(N),&
                AYX(N),AYY(N)
        END DO
        CLOSE(10)
    END IF

    RETURN
    END SUBROUTINE

    !****************判断周围点是否为可视流体点（简易版）***************!
    SUBROUTINE DETERMINE_IF_ELIGIBLE_U(IAU,JAU,ELIGIBILITY_U)
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE

    INTEGER::IAU,JAU
    INTEGER::I,J
    LOGICAL::ELIGIBILITY_U(-2:2,-2:2)

    DO J=-2,+2,1
        DO I=-2,+2,1

            IF( TYPEUX(IAU+I,JAU+J)==-10 .AND. &
                TYPEUY(IAU+I,JAU+J)==-10 )THEN
                ELIGIBILITY_U(I,J)=.FALSE.
            ELSE
                ELIGIBILITY_U(I,J)=.TRUE.
            END IF

        END DO
    END DO

    RETURN
    END SUBROUTINE

    !****************判断周围点是否为可视流体点（简易版）***************!
    SUBROUTINE DETERMINE_IF_ELIGIBLE_V(IAV,JAV,ELIGIBILITY_V)
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE

    INTEGER::IAV,JAV
    INTEGER::I,J
    LOGICAL::ELIGIBILITY_V(-2:2,-2:2)

    DO J=-2,+2,1
        DO I=-2,+2,1

            IF( TYPEVX(IAV+I,JAV+J)==-10 .AND. &
                TYPEVY(IAV+I,JAV+J)==-10 )THEN
                ELIGIBILITY_V(I,J)=.FALSE.
            ELSE
                ELIGIBILITY_V(I,J)=.TRUE.
            END IF

        END DO
    END DO

    RETURN
    END SUBROUTINE

    !***************找出离边界点最近的两个可视流体点********************!
    SUBROUTINE FIND_2SMALLEST_LOC(DISTANCE,INDEX_1)
    IMPLICIT NONE

    INTEGER::INDEX_1(2,2)
    INTEGER::I,J
    REAL(KIND=8)::DISTANCE(-2:2,-2:2)
    REAL(KIND=8)::SMALLEST_1ST,SMALLEST_2ND

    SMALLEST_1ST=MAXVAL(DISTANCE)
    SMALLEST_2ND=MAXVAL(DISTANCE)
    INDEX_1(1,:)=MAXLOC(DISTANCE)
    INDEX_1(2,:)=MAXLOC(DISTANCE)

    DO J=-2,+2,1
        DO I=-2,+2,1

            IF(DISTANCE(I,J)<SMALLEST_1ST)THEN
                SMALLEST_2ND=SMALLEST_1ST
                SMALLEST_1ST=DISTANCE(I,J)
                INDEX_1(2,1)=INDEX_1(1,1)
                INDEX_1(2,2)=INDEX_1(1,2)
                INDEX_1(1,1)=I
                INDEX_1(1,2)=J
            ELSE IF(DISTANCE(I,J)<SMALLEST_2ND)THEN
                SMALLEST_2ND=DISTANCE(I,J)
                INDEX_1(2,1)=I
                INDEX_1(2,2)=J
            END IF

        END DO
    END DO

    RETURN
    END SUBROUTINE

    !*********************求解升力推力（表面力积分法）**********************!
    SUBROUTINE CAL_CLCT_2DCURVE
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !---------具有算法特殊性的部分---------!
    !运算区域
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CYLINDER_KINETIC(9)
    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !---------具有算法特殊性的部分---------!

    !网格密度
    REAL(KIND=8)::DN1!
    REAL(KIND=8),ALLOCATABLE::DN2(:)!
    !坐标值
    REAL(KIND=8),ALLOCATABLE::XBF(:,:),YBF(:,:)
    REAL(KIND=8),ALLOCATABLE::N1X(:),N1Y(:),N2X(:),N2Y(:)
    REAL(KIND=8)::XA,YA
    !网格数
    INTEGER::BFIM
    INTEGER::BFJM
    INTEGER::BFI,BFJ
    !插值涉及的绝对网格角点信息
    INTEGER::IAU,JAU!左下角脚标值
    INTEGER::IAV,JAV!左下角脚标值
    INTEGER::IAP,JAP!左下角脚标值

    !贴体速度压力场
    REAL(KIND=8),ALLOCATABLE::UBFN1(:,:),UBFN2(:,:),PBF(:,:)!曲线坐标系速度
    REAL(KIND=8),ALLOCATABLE::UBF(:,:),VBF(:,:)!绝对坐标系速度
    !坐标转换矩阵
    REAL(KIND=8)::MAT_ABS2CUR(2,2),MAT_CUR2ABS(2,2)

    !---------力系数---------!
    REAL(KIND=8),ALLOCATABLE::CFV1(:),CFV2(:),CFP(:)
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!计算域绝对坐标系下的力
    REAL(KIND=8)::CXT_TOTAL,CYT_TOTAL!真实世界坐标系下的力
    REAL(KIND=8)::CXF_TOTAL,CYF_TOTAL!扑翼坐标系下的力

    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------边界1-----------------------------------------!
    CYLINDER_KINETIC=QUADRIC_KINETIC_1
    CEN_CRC(1)=QUADRIC_KINETIC_1(4)
    CEN_CRC(2)=QUADRIC_KINETIC_1(5)
    COX2=QUADRIC_GEOMETRICAL_1(7)
    COY2=QUADRIC_GEOMETRICAL_1(8)
    COXY=QUADRIC_GEOMETRICAL_1(9)
    COX=QUADRIC_GEOMETRICAL_1(10)
    COY=QUADRIC_GEOMETRICAL_1(11)
    COM=QUADRIC_GEOMETRICAL_1(13)
    !认为COXY==0，COX2==COY2
    RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

    BFJM=3
    BFIM=IDNINT( 2.0D0*RADIUS/DX3 )

    ALLOCATE( DN2(BFJM) )

    DN1=DSQRT(2.0D0)*DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "升阻力系数格点数1：",BFIM,BFJM

    ALLOCATE( XBF(BFIM,BFJM),YBF(BFIM,BFJM) )
    ALLOCATE( N1X(BFIM),N1Y(BFIM),N2X(BFIM),N2Y(BFIM) )
    ALLOCATE( UBFN1(BFIM,BFJM),UBFN2(BFIM,BFJM),PBF(BFIM,BFJM) )
    ALLOCATE( UBF(BFIM,BFJM),VBF(BFIM,BFJM) )
    ALLOCATE( CFV1(BFIM),CFV2(BFIM),CFP(BFIM) )

    !------绘制贴体网格------!
    DO BFJ=1,BFJM,1
        DO BFI=1,BFIM,1
            !从X轴逆时针旋转
            XBF(BFI,BFJ)=CEN_CRC(1)+(RADIUS+DBLE(BFJM-BFJ)*DN1)*DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            YBF(BFI,BFJ)=CEN_CRC(2)+(RADIUS+DBLE(BFJM-BFJ)*DN1)*DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N1X(BFI)=-DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N1Y(BFI)=-DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N2X(BFI)= DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N2Y(BFI)=-DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
        END DO
    END DO

    !------初始化------!
    UBFN1=0.0D0
    UBFN2=0.0D0
    PBF=0.0D0
    UBF=0.0D0
    VBF=0.0D0
    CFV1=0.0D0
    CFV2=0.0D0
    CFP=0.0D0

    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0
    CXT_TOTAL=0.0D0
    CYT_TOTAL=0.0D0
    CXF_TOTAL=0.0D0
    CYF_TOTAL=0.0D0

    !------外层点直接插值------!
    DO BFJ=1,BFJM-1,1
        DO BFI=1,BFIM,1
            XA=XBF(BFI,BFJ)
            YA=YBF(BFI,BFJ)

            IF( XA>LEIN-CRITERIA .AND. XA<RIIN+CRITERIA )THEN
                IAU=IL+FLOOR( (XA-LEIN) / DX3  )
                IAV=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
                IAP=IL+FLOOR( (XA-LEIN) / DX3 - 0.5D0  )
            ELSE
                WRITE(*,*)"XA ERROR"
                STOP
            END IF
            IF( YA>BOIN-CRITERIA .AND. YA<TOIN+CRITERIA )THEN
                JAU=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
                JAV=JB+FLOOR( (YA-BOIN) / DX3  )
                JAP=JB+FLOOR( (YA-BOIN) / DX3 - 0.5D0  )
            ELSE
                WRITE(*,*)"YA ERROR"
                STOP
            END IF

            CALL BILINEAR_INTERPOLATION(IAU,JAU,IAV,JAV,IAP,JAP,XA,YA,UBF(BFI,BFJ),VBF(BFI,BFJ),PBF(BFI,BFJ))

        END DO
    END DO

    !------翅膀层赋予动边界运动速度------!
    BFJ=BFJM
    DO BFI=1,BFIM,1
        XA=XBF(BFI,BFJ)
        YA=YBF(BFI,BFJ)
        CALL VELOCITY_LB(CYLINDER_KINETIC,XA,YA,UBF(BFI,BFJ),VBF(BFI,BFJ))
    END DO

    !------速度场坐标变换------!
    DO BFI=1,BFIM,1

        MAT_CUR2ABS(1,1)=N1X(BFI)
        MAT_CUR2ABS(2,1)=N1Y(BFI)
        MAT_CUR2ABS(1,2)=N2X(BFI)
        MAT_CUR2ABS(2,2)=N2Y(BFI)
        MAT_ABS2CUR=TRANSPOSE(MAT_CUR2ABS)

        DO BFJ=1,BFJM,1

            CALL CRDNT_TRANSFORM(MAT_ABS2CUR,UBF(BFI,BFJ),VBF(BFI,BFJ),UBFN1(BFI,BFJ),UBFN2(BFI,BFJ))

        END DO
    END DO

    !------求解气动力------!
    BFJ=BFJM-1
    BFI=1
    CFV1(BFI)=2.0D0/Re*2.0D0*(UBFN1(BFI,BFJ+1)-UBFN1(BFI,BFJ-1))/(2.0D0*DN1)
    CFV2(BFI)=2.0D0/Re*( (UBFN1(2,BFJ)-UBFN1(BFIM,BFJ))/(2.0D0*DN2(BFJ)) &
        & + (UBFN2(BFI,BFJ+1)-UBFN2(BFI,BFJ-1))/(2.0D0*DN1) )
    CFP(BFI) =-2.0D0*PBF(BFI,BFJ)

    CXC_TOTAL=CXC_TOTAL+(CFV1(BFI)*N1X(BFI)+CFV2(BFI)*N2X(BFI)+CFP(BFI)*N1X(BFI))*DN2(BFJ)
    CYC_TOTAL=CYC_TOTAL+(CFV1(BFI)*N1Y(BFI)+CFV2(BFI)*N2Y(BFI)+CFP(BFI)*N1Y(BFI))*DN2(BFJ)

    BFI=BFIM
    CFV1(BFI)=2.0D0/Re*2.0D0*(UBFN1(BFI,BFJ+1)-UBFN1(BFI,BFJ-1))/(2.0D0*DN1)
    CFV2(BFI)=2.0D0/Re*( (UBFN1(1,BFJ)-UBFN1(BFIM-1,BFJ))/(2.0D0*DN2(BFJ)) &
        & + (UBFN2(BFI,BFJ+1)-UBFN2(BFI,BFJ-1))/(2.0D0*DN1) )
    CFP(BFI) =-2.0D0*PBF(BFI,BFJ)

    CXC_TOTAL=CXC_TOTAL+(CFV1(BFI)*N1X(BFI)+CFV2(BFI)*N2X(BFI)+CFP(BFI)*N1X(BFI))*DN2(BFJ)
    CYC_TOTAL=CYC_TOTAL+(CFV1(BFI)*N1Y(BFI)+CFV2(BFI)*N2Y(BFI)+CFP(BFI)*N1Y(BFI))*DN2(BFJ)
    DO BFI=2,BFIM-1,1
        CFV1(BFI)=2.0D0/Re*2.0D0*(UBFN1(BFI,BFJ+1)-UBFN1(BFI,BFJ-1))/(2.0D0*DN1)
        CFV2(BFI)=2.0D0/Re*( (UBFN1(BFI+1,BFJ)-UBFN1(BFI-1,BFJ))/(2.0D0*DN2(BFJ)) &
            & + (UBFN2(BFI,BFJ+1)-UBFN2(BFI,BFJ-1))/(2.0D0*DN1) )
        CFP(BFI) =-2.0D0*PBF(BFI,BFJ)

        !CXC_TOTAL=CXC_TOTAL+(CFV1(BFI)*N1X(BFI)+CFV2(BFI)*N2X(BFI)+CFP(BFI)*N1X(BFI))*DN2(BFJM)
        !CYC_TOTAL=CYC_TOTAL+(CFV1(BFI)*N1Y(BFI)+CFV2(BFI)*N2Y(BFI)+CFP(BFI)*N1Y(BFI))*DN2(BFJM)

        CXC_TOTAL=CXC_TOTAL+(CFV1(BFI)*N1X(BFI)+CFV2(BFI)*N2X(BFI)+CFP(BFI)*N1X(BFI))*DN2(BFJ)
        CYC_TOTAL=CYC_TOTAL+(CFV1(BFI)*N1Y(BFI)+CFV2(BFI)*N2Y(BFI)+CFP(BFI)*N1Y(BFI))*DN2(BFJ)
    END DO

    CXC_TOTAL=-CXC_TOTAL
    CYC_TOTAL=-CYC_TOTAL

    WRITE(322,"( I6,2(1X,F9.5))")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(323,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE

    !*********************求解升力推力（控制体方法）**********************!
    SUBROUTINE CAL_CLCT_CONTROL_VOLUME
    USE DECLARATION
    IMPLICIT NONE

    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !控制体区域
    REAL(KIND=8)::CV_LEFT,CV_RIGH,CV_BOTT,CV_TOPP

    !寻址
    REAL(KIND=8),ALLOCATABLE::DISTANCE_X(:),DISTANCE_Y(:)
    INTEGER::INDEX_XL,INDEX_XR
    INTEGER::INDEX_YB,INDEX_YT

    !界面信息
    !REAL(KIND=8),ALLOCATABLE::U_LEFT_BOUNDARY(:),U_RIGH_BOUNDARY(:)
    !REAL(KIND=8),ALLOCATABLE::P_LEFT_BOUNDARY(:),P_RIGH_BOUNDARY(:)
    !REAL(KIND=8),ALLOCATABLE::BOUNDARY_AREA(:)

    !------力系数------!
    !计算域绝对坐标系下的力
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL

    !非定常
    REAL(KIND=8)::CFTX,CFTY

    !对流
    REAL(KIND=8)::CFCLX,CFCLY
    REAL(KIND=8)::CFCRX,CFCRY
    REAL(KIND=8)::CFCBX,CFCBY
    REAL(KIND=8)::CFCTX,CFCTY

    !压力
    REAL(KIND=8)::CFPL
    REAL(KIND=8)::CFPR
    REAL(KIND=8)::CFPB
    REAL(KIND=8)::CFPT

    !粘性应力
    REAL(KIND=8)::CFVLX,CFVLY
    REAL(KIND=8)::CFVRX,CFVRY
    REAL(KIND=8)::CFVBX,CFVBY
    REAL(KIND=8)::CFVTX,CFVTY

    !合力
    REAL(KIND=8)::CFCX,CFCY
    REAL(KIND=8)::CFPX,CFPY
    REAL(KIND=8)::CFVX,CFVY


    ALLOCATE( DISTANCE_X(IM),DISTANCE_Y(JM) )

    IF(IB_LOCOMOTION==-1)THEN
        CV_LEFT=-0.6D0
        CV_RIGH= 0.6D0
        CV_BOTT=-0.6D0
        CV_TOPP= 0.6D0
    ELSE
        CV_LEFT=LEIN+2.0D0*DX3
        CV_RIGH=RIIN-2.0D0*DX3
        CV_BOTT=BOIN+2.0D0*DX3
        CV_TOPP=TOIN-2.0D0*DX3
    END IF

    !求解最近坐标
    DISTANCE_X=X(:)-CV_LEFT
    INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    DISTANCE_X=X(:)-CV_RIGH
    INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    DISTANCE_Y=Y(:)-CV_BOTT
    INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    IF( DABS(Y(INDEX_YB)-CV_BOTT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    DISTANCE_Y=Y(:)-CV_TOPP
    INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    IF( DABS(Y(INDEX_YT)-CV_TOPP)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"

    !力初始化
    CXC_TOTAL=0.0D0
    CYC_TOTAL=0.0D0

    CFTX=0.0D0
    CFTY=0.0D0

    CFCLX=0.0D0
    CFCLY=0.0D0
    CFCRX=0.0D0
    CFCRY=0.0D0
    CFCBX=0.0D0
    CFCBY=0.0D0
    CFCTX=0.0D0
    CFCTY=0.0D0

    CFPL=0.0D0
    CFPR=0.0D0
    CFPB=0.0D0
    CFPT=0.0D0

    CFVLX=0.0D0
    CFVLY=0.0D0
    CFVRX=0.0D0
    CFVRY=0.0D0
    CFVBX=0.0D0
    CFVBY=0.0D0
    CFVTX=0.0D0
    CFVTY=0.0D0

    CFCX=0.0D0
    CFCY=0.0D0
    CFPX=0.0D0
    CFPY=0.0D0
    CFVX=0.0D0
    CFVY=0.0D0

    !求解各项面力
    I=INDEX_XL
    DO J=INDEX_YB,INDEX_YT-1,1
        CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
        CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
        CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
        CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
        CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
            & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    END DO
    I=INDEX_XR
    DO J=INDEX_YB,INDEX_YT-1,1
        CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
        CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
        CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
        CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
        CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
            & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    END DO
    J=INDEX_YB
    DO I=INDEX_XL,INDEX_XR-1,1
        CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
        CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
        CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
        CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
            & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
        CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    END DO
    J=INDEX_YT
    DO I=INDEX_XL,INDEX_XR-1,1
        CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
        CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
        CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
        CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
            & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
        CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    END DO
    !求解非定常
    DO J=INDEX_YB,INDEX_YT-1,1
        DO I=INDEX_XL,INDEX_XR-1,1
            CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
            CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
        END DO
    END DO

    !合力
    CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    CFPX=-CFPL+CFPR
    CFPY=-CFPB+CFPT
    CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    CFVY=-CFVLY+CFVRY-CFVBY+CFVTY

    !计算域绝对坐标系下的力
    CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)

    WRITE(320,"( I6,2(1X,F9.5))")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(321,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE