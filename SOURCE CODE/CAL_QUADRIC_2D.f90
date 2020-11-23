    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !*************************************求解蜻蜓翅膀某一时刻下欧拉角等角度信息，确定转换矩阵系数******************************************!
    SUBROUTINE CAL_QUADRIC_2D
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE

    QUADRIC_GEOMETRICAL_N1=QUADRIC_GEOMETRICAL_1
    QUADRIC_GEOMETRICAL_N2=QUADRIC_GEOMETRICAL_2
    QUADRIC_KINETIC_N1=QUADRIC_KINETIC_1
    QUADRIC_KINETIC_N2=QUADRIC_KINETIC_2

    !--------------------边界1---------------------!
    IF( BOUNDARY_EXISTENCE_1==1 )THEN
        !----------转动中心位置；平转动速度；坐标转换所需-----------!
        IF(IB_LOCOMOTION==0)THEN
            PHASE_DIFFERENCE=0.0D0!静置时无相位差值
            PHASE_INITIATION=90.0D0!0.0D0!
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=-0.6D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0!转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_STATIC
        ELSE IF(IB_LOCOMOTION==-1)THEN
            PHASE_DIFFERENCE=0.0D0!静置时无相位差值
            PHASE_INITIATION=0.0D0!0.0D0!
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0!转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_STATIC
        ELSE IF(IB_LOCOMOTION==-4)THEN
            PHASE_DIFFERENCE=0.0D0!静置时无相位差值
            PHASE_INITIATION=0.0D0!0.0D0!
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0!转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_STATIC
        ELSE IF(IB_LOCOMOTION==1)THEN!模拟1请确认符合模拟目标
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=-0.6D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC(T)
        ELSE IF(IB_LOCOMOTION==11)THEN!模拟1请确认符合模拟目标
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC_WANG(T)
        ELSE IF(IB_LOCOMOTION==2)THEN
            PHASE_DIFFERENCE=-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=-0.6D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            INTERMITTENT_INITIATION=5.0D0
            CALL POSE_VELO_QUADRIC_2D_INTERMITTENT_FORE(T)
        ELSE IF(IB_LOCOMOTION==3)THEN
            CEN_TRANSLATION(1)=3.0D0+0.5D0*1.0D0*0.01D0
            CEN_TRANSLATION(2)=0.0D0
            CALL POSE_VELO_QUADRIC_2D_IMPULSIVE_START
        ELSE IF(IB_LOCOMOTION==41)THEN!暂时用于lock-in研究
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0  !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0  !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_STATIC_FOR_HEAVING_PLUNGING
        ELSE IF(IB_LOCOMOTION==42)THEN!暂时用于lock-in研究
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0 !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC_HEAVING_PLUNGING(T)
        ELSE IF(IB_LOCOMOTION==4)THEN
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0!转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC_HEAVING_PLUNGING(T)
        ELSE IF(IB_LOCOMOTION==5)THEN
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=-90.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0!转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0!转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC_HEAVING_PLUNGING(T)
        ELSE IF(IB_LOCOMOTION==6)THEN
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=-90.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0 !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC_OSCILLATING_X(T)
        ELSE IF(IB_LOCOMOTION==7)THEN
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=0.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0    !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0    !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PURE_ROTATING(T)
        ELSE IF(IB_LOCOMOTION==8)THEN
            PHASE_DIFFERENCE=0.0D0!-75.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=0.0D0
            CEN_DEVIATION(1)=0.5D0-0.5D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.0D0 !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0 !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PURE_ROTATING_STEADY(T)
        END IF

        WRITE(180,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,PHIW/PI*180.0D0,PSIW/PI*180.0D0
        IF(IB_LOCOMOTION==1 .OR. IB_LOCOMOTION==2)THEN
            WRITE(190,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,VELO_TRAN_R(2),VELO_ANGL
        ELSE
            WRITE(190,"( I6,(1X,F9.5),3(1X,F14.10))") NSTEP,TAU/TAUC,VELO_TRAN_A,VELO_ANGL
        END IF
        WRITE(220,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,CEN(1),CEN(2)

        !----------确定二次曲线表达式-----------!
        CALL CAL_QUADRIC_COEFFICIENT

        QUADRIC_KINETIC_1(1)=VELO_TRAN_A(1)
        QUADRIC_KINETIC_1(2)=VELO_TRAN_A(2)
        QUADRIC_KINETIC_1(4)=CEN(1)
        QUADRIC_KINETIC_1(5)=CEN(2)
        QUADRIC_KINETIC_1(7)=VELO_ANGL

        QUADRIC_GEOMETRICAL_1(1)=T11
        QUADRIC_GEOMETRICAL_1(2)=T12
        QUADRIC_GEOMETRICAL_1(4)=T21
        QUADRIC_GEOMETRICAL_1(5)=T22

        QUADRIC_GEOMETRICAL_1(7)=COX2
        QUADRIC_GEOMETRICAL_1(8)=COY2
        QUADRIC_GEOMETRICAL_1(9)=COXY
        QUADRIC_GEOMETRICAL_1(10)=COX
        QUADRIC_GEOMETRICAL_1(11)=COY
        QUADRIC_GEOMETRICAL_1(13)=COM

        QUADRIC_GEOMETRICAL_1(15)=XMAXT
        QUADRIC_GEOMETRICAL_1(16)=XMINT
        QUADRIC_GEOMETRICAL_1(17)=YMAXT
        QUADRIC_GEOMETRICAL_1(18)=YMINT

        QUADRIC_GEOMETRICAL_1(21)=XMAXB
        QUADRIC_GEOMETRICAL_1(22)=XMINB
        QUADRIC_GEOMETRICAL_1(23)=YMAXB
        QUADRIC_GEOMETRICAL_1(24)=YMINB

        QUADRIC_GEOMETRICAL_1(29)=RYTOP
        QUADRIC_GEOMETRICAL_1(30)=RYBOT
        QUADRIC_GEOMETRICAL_1(31)=RXTRL
        QUADRIC_GEOMETRICAL_1(32)=RXLED

        QUADRIC_GEOMETRICAL_1(33)=LAXIS
        QUADRIC_GEOMETRICAL_1(34)=SAXIS
        QUADRIC_GEOMETRICAL_1(35)=CEN_DEVIATION(1)
        QUADRIC_GEOMETRICAL_1(36)=CEN_DEVIATION(2)

    END IF

    !--------------------边界2---------------------!
    IF( BOUNDARY_EXISTENCE_2==1 )THEN
        !----------转动中心位置；平转动速度；坐标转换所需-----------!
        IF(IB_LOCOMOTION==0)THEN
            PHASE_DIFFERENCE=0.0D0!静置时无相位差值
            PHASE_INITIATION=90.0D0!0.0D0!
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.6D0    !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0    !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_STATIC
        ELSE IF(IB_LOCOMOTION==1)THEN!模拟1请确认符合模拟目标
            PHASE_DIFFERENCE=0.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.6D0      !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0      !转动中心绝对坐标系下振荡中心Y
            CALL POSE_VELO_QUADRIC_2D_PERIODIC(T)
        ELSE IF(IB_LOCOMOTION==2)THEN
            PHASE_DIFFERENCE=0.0D0!此时以后翼为基准，相应地前翼有一个负的相位差值
            PHASE_INITIATION=90.0D0
            CEN_DEVIATION(1)=0.5D0-0.3D0!转动中心在弦向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_DEVIATION(2)=0.5D0-0.5D0!转动中心在拍动向相对二次图形几何图形中心偏移量，更改第二个数字为转动中心相对位置
            CEN_TRANSLATION(1)=0.6D0        !转动中心绝对坐标系下振荡中心X
            CEN_TRANSLATION(2)=0.0D0        !转动中心绝对坐标系下振荡中心Y
            INTERMITTENT_INITIATION=5.0D0
            CALL POSE_VELO_QUADRIC_2D_INTERMITTENT_HIND(T)
        END IF

        WRITE(200,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,PHIW/PI*180.0D0,PSIW/PI*180.0D0
        WRITE(210,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,VELO_TRAN_R(2),VELO_ANGL
        WRITE(230,"( I6,3(1X,F14.10))") NSTEP,TAU/TAUC,CEN(1),CEN(2)

        !----------确定二次曲线表达式-----------!
        CALL CAL_QUADRIC_COEFFICIENT

        QUADRIC_KINETIC_2(1)=VELO_TRAN_A(1)
        QUADRIC_KINETIC_2(2)=VELO_TRAN_A(2)
        QUADRIC_KINETIC_2(4)=CEN(1)
        QUADRIC_KINETIC_2(5)=CEN(2)
        QUADRIC_KINETIC_2(7)=VELO_ANGL

        QUADRIC_GEOMETRICAL_2(1)=T11
        QUADRIC_GEOMETRICAL_2(2)=T12
        QUADRIC_GEOMETRICAL_2(4)=T21
        QUADRIC_GEOMETRICAL_2(5)=T22

        QUADRIC_GEOMETRICAL_2(7)=COX2
        QUADRIC_GEOMETRICAL_2(8)=COY2
        QUADRIC_GEOMETRICAL_2(9)=COXY
        QUADRIC_GEOMETRICAL_2(10)=COX
        QUADRIC_GEOMETRICAL_2(11)=COY
        QUADRIC_GEOMETRICAL_2(13)=COM

        QUADRIC_GEOMETRICAL_2(15)=XMAXT
        QUADRIC_GEOMETRICAL_2(16)=XMINT
        QUADRIC_GEOMETRICAL_2(17)=YMAXT
        QUADRIC_GEOMETRICAL_2(18)=YMINT

        QUADRIC_GEOMETRICAL_2(21)=XMAXB
        QUADRIC_GEOMETRICAL_2(22)=XMINB
        QUADRIC_GEOMETRICAL_2(23)=YMAXB
        QUADRIC_GEOMETRICAL_2(24)=YMINB

        QUADRIC_GEOMETRICAL_2(29)=RYTOP
        QUADRIC_GEOMETRICAL_2(30)=RYBOT
        QUADRIC_GEOMETRICAL_2(31)=RXTRL
        QUADRIC_GEOMETRICAL_2(32)=RXLED

        QUADRIC_GEOMETRICAL_2(33)=LAXIS
        QUADRIC_GEOMETRICAL_2(34)=SAXIS
        QUADRIC_GEOMETRICAL_2(35)=CEN_DEVIATION(1)
        QUADRIC_GEOMETRICAL_2(36)=CEN_DEVIATION(2)

    END IF

    RETURN
    END SUBROUTINE


    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************转动中心位置；平转动速度；坐标转换所需（前翅间歇性拍动）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_INTERMITTENT_FORE(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME


    !--------------------本子函数根据不同扑翼规律需改变下方---------------------!
    !--------------------周期性拍动基本参数---------------------!
    TAUC=8.296510765D0!周期时长
    !旋转时间/刻
    TAU_R1=0.0625D0*TAUC
    TAU_R2=0.5625D0*TAUC!两个翻转开始的时刻
    DTAUR=0.25D0*TAUC!单次翻转时长
    !几何信息
    PHIM=32.5D0/180.0D0*PI!拍动振幅
    PSIM=50.0D0/180.0D0*PI!翻转振幅

    SPAN=4.556701031D0*0.8D0!展长
    !--------------------本子函数根据不同扑翼规律需改变上方---------------------!

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------间歇性拍动基本参数---------------------!
    !间歇性特征时长
    DTAUH=5.0D0/48.0D0*TAUC!HIGH STROKE
    DTAUI=1.0D0/2.0D0*TAUC!0.0D0!INTERMISSION
    DTAUL=7.0D0/48.0D0*TAUC!LOW STROKE
    !间歇性特征时刻
    TSTART=INTERMITTENT_INITIATION*TAUC-PHASE_DIFFERENCE/360.0D0*TAUC-TAU_INIT!0.0
    TSTOP=TSTART+DTAUL
    TRESTART=TSTOP+DTAUI
    TEND=TRESTART+DTAUH
    !几何信息
    PHII= PHIM*DSIN(2.0D0*PI*DTAUL/TAUC)!间隙拍动角
    PSII=-PSIM*DCOS(PI*(DTAUL-TAU_R1)/DTAUR)!间隙翻转角

    !--------------------周期性拍动分段---------------------!
    IF(TIME<TSTART)THEN
        CALL POSE_VELO_QUADRIC_2D_PERIODIC(TIME)
        RETURN
    ELSE IF(TIME>TEND)THEN
        CALL POSE_VELO_QUADRIC_2D_PERIODIC(TIME-TEND+TSTART+0.25D0*TAUC)
        RETURN
    END IF

    !--------------------间歇性拍动周期内时刻TAU---------------------!
    IF(TIME<=TSTOP)THEN
        TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)
    ELSE IF(TIME>=TRESTART)THEN
        TAU=DMOD(TIME-TRESTART+TSTOP+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)
    ELSE
        TAU=0.00D0
    END IF

    !--------------------根据TAU/TAUI确定各角度大小（间歇性时附加）---------------------!
    PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
    !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
    !ϕw拍动角
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        TAUI=TIME-TSTART
        PHIW=PHIM*( (1.0D0-TAUI/DTAUL)*DSIN(2.0D0*PI*TAUI/TAUC)+TAUI/DTAUL*DSIN(2.0D0*PI*DTAUL/TAUC)*DSIN(2.0D0*PI*TAUI/4.0D0/DTAUL) )
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        PHIW=PHII
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        TAUI=(TIME-TRESTART)
        PHIW=PHII+(PHIM-PHII)/2.0D0*(1.0D0+DSIN(2.0D0*PI*TAUI/2.0D0/DTAUH-0.5D0*PI))
    END IF
    !ψw翻转角
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        IF(TAU>=0.0D0 .AND. TAU<=TAU_R1)THEN
            PSIW=-PSIM
        ELSE IF(TAU>TAU_R1 .AND. TAU<TSTOP-TSTART)THEN
            PSIW=-PSIM*DCOS(PI*(TAU-TAU_R1)/DTAUR)
        END IF
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        PSIW=-PSIM*DCOS(PI*( (TSTOP-TSTART)-TAU_R1)/DTAUR)
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        PSIW=-PSIM*DCOS(PI*(TAU-(TRESTART-TSTOP)-TAU_R1)/DTAUR)
    END IF

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（间歇性时附加）---------------------!
    CEN_P(1)=CEN_TRANSLATION(1)+0.0D0!-0.8D0
    CEN_P(2)=CEN_TRANSLATION(2)+0.0D0+PHIW*SPAN
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    !间歇性平动速度
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        TAUI=TIME-TSTART
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)= SPAN*( PHIM*( -DSIN(2.0D0*PI*TAUI/TAUC)/DTAUL + (1.0D0-TAUI/DTAUL)*2.0D0*PI/TAUC*DCOS(2.0D0*PI*TAUI/TAUC) ) + PHII/DTAUL*DSIN(2.0D0*PI*TAUI/4.0D0/DTAUL) + PHII*TAUI/DTAUL*2.0D0*PI/4.0D0/DTAUL*DCOS(2.0D0*PI*TAUI/4.0D0/DTAUL) )
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        TAUI=TIME-TRESTART
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=SPAN*(PHIM-PHII)/2.0D0*2.0D0*PI/2.0D0/DTAUH*DCOS(2.0D0*PI*TAUI/2.0D0/DTAUH-0.5D0*PI)
    END IF
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !间歇性转动速度
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        IF(TAU>=0 .AND. TAU<=TAU_R1)THEN
            VELO_ANGL=0.0D0
        ELSE IF(TAU>TAU_R1 .AND. TAU<TSTOP-TSTART)THEN
            VELO_ANGL=PI*PSIM/DTAUR*DSIN(PI*(TAU-TAU_R1)/DTAUR)
        END IF
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        VELO_ANGL=0.0D0
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        VELO_ANGL=PI*PSIM/DTAUR*DSIN(PI*(TAU-(TRESTART-TSTOP)-TAU_R1)/DTAUR)
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（后翅间歇性拍动）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_INTERMITTENT_HIND(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !--------------------本子函数根据不同扑翼规律需改变下方---------------------!
    !--------------------周期性拍动基本参数---------------------!
    TAUC=8.296510765D0!周期时长
    !旋转时间/刻
    TAU_R1=0.0625D0*TAUC
    TAU_R2=0.5625D0*TAUC!两个翻转开始的时刻
    DTAUR=0.25D0*TAUC!单次翻转时长
    !几何信息
    PHIM=32.5D0/180.0D0*PI!拍动振幅
    PSIM=50.0D0/180.0D0*PI!翻转振幅

    SPAN=4.556701031D0*0.8D0!展长
    !--------------------本子函数根据不同扑翼规律需改变上方---------------------!

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------间歇性拍动基本参数---------------------!
    !间歇性特征时长
    DTAUH=5.0D0/48.0D0*TAUC!HIGH STROKE
    DTAUI=1.0D0/2.0D0*TAUC!0.0D0!INTERMISSION
    DTAUL=7.0D0/48.0D0*TAUC!LOW STROKE
    !间歇性特征时刻
    TSTART=(INTERMITTENT_INITIATION+0.25D0)*TAUC-PHASE_DIFFERENCE/360.0D0*TAUC-TAU_INIT!0.0
    TSTOP=TSTART+DTAUH
    TRESTART=TSTOP+DTAUI
    TEND=TRESTART+DTAUL
    !几何信息
    PHII= PHIM*DSIN(2.0D0*PI*DTAUL/TAUC)!间隙拍动角
    PSII=-PSIM*DCOS(PI*(DTAUL-TAU_R1)/DTAUR)!间隙翻转角

    !--------------------周期性拍动分段---------------------!
    IF(TIME<TSTART)THEN
        CALL POSE_VELO_QUADRIC_2D_PERIODIC(TIME)
        RETURN
    ELSE IF(TIME>TEND)THEN
        CALL POSE_VELO_QUADRIC_2D_PERIODIC(TIME-TEND+TSTART+0.25D0*TAUC)
        RETURN
    END IF

    !--------------------间歇性拍动周期内时刻TAU---------------------!
    IF(TIME<=TSTOP)THEN
        TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)
    ELSE IF(TIME>=TRESTART)THEN
        TAU=DMOD(TIME-TRESTART+TSTOP+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)
    ELSE
        TAU=0.00D0
    END IF

    !--------------------根据TAU/TAUI确定各角度大小（间歇性时附加）---------------------!
    PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
    !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
    !ϕw拍动角
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        TAUI=(TIME-TSTART)/TAUC
        PHIW=PHII+(PHIM-PHII)/2.0D0*(1.0D0+DSIN(2.0D0*PI*TAUI*24.0D0/5.0D0+0.5D0*PI))
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        PHIW=PHII
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        TAUI=(TEND-TIME)/TAUC
        PHIW=PHIM*( (1.0D0-TAUI*48.0D0/7.0D0)*DSIN(2.0D0*PI*TAUI)+TAUI*48.0D0/7.0D0*DSIN(PI*7.0D0/24.0D0)*DSIN(2.0D0*PI*TAUI/7.0D0*12.0D0) )
    END IF
    !ψw翻转角
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        IF(TAU<=TAU_R1+DTAUR)THEN!TAU>=0.25D0*TAUC .AND.
            PSIW=-PSIM*DCOS(PI*(TAU-TAU_R1)/DTAUR)
        ELSE IF(TAU>TAU_R1+DTAUR)THEN! .AND. TAU<TSTOP-TSTART
            PSIW=PSIM
        END IF
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        PSIW=PSIM
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        PSIW=PSIM
    END IF

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（间歇性时附加）---------------------!
    CEN_P(1)=CEN_TRANSLATION(1)+0.0D0!-0.8D0
    CEN_P(2)=CEN_TRANSLATION(2)+0.0D0+PHIW*SPAN
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    !间歇性平动速度
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        TAUI=TIME-TSTART
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=SPAN*(PHIM-PHII)/2.0D0*2.0D0*PI/2.0D0/DTAUH*DCOS(2.0D0*PI*TAUI/2.0D0/DTAUH+0.5D0*PI)
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        TAUI=TEND-TIME
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=SPAN*( PHIM*( DSIN(2.0D0*PI*TAUI/TAUC)/DTAUL - (1.0D0-TAUI/DTAUL)*2.0D0*PI/TAUC*DCOS(2.0D0*PI*TAUI/TAUC) ) - PHII/DTAUL*DSIN(2.0D0*PI*TAUI/4.0D0/DTAUL) - PHII*TAUI/DTAUL*2.0D0*PI/4.0D0/DTAUL*DCOS(2.0D0*PI*TAUI/4.0D0/DTAUL) )
    END IF
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !间歇性转动速度
    IF( TIME>=TSTART .AND. TIME<=TSTOP )THEN
        IF(TAU<=TAU_R1+DTAUR)THEN!TAU>=0.25D0*TAUC .AND.
            VELO_ANGL=PI*PSIM/DTAUR*DSIN(PI*(TAU-TAU_R1)/DTAUR)
        ELSE IF(TAU>TAU_R1+DTAUR)THEN! .AND. TAU<TSTOP-TSTART
            VELO_ANGL=0.0D0
        END IF
    ELSE IF( TIME>TSTOP .AND. TIME<=TRESTART )THEN
        VELO_ANGL=0.0D0
    ELSE IF( TIME>TRESTART .AND. TIME<=TEND )THEN
        VELO_ANGL=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（周期性拍动）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PERIODIC(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !模拟1请确认符合模拟目标
    !--------------------本子函数根据不同扑翼规律需改变下方---------------------!
    !--------------------周期性拍动基本参数---------------------!
    TAUC=8.623671834D0!7.725054831D0!周期时长
    !旋转时间/刻
    TAU_R1=0.125D0*TAUC!0.0625D0*TAUC
    TAU_R2=0.625D0*TAUC!0.5625D0*TAUC!两个翻转开始的时刻
    DTAUR=0.25D0*TAUC!单次翻转时长
    !几何信息
    PHIM=45.0D0/180.0D0*PI!45.0D0/180.0D0*PI!拍动振幅
    PSIM=45.0D0/180.0D0*PI!翻转振幅
    PSI0=-15.0D0/180.0D0*PI!初始翻转角

    SPAN=2.745D0!展长
    !--------------------本子函数根据不同扑翼规律需改变上方---------------------!

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据TAU确定各角度大小（周期性）---------------------!
    IF( (TIME+PHASE_DIFFERENCE/360.0D0*TAUC)<-CRITERIA)THEN
        PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
        !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
        THETAW=0.0D0!θw偏离角/偏移角
        PHIW=PHIM*DSIN(2.0D0*PI*TAU_INIT/TAUC)!ϕw拍动角
        IF     (TAU_INIT>=0.0D0 .AND. TAU_INIT<=TAU_R1)THEN!ψw翻转角
            PSIW=PSI0-PSIM
        ELSE IF(TAU_INIT>TAU_R1 .AND. TAU_INIT<TAU_R1+DTAUR)THEN
            PSIW=PSI0-PSIM*DCOS(PI*(TAU_INIT-TAU_R1)/DTAUR)
        ELSE IF(TAU_INIT>=TAU_R1+DTAUR .AND. TAU_INIT<=TAU_R2)THEN
            PSIW=PSI0+PSIM
        ELSE IF(TAU_INIT>TAU_R2 .AND. TAU_INIT<TAU_R2+DTAUR)THEN
            PSIW=PSI0+PSIM*DCOS(PI*(TAU_INIT-TAU_R2)/DTAUR)
        ELSE IF(TAU_INIT>=TAU_R2+DTAUR .AND. TAU_INIT<TAUC)THEN
            PSIW=PSI0-PSIM
        END IF
    ELSE
        PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
        !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
        THETAW=0.0D0!θw偏离角/偏移角
        PHIW=PHIM*DSIN(2.0D0*PI*TAU/TAUC)!ϕw拍动角
        IF(TAU>=0.0D0 .AND. TAU<=TAU_R1)THEN!ψw翻转角
            PSIW=PSI0-PSIM
        ELSE IF(TAU>TAU_R1 .AND. TAU<TAU_R1+DTAUR)THEN
            PSIW=PSI0-PSIM*DCOS(PI*(TAU-TAU_R1)/DTAUR)
        ELSE IF(TAU>=TAU_R1+DTAUR .AND. TAU<=TAU_R2)THEN
            PSIW=PSI0+PSIM
        ELSE IF(TAU>TAU_R2 .AND. TAU<TAU_R2+DTAUR)THEN
            PSIW=PSI0+PSIM*DCOS(PI*(TAU-TAU_R2)/DTAUR)
        ELSE IF(TAU>=TAU_R2+DTAUR .AND. TAU<TAUC)THEN
            PSIW=PSI0-PSIM
        END IF
    END IF

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    IF(TIME+PHASE_DIFFERENCE/360.0D0*TAUC<-CRITERIA)THEN
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=0.0D0
        CEN_P(2)=0.0D0+PHIW*SPAN
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        CEN(1)=CEN(1)+CEN_TRANSLATION(1)
        CEN(2)=CEN(2)+CEN_TRANSLATION(2)
        VELO_ANGL=0.0D0
    ELSE
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=2.0D0*PI*PHIM/TAUC*DCOS(2.0D0*PI*TAU/TAUC)*SPAN
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        !CEN_P(1)=CEN_TRANSLATION(1)+0.0D0!-0.8D0
        !CEN_P(2)=CEN_TRANSLATION(2)+0.0D0+PHIW*SPAN
        CEN_P(1)=0.0D0
        CEN_P(2)=0.0D0+PHIW*SPAN
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        CEN(1)=CEN(1)+CEN_TRANSLATION(1)
        CEN(2)=CEN(2)+CEN_TRANSLATION(2)
        IF(TAU>=0.0D0 .AND. TAU<=TAU_R1)THEN!ψw翻转角
            VELO_ANGL=0.0D0
        ELSE IF(TAU>TAU_R1 .AND. TAU<TAU_R1+DTAUR)THEN
            VELO_ANGL= PI*PSIM/DTAUR*DSIN(PI*(TAU-TAU_R1)/DTAUR)
        ELSE IF(TAU>=TAU_R1+DTAUR .AND. TAU<=TAU_R2)THEN
            VELO_ANGL=0.0D0
        ELSE IF(TAU>TAU_R2 .AND. TAU<TAU_R2+DTAUR)THEN
            VELO_ANGL=-PI*PSIM/DTAUR*DSIN(PI*(TAU-TAU_R2)/DTAUR)
        ELSE IF(TAU>=TAU_R2+DTAUR .AND. TAU<TAUC)THEN
            VELO_ANGL=0.0D0
        END IF
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（周期性拍动）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PERIODIC_WANG(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !模拟1请确认符合模拟目标
    !--------------------本子函数根据不同扑翼规律需改变下方---------------------!
    !--------------------周期性拍动基本参数---------------------!
    TAUC=2.8D0*PI!7.725054831D0!周期时长
    !旋转时间/刻
    TAU_R=-TAUC/8.0D0
    !几何信息
    PHIM=1.0!拍动振幅
    PSIM=PI/4.0D0!翻转振幅
    PSI0=0.0D0*PI!初始翻转角

    SPAN=1.4D0!展长
    !--------------------本子函数根据不同扑翼规律需改变上方---------------------!

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据TAU确定各角度大小（周期性）---------------------!
    PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
    !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
    THETAW=0.0D0!θw偏离角/偏移角
    PHIW=PHIM*DSIN(2.0D0*PI*TAU/TAUC)!ϕw拍动角
    PSIW=PSI0-PSIM*DCOS(2.0D0*PI*(TAU-TAU_R)/TAUC)!ψw翻转角

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    !平动
    VELO_TRAN_R(1)=0.0D0
    VELO_TRAN_R(2)=2.0D0*PI*PHIM/TAUC*DCOS(2.0D0*PI*TAU/TAUC)*SPAN
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !转动
    CEN_P(1)=CEN_TRANSLATION(1)+0.0D0!-0.8D0
    CEN_P(2)=CEN_TRANSLATION(2)+0.0D0+PHIW*SPAN
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    VELO_ANGL= 2.0D0*PI*PSIM/TAUC*DSIN(2.0D0*PI*(TAU-TAU_R)/TAUC)

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（Heaving and Plunging）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PERIODIC_HEAVING_PLUNGING(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !--------------------周期性拍动基本参数---------------------!
    TAUC=2.0D0*PI/K!周期时长

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    PSI=0.0D0
    PSIW=0.0D0
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    IF(TIME+PHASE_DIFFERENCE/360.0D0*TAUC<-CRITERIA)THEN
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)+H*DSIN(K*TAU)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    ELSE
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=KH*DCOS(K*TAU)
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)+H*DSIN(K*TAU)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（OSCILLATING_X）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PERIODIC_OSCILLATING_X(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !--------------------周期性拍动基本参数---------------------!
    TAUC=2.0D0*PI/K!周期时长

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    PSI=0.0D0
    PSIW=0.0D0
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    IF(TIME+PHASE_DIFFERENCE/360.0D0*TAUC<-CRITERIA)THEN
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)-H*DSIN(K*TAU)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    ELSE
        !平动
        VELO_TRAN_R(1)=-KH*DCOS(K*TAU)
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)-H*DSIN(K*TAU)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（PURE_ROTATING）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PURE_ROTATING(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !--------------------周期性拍动基本参数---------------------!
    TAUC=PI**2.0D0!周期时长

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    PSI=0.0D0
    PSIW=0.0D0
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    IF(TIME+PHASE_DIFFERENCE/360.0D0*TAUC<-CRITERIA)THEN
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    ELSE
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=-2.0D0*PI**2.0D0/TAUC*DSIN(2.0D0*PI*TAU/TAUC)
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（PURE_ROTATING）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_PURE_ROTATING_STEADY(TIME)
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::TIME

    !--------------------周期性拍动基本参数---------------------!
    TAUC=PI!周期时长

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(TIME+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    PSI=0.0D0
    PSIW=0.0D0
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    IF(TIME+PHASE_DIFFERENCE/360.0D0*TAUC<-CRITERIA)THEN
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=0.0D0
    ELSE
        !平动
        VELO_TRAN_R(1)=0.0D0
        VELO_TRAN_R(2)=0.0D0
        VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
        !转动
        CEN_P(1)=CEN_TRANSLATION(1)
        CEN_P(2)=CEN_TRANSLATION(2)
        CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
        VELO_ANGL=-2.0D0
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（圆柱突然启动）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_IMPULSIVE_START
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE

    !--------------------突然启动基本参数---------------------!
    TAUC=0.01D0!0.0D0!加速时长
    DT=1.0D0/10000.0D0
    IF(NSTEP<=TAUC/DT)THEN
        T=DT*NSTEP
    ELSE
        T=1.0D0/DBLE(NCYCLE)*(NSTEP-TAUC/DT)+TAUC
        DT=1.0D0/DBLE(NCYCLE)
    END IF

    !--------------------根据TAU确定各角度大小（周期性）---------------------!
    PSI=0.0D0!此类算例中应为0°
    THETAW=0.0D0!θw偏离角/偏移角
    PHIW=0.0D0!ϕw拍动角
    PSIW=0.0D0!ψw翻转角

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    !平动
    IF(T<=TAUC)THEN
        VELO_TRAN_R(1)=-1.0D0/TAUC*T
    ELSE
        VELO_TRAN_R(1)=-1.0D0
    END IF
    VELO_TRAN_R(2)=0.0D0
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !转动
    IF(T<=TAUC)THEN
        CEN_P(1)=CEN_TRANSLATION(1)+0.5D0*VELO_TRAN_R(1)*TAUC
    ELSE
        CEN_P(1)=CEN_TRANSLATION(1)-0.5D0*1.0D0*TAUC-1.0D0*(T-TAUC)
    END IF
    CEN_P(2)=CEN_TRANSLATION(2)+0.0D0
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    VELO_ANGL=0.0D0

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（静止）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_STATIC
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE

    !--------------------本子函数根据不同扑翼规律需改变下方---------------------!
    !--------------------周期性拍动基本参数---------------------!
    TAUC=8.296510765D0!周期时长
    !旋转时间/刻
    TAU_R1=0.0625D0*TAUC
    TAU_R2=0.5625D0*TAUC!两个翻转开始的时刻
    DTAUR=0.25D0*TAUC!单次翻转时长
    !几何信息
    PHIM=32.5D0/180.0D0*PI!拍动振幅
    PSIM=50.0D0/180.0D0*PI!翻转振幅

    SPAN=4.556701031D0*0.8D0!展长
    !--------------------本子函数根据不同扑翼规律需改变上方---------------------!

    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间
    !!--------------------周期内时刻TAU---------------------!
    !TAU=DMOD(T+PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------给定各角度大小（恒定）---------------------!
    PSI=ABSX_UPSTROKE_ANGLE-90.0D0/180.0D0*PI!ψ拍动平面夹角，默认0°时上拍方向与绝对坐标系，
    !即计算坐标系Y轴正方向重合，即为ABSX_UPSTROKE_ANGLE-90°
    THETAW=0.0D0!θw偏离角/偏移角
    PHIW=PHIM*DSIN(2.0D0*PI*TAU_INIT/TAUC)!ϕw拍动角
    IF(TAU_INIT>=0.0D0 .AND. TAU_INIT<=TAU_R1)THEN!ψw翻转角
        PSIW=-PSIM
    ELSE IF(TAU_INIT>TAU_R1 .AND. TAU_INIT<TAU_R1+DTAUR)THEN
        PSIW=-PSIM*DCOS(PI*(TAU_INIT-TAU_R1)/DTAUR)
    ELSE IF(TAU_INIT>=TAU_R1+DTAUR .AND. TAU_INIT<=TAU_R2)THEN
        PSIW= PSIM
    ELSE IF(TAU_INIT>TAU_R2 .AND. TAU_INIT<TAU_R2+DTAUR)THEN
        PSIW= PSIM*DCOS(PI*(TAU_INIT-TAU_R2)/DTAUR)
    ELSE IF(TAU_INIT>=TAU_R2+DTAUR .AND. TAU_INIT<TAUC)THEN
        PSIW=-PSIM
    END IF

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    !平动
    VELO_TRAN_R(1)=0.0D0
    VELO_TRAN_R(2)=0.0D0
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !转动
    CEN_P(1)=CEN_TRANSLATION(1)+0.0D0!-0.8D0
    CEN_P(2)=CEN_TRANSLATION(2)+0.0D0+PHIW*SPAN
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    VELO_ANGL=0.0D0

    RETURN
    END SUBROUTINE

    !***************************************************转动中心位置；平转动速度；坐标转换所需（静止）******************************************************!
    SUBROUTINE POSE_VELO_QUADRIC_2D_STATIC_FOR_HEAVING_PLUNGING
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE

    !--------------------周期性拍动基本参数---------------------!
    TAUC=2.0D0*PI/K!周期时长
    TAU_INIT=PHASE_INITIATION/360.0D0*TAUC!由起始相位角推得起始时间

    !--------------------周期内时刻TAU---------------------!
    TAU=DMOD(PHASE_DIFFERENCE/360.0D0*TAUC+TAU_INIT,TAUC)

    !--------------------根据各角度确定坐标转换矩阵---------------------!
    PSI=0.0D0
    PSIW=0.0D0
    CALL CAL_TRANMAT(PSI,MATP)
    CALL CAL_TRANMAT(PSIW,MATW)

    TRANMAT=MATMUL( MATW,MATP )
    TRANMAT_INVERSE=TRANSPOSE(TRANMAT)

    !坐标转换矩阵系数
    T11=TRANMAT(1,1)
    T12=TRANMAT(1,2)
    T21=TRANMAT(2,1)
    T22=TRANMAT(2,2)

    !--------------------确定平动转动速度和中心（周期性）---------------------!
    !平动
    VELO_TRAN_R(1)=0.0D0
    VELO_TRAN_R(2)=0.0D0
    VELO_TRAN_A=MATMUL( TRANSPOSE(MATP),VELO_TRAN_R )
    !转动
    CEN_P(1)=CEN_TRANSLATION(1)
    CEN_P(2)=CEN_TRANSLATION(2)+H*DSIN(K*TAU)
    CEN=MATMUL( TRANSPOSE(MATP),CEN_P )
    VELO_ANGL=0.0D0

    RETURN
    END SUBROUTINE

    !***************************************************确定二次曲线表达式******************************************************!
    SUBROUTINE CAL_QUADRIC_COEFFICIENT
    USE QUADRIC_PARAMETER
    USE DECLARATION

    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE

    !给定相对坐标系下二次曲面的数学表达式系数
    IF(IB_SHAPE==1)THEN
        !圆
        RCOX2=1.0D0
        RCOY2=1.0D0
        RCOXY=0.0D0
        RCOX=0.0D0
        RCOY=0.0D0
        RCOM=-0.25D0!-0.0625D0!
    ELSE IF(IB_SHAPE==2)THEN
        !椭圆
        RCOX2=1.0D0
        RCOY2=100.0D0
        RCOXY=0.0D0
        RCOX=-2.0D0*CEN_DEVIATION(1)-0.0D0
        RCOY=-100.0D0*2.0D0*CEN_DEVIATION(2)-0.0D0
        RCOM=CEN_DEVIATION(1)**2.0D0+100D0*CEN_DEVIATION(2)**2.0D0-0.25D0!
        LAXIS=0.5D0
        SAXIS=0.05D0
        IF(IB_LOCOMOTION==11)THEN!WANG
            LAXIS=1.0D0/2.0D0
            SAXIS=0.125D0/2.0D0
            RCOX2=1.0D0/LAXIS**2.0D0
            RCOY2=1.0D0/SAXIS**2.0D0
            RCOXY=0.0D0
            RCOX=-2.0D0*CEN_DEVIATION(1)/LAXIS**2.0D0
            RCOY=-2.0D0*CEN_DEVIATION(2)/SAXIS**2.0D0
            RCOM=CEN_DEVIATION(1)**2.0D0/LAXIS**2.0D0&
                +CEN_DEVIATION(2)**2.0D0/SAXIS**2.0D0-1.0D0
        END IF
    END IF

    !根据坐标转换系数确定绝对坐标系下二次曲面的数学表达式系数
    PCOX2= RCOX2*T11**2.0D0 + RCOY2*T21**2.0D0 + RCOXY*T11*T21
    PCOY2= RCOX2*T12**2.0D0 + RCOY2*T22**2.0D0 + RCOXY*T12*T22
    PCOXY= 2.0D0*RCOX2*T11*T12 + 2.0D0*RCOY2*T21*T22 + RCOXY*(T11*T22+T21*T12)
    PCOX=  RCOX*T11 + RCOY*T21
    PCOY=  RCOX*T12 + RCOY*T22
    PCOM=  RCOM !- (RCOX*T11 + RCOY*T21)*CEN(1) - (RCOX*T12 + RCOY*T22)*CEN(2)

    !在绝对坐标系下对二次曲面的数学表达式进行平移
    COX2=PCOX2
    COY2=PCOY2
    COXY=PCOXY
    COX= PCOX-2.0D0*PCOX2*CEN(1)-PCOXY*CEN(2)
    COY= PCOY-2.0D0*PCOY2*CEN(2)-PCOXY*CEN(1)
    COM= PCOM+PCOX2*CEN(1)**2.0D0+PCOY2*CEN(2)**2.0D0+PCOXY*CEN(1)*CEN(2)-PCOX*CEN(1)-PCOY*CEN(2)

    !!给定上下表面XY大小范围，大者在上，XZ上下表面共用，Y一面一个
    !RXTRL=0.7D0
    !RXLED=-0.3D0
    !RYTOP=RCOMT/RCOY
    !RYBOT=RCOMB/RCOY

    !!上平面两端绝对坐标值
    !X1=RXLED*T11+RYTOP*T21+CEN(1)
    !X2=RXTRL*T11+RYTOP*T21+CEN(1)
    !
    !Y1=RXLED*T12+RYTOP*T22+CEN(2)
    !Y2=RXTRL*T12+RYTOP*T22+CEN(2)
    !
    !!上平面X,Y范围
    !XMAXT=DMAX1(X1,X2)
    !XMINT=DMIN1(X1,X2)
    !YMAXT=DMAX1(Y1,Y2)
    !YMINT=DMIN1(Y1,Y2)
    !
    !!下平面四角绝对坐标值
    !X1=RXLED*T11+RYBOT*T21+CEN(1)
    !X2=RXTRL*T11+RYBOT*T21+CEN(1)
    !
    !Y1=RXLED*T12+RYBOT*T22+CEN(2)
    !Y2=RXTRL*T12+RYBOT*T22+CEN(2)
    !
    !!下平面X,Y,Z范围
    !XMAXB=DMAX1(X1,X2)
    !XMINB=DMIN1(X1,X2)
    !YMAXB=DMAX1(Y1,Y2)
    !YMINB=DMIN1(Y1,Y2)

    RETURN
    END SUBROUTINE