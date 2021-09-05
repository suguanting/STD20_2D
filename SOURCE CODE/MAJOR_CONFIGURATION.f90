    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !**********************************************主要算例设置（母版）**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.02D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------内边界运动情况------!0-静止边界，1-周期运动，2-间歇性飞行，3-圆柱突然启动,4-HEAVING&PLUNGING,5-CAVITY_OSCILLATING
    IF(IB_LOCOMOTION>=4)THEN
        KH=0.8D0
        K=8.0D0
        H=KH/K
    END IF
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=500.0D0

    !------计算域尺度量------!
    LEFT=-7.5D0
    RIGH=15.0D0
    BOTT=-7.5D0
    TOPP= 7.5D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-1.2D0!-0.7D0!-1.2D0
    RIIN= 0.0D0! 0.7D0! 1.2D0
    BOIN=-0.7D0!-1.0D0!-2.5D0
    TOIN= 0.7D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/100.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=1600!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=500*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=2.0D0*PI/K/DBLE(NCYCLE)!1.0D0/DBLE(NCYCLE)!

    !------调用输出次数控制------!
    NPROBE=200
    NCLCT=200
    NPLT=20
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0!23.0D0/180.0D0*PI!以目前的算法只能为0
    ABSX_UPSTROKE_ANGLE=67.0D0/180.0D0*PI!0.0D0!90.0D0/180.0D0*PI!
    TRUX_FLIGHT_ANGLE=113.0D0/180.0D0*PI!90.0D0/180.0D0*PI!23.0D0/180.0D0*PI!
    ABSX_TRUX_ANGLE=67.0D0/180.0D0*PI!

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------边界条件系数------!
    !BC_A*U+BC_B
    !BC_A*DU
    !速度左进口：BC_A=0，BC_B=U_FREESTREAM；无剪切（自由滑移）：上下BC_A= 1，BC_B=0/左右BC_A=0，BC_B=0
    !速度右出口：BC_A=1，BC_B=0           ；固壁（无滑移）    ：上下BC_A=-1，BC_B=0/左右BC_A=0，BC_B=0
    !（上下左右为远场条件）
    BCU_AL=1.0D0
    BCU_BL=0.0D0
    BCU_AR=1.0D0
    BCU_BR=0.0D0
    BCU_AB=1.0D0
    BCU_BB=0.0D0
    BCU_AT=1.0D0
    BCU_BT=0.0D0
    !（左进口，上下无剪切流，右出口）
    !BCU_AL=0.0D0
    !BCU_BL=U_FREESTREAM
    !BCU_AR=1.0D0
    !BCU_BR=0.0D0
    !BCU_AB=1.0D0
    !BCU_BB=0.0D0
    !BCU_AT=1.0D0
    !BCU_BT=0.0D0
    !（上下左右为无滑移固壁边条）
    !BCU_AL=0.0D0
    !BCU_BL=0.0D0
    !BCU_AR=0.0D0
    !BCU_BR=0.0D0
    !BCU_AB=-1.0D0
    !BCU_BB=0.0D0
    !BCU_AT=-1.0D0
    !BCU_BT=0.0D0
    !BC_A*V+BC_B
    !BC_A*DV
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM；无剪切（自由滑移）：左右BC_A= 1，BC_B=0/上下BC_A=0，BC_B=V_FREESTREAM/0
    !速度右出口：BC_A= 1，BC_B=0             ；固壁（无滑移）    ：左右BC_A=-1，BC_B=0/上下BC_A=0，BC_B=0
    !（上下左右为远场条件）
    BCV_AL=1.0D0
    BCV_BL=0.0D0
    BCV_AR=1.0D0
    BCV_BR=0.0D0
    BCV_AB=1.0D0
    BCV_BB=0.0D0
    BCV_AT=1.0D0
    BCV_BT=0.0D0
    !（左进口，上下无剪切流，右出口）
    !BCV_AL=-1.0D0
    !BCV_BL=2.0D0*V_FREESTREAM
    !BCV_AR=1.0D0
    !BCV_BR=0.0D0
    !BCV_AB=0.0D0
    !BCV_BB=V_FREESTREAM
    !BCV_AT=0.0D0
    !BCV_BT=V_FREESTREAM
    !（上下左右为无滑移固壁边条）
    !BCV_AL=-1.0D0
    !BCV_BL=0.0D0
    !BCV_AR=-1.0D0
    !BCV_BR=0.0D0
    !BCV_AB=0.0D0
    !BCV_BB=0.0D0
    !BCV_AT=0.0D0
    !BCV_BT=0.0D0

    !BC_A*PHI+BC_B
    !压力进口/无剪切/固壁:BC_A=1，BC_B=0；压力出口:BC_A=0，BC_B=0
    !（上下左右为远场条件）
    BCPHI_AL=0.0D0
    BCPHI_BL=0.0D0
    BCPHI_AR=0.0D0
    BCPHI_BR=0.0D0
    BCPHI_AB=0.0D0
    BCPHI_BB=0.0D0
    BCPHI_AT=0.0D0
    BCPHI_BT=0.0D0
    !（左进口，上下无剪切流，右出口）
    !BCPHI_AL=1.0D0
    !BCPHI_BL=0.0D0
    !BCPHI_AR=0.0D0
    !BCPHI_BR=0.0D0
    !BCPHI_AB=1.0D0
    !BCPHI_BB=0.0D0
    !BCPHI_AT=1.0D0
    !BCPHI_BT=0.0D0
    !（上下左右为无滑移固壁边条）
    !BCPHI_AL=1.0D0
    !BCPHI_BL=0.0D0
    !BCPHI_AR=1.0D0
    !BCPHI_BR=0.0D0
    !BCPHI_AB=1.0D0
    !BCPHI_BB=0.0D0
    !BCPHI_AT=1.0D0
    !BCPHI_BT=0.0D0


    RETURN
    END SUBROUTINE

    !********************************************主要算例设置(3-圆柱突然启动)****************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_IMPULSIVE_STARTED_CIRCULAR_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=2
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差，5输出一列速度分布
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=1
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=3000.0D0

    !------计算域尺度量------!
    LEFT=-5.0D0!-7.5D0!-5.0D0!-10.0D0
    RIGH=10.0D0!15.0D0!10.0D0! 20.0D0
    BOTT=-5.0D0!-7.5D0!-5.0D0!-10.0D0
    TOPP= 5.0D0! 7.5D0! 5.0D0! 10.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-0.6D0!-1.2D0!-0.7D0!-1.2D0
    RIIN= 3.6D0! 0.0D0! 0.7D0! 1.2D0
    BOIN=-0.7D0!-0.7D0!-1.0D0!-2.5D0
    TOIN= 0.7D0! 0.7D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/120.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1 .OR. TASK_TYPE==5)THEN
        NCYCLE=600!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=2000
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=3*NCYCLE+100
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=1.0D0/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00550N000100.PLT"

    !------调用输出次数控制------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=20
    NIB=200

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=2
    BCTYPE_R=2
    BCTYPE_B=2
    BCTYPE_T=2

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !**********************主要算例设置(5-CAVITY_OSCILLATING,7-PURE_ROTATING(UNSTEADY),8-PURE_ROTATING(STEADY))*******************!
    SUBROUTINE MAJOR_CONFIGURATION_CAVITY_CIRCULAR_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=0
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.01D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=2
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=4
    !------内边界运动情况------!5-CAVITY_OSCILLATING,7-PURE_ROTATING
    IF(IB_LOCOMOTION==5)THEN
        KH=1.0D0
        H=0.125D0
        K=KH/H
        !K=4.0D0/PI
    END IF
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=100.0D0

    !------计算域尺度量------!
    LEFT=-1.0D0
    RIGH= 1.0D0
    BOTT=-1.0D0
    TOPP= 1.0D0

    LEM1=-1.0D0
    RIM1= 1.0D0
    BOM1=-1.0D0
    TOM1= 1.0D0

    LEM2=-1.0D0
    RIM2= 1.0D0
    BOM2=-1.0D0
    TOM2= 1.0D0

    LEM3=-1.0D0
    RIM3= 1.0D0
    BOM3=-1.0D0
    TOM3= 1.0D0

    LEIN=-1.0D0
    RIIN= 1.0D0
    BOIN=-1.0D0
    TOIN= 1.0D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/225.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/100.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=3600!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        IF(IB_LOCOMOTION==7)THEN
            NDURATION=1.5D0*NCYCLE
        ELSE IF(IB_LOCOMOTION==8)THEN
            NDURATION=5.0D0*NCYCLE
        END IF
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    IF(IB_LOCOMOTION==5) DT=2.0D0*PI/K/DBLE(NCYCLE)
    IF(IB_LOCOMOTION==7) DT=PI**2.0D0/DBLE(NCYCLE)
    IF(IB_LOCOMOTION==8) DT=PI/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00100N008000.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    IF(IB_LOCOMOTION==7)THEN
        NPLT=4
    ELSE IF(IB_LOCOMOTION==8)THEN
        NPLT=5
    END IF
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（上下左右为固壁）
    BCTYPE_L=3
    BCTYPE_R=3
    BCTYPE_B=3
    BCTYPE_T=3

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !*********************************主要算例设置(-3-DRIVEN_CAVITY)*****************************************!
    SUBROUTINE MAJOR_CONFIGURATION_DRIVEN_CAVITY
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=0
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=2
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差，5输出一列速度分布
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=20.0D0

    !------计算域尺度量------!
    LEFT= 0.0D0
    RIGH= 1.0D0
    BOTT= 0.0D0
    TOPP= 1.0D0

    LEM1= 0.0D0
    RIM1= 1.0D0
    BOM1= 0.0D0
    TOM1= 1.0D0

    LEM2= 0.0D0
    RIM2= 1.0D0
    BOM2= 0.0D0
    TOM2= 1.0D0

    LEM3= 0.0D0
    RIM3= 1.0D0
    BOM3= 0.0D0
    TOM3= 1.0D0

    LEIN= 0.0D0
    RIIN= 1.0D0
    BOIN= 0.0D0
    TOIN= 1.0D0

    IF(IB_LOCOMOTION==-4)THEN
        LEFT= -1.0D0
        BOTT= -1.0D0
        LEM1= -1.0D0
        BOM1= -1.0D0
        LEM2= -1.0D0
        BOM2= -1.0D0
        LEM3= -1.0D0
        BOM3= -1.0D0
        LEIN= -1.0D0
        BOIN= -1.0D0
    END IF

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/25.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=1000!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=800*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=1.0D0/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00400N012608.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=10
    IF(IB_LOCOMOTION==-4) NPLT=10
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=0
    BOUNDARY_EXISTENCE_2=0
    IF(IB_LOCOMOTION==-4) BOUNDARY_EXISTENCE_1=1

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（下左右为无滑移固壁边条,上为驱动边界条件）
    BCTYPE_L=3
    BCTYPE_R=3
    BCTYPE_B=3
    BCTYPE_T=3

    !!!!!!!!!!!!!!!BCU_AT=-1.0D0
    !!!!!!!!!!!!!!!BCU_BT=2.0D0

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !*********************************主要算例设置(-2-ANALYTICAL_CASES)*****************************************!
    SUBROUTINE MAJOR_CONFIGURATION_ANALYTICAL_CASES
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=0
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=1
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=256.0D0

    !------计算域尺度量------!
    LEFT=-PI
    RIGH= PI
    BOTT=-PI
    TOPP= PI

    LEM1=-PI
    RIM1= PI
    BOM1=-PI
    TOM1= PI

    LEM2=-PI
    RIM2= PI
    BOM2=-PI
    TOM2= PI

    LEM3=-PI
    RIM3= PI
    BOM3=-PI
    TOM3= PI

    LEIN=-PI
    RIIN= PI
    BOIN=-PI
    TOIN= PI

    !------网格密度量------!
    DX1 =1.0D0/ 3.0D0!外层
    DX21=1.0D0/10.0D0!中外层
    DX22=1.0D0/30.0D0!中中层
    DX23=1.0D0/60.0D0!中内层
    DX3 =PI/16.0D0!内层

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=64!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=2*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=PI**2.0D0/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00100N008000.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=4!NCYCLE
    NIB=NCYCLE

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=0
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（上下左右无剪切）
    BCTYPE_L=4
    BCTYPE_R=4
    BCTYPE_B=4
    BCTYPE_T=4

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !***********************************************主要算例设置(HEAVING&PLUNGING)************************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_HEAVING_PLUNGING
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.02D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    IF(IB_LOCOMOTION==41)THEN
        CASE_TYPE=1
    ELSE IF(IB_LOCOMOTION==42)THEN
        CASE_TYPE=2
    END IF
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差，5输出一列速度分布
    TASK_TYPE=1
    !------内边界运动情况------!
    !KH=1.0D0!表征速度最大值
    !H=1.0D0!5.0D0/(2.0D0*PI)!表征无量纲化振幅KH/K
    !K=KH/H!表征无量纲化频率
    H=1.0D0!5.0D0/(2.0D0*PI)!表征无量纲化振幅KH/K
    K=1.0D0!表征无量纲化频率
    KH=K*H!表征速度最大值
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=1
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=500.0D0

    !------计算域尺度量------!
    LEFT=-7.5D0!-10.0D0!-15.0D0!-10.0D0!-20.0D0!
    RIGH=15.0D0! 10.0D0! 15.0D0! 10.0D0! 20.0D0!
    BOTT=-7.5D0!-10.0D0!-15.0D0!-10.0D0!-20.0D0!
    TOPP= 7.5D0! 10.0D0! 15.0D0! 10.0D0! 20.0D0!

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-0.8D0!-1.2D0!-0.7D0!-1.2D0
    RIIN= 1.8D0! 0.0D0! 0.7D0! 1.2D0
    BOIN=-1.8D0!-0.7D0!-1.0D0!-2.5D0
    TOIN= 1.8D0! 0.7D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1 .OR. TASK_TYPE==5)THEN
        NCYCLE=1000
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=500*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=2.0D0*PI/K/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00500N029775.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=20
    NIB=NCYCLE

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=1
    BCTYPE_R=2
    BCTYPE_B=2
    BCTYPE_T=2

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !***********************************************主要算例设置************************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_X_OSCILLATING_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.02D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=2
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差，5输出一列速度分布
    TASK_TYPE=1
    !------内边界运动情况------!
    KH=1.0D0!0.8D0
    H=5.0D0/(2.0D0*PI)!KH/K
    K=KH/H!0.4D0*PI!8.0D0
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=1
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=100.0D0

    !------计算域尺度量------!
    LEFT=-25.0D0!-15.0D0!-10.0D0!-7.5D0!-20.0D0!
    RIGH= 25.0D0! 15.0D0! 10.0D0!15.0D0! 20.0D0!
    BOTT=-15.0D0!-15.0D0!-10.0D0!-7.5D0!-20.0D0!
    TOPP= 15.0D0! 15.0D0! 10.0D0! 7.5D0! 20.0D0!

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-3.0D0!-1.2D0!-0.7D0!-1.2D0
    RIIN= 3.0D0! 0.0D0! 0.7D0! 1.2D0
    BOIN=-2.0D0!-0.7D0!-1.0D0!-2.5D0
    TOIN= 2.0D0! 0.7D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1 .OR. TASK_TYPE==5)THEN
        NCYCLE=600!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=10*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=2.0D0*PI/K/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00100N015300.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=60
    NIB=NCYCLE

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0!1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !**********************************************主要算例设置(1-静止初场，2,8-间歇性飞行（扑翼）)**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_WING_FLAPPING
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.02D0
        BB=1.03D0
        BT=1.01D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=2
    !------各种无量纲参数和准则数------!
    Re=1167.292596!1750.938894D0!150.0D0!模拟1请确认符合模拟目标

    !------计算域尺度量------!
    LEFT=-15.0D0
    RIGH= 15.0D0
    BOTT=-15.0D0
    TOPP= 15.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-2.5D0!-0.7D0!-1.2D0
    RIIN= 2.5D0! 0.7D0! 1.2D0
    BOIN=-1.0D0!-1.0D0!-2.5D0
    TOIN= 1.5D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=2400!1600!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=400
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=0.5*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=5.749114556D0/DBLE(NCYCLE)!8.623671834D0/DBLE(NCYCLE)!7.725054831D0/DBLE(NCYCLE)!模拟1请确认符合模拟目标

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe01580N040000.PLT"

    !------调用输出次数控制------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=20
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1!模拟1请确认符合模拟目标
    BOUNDARY_EXISTENCE_2=0!模拟1请确认符合模拟目标

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=90.0D0/180.0D0*PI!只能为0
    ABSX_UPSTROKE_ANGLE=0.0D0/180.0D0*PI!60.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    TRUX_FLIGHT_ANGLE=180.0D0/180.0D0*PI!113.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    ABSX_TRUX_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!模拟1请确认符合模拟目标

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0!0.399262959D0!模拟1请确认符合模拟目标

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=-1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !**********************************************主要算例设置(12-极限负载飞行)**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_WING_FLAPPING_MAXIMUM
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.03D0
        BB=1.03D0
        BT=1.02D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=2
    !------各种无量纲参数和准则数------!
    Re=1348.087485D0!模拟1请确认符合模拟目标

    !------计算域尺度量------!
    LEFT=-15.0D0
    RIGH= 15.0D0
    BOTT=-15.0D0
    TOPP= 15.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-2.5D0!
    RIIN= 2.5D0!
    BOIN=-2.0D0!
    TOIN= 1.5D0!

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=2400!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=400
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=20*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=5.412684889D0/DBLE(NCYCLE)!模拟1请确认符合模拟目标

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe01580N040000.PLT"

    !------调用输出次数控制------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=4
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1!模拟1请确认符合模拟目标
    BOUNDARY_EXISTENCE_2=1!模拟1请确认符合模拟目标

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=-19.03443435D0/180.0D0*PI!
    ABSX_UPSTROKE_ANGLE=0.0D0/180.0D0*PI!60.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    TRUX_FLIGHT_ANGLE=180.0D0/180.0D0*PI!113.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    ABSX_TRUX_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!模拟1请确认符合模拟目标

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=1.062158495D0!模拟1请确认符合模拟目标

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=-1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=1
    BCTYPE_R=5
    BCTYPE_B=1
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF

    RETURN
    END SUBROUTINE
    
    !**********************************************主要算例设置(1-静止初场，2,8-间歇性飞行（扑翼）)**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_WING_FLAPPING_WANG
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.03D0
        BT=1.03D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=2
    !------各种无量纲参数和准则数------!
    Re=75.0D0!1954.616861D0!模拟1请确认符合模拟目标

    !------计算域尺度量------!
    LEFT=-15.0D0
    RIGH= 15.0D0
    BOTT=-15.0D0
    TOPP= 15.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-2.5D0!-0.7D0!-1.2D0
    RIIN= 2.5D0! 0.7D0! 1.2D0
    BOIN=-2.0D0!-1.0D0!-2.5D0
    TOIN= 1.0D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/80.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/50.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=800!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=20*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    DT=2.8D0*PI/DBLE(NCYCLE)!7.725054831D0/DBLE(NCYCLE)!模拟1请确认符合模拟目标

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe01580N040000.PLT"

    !------调用输出次数控制------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=4
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1!模拟1请确认符合模拟目标
    BOUNDARY_EXISTENCE_2=0!模拟1请确认符合模拟目标

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0!只能为0
    ABSX_UPSTROKE_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    TRUX_FLIGHT_ANGLE=180.0D0/180.0D0*PI!113.0D0/180.0D0*PI!模拟1请确认符合模拟目标
    ABSX_TRUX_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!模拟1请确认符合模拟目标

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=0.0D0!0.399262959D0!模拟1请确认符合模拟目标

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF

    RETURN
    END SUBROUTINE

    !**********************************************主要算例设置（静止圆柱）**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_STATIC_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.01D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=25.0D0

    !------计算域尺度量------!
    LEFT=-10.0D0
    RIGH= 10.0D0
    BOTT=-10.0D0
    TOPP= 10.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-1.0D0!-0.7D0!-1.2D0
    RIIN= 1.0D0! 0.7D0! 1.2D0
    BOIN=-1.0D0!-1.0D0!-2.5D0
    TOIN= 1.0D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/50.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/50.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=100!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=1000*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=1.0D0/DBLE(NCYCLE)

    !------续算文件名------!
    FILENAME_RESTART="2DXYRe00020N005607.PLT"

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=1.0D0/20.0D0
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=1
    BCTYPE_R=5
    BCTYPE_B=4
    BCTYPE_T=4

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE

    !**********************************************主要算例设置（管道静止圆柱）**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_STATIC_CYLINDER_IN_CHANNEL
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------网格种类------!1-连续网格；0-分级网格（某些只使用均匀网格算例，置0并更改LEM1等数值）
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.03D0
        BT=1.03D0
    END IF
    !------算例种类------!1初场为均匀流场的算例，2初场为真实流场的算例
    CASE_TYPE=1
    !------任务种类------!1正常计算，0生成intersection点分布视频，2输出运动规律，3转换相对流场，4求解结果误差
    TASK_TYPE=1
    !------粘性项计算方法------!1正常，2根据上一时间步决定其离散方式
    VISCOUS_TERM_METHOD=2
    !------动边界形状------!1圆，2椭圆
    IB_SHAPE=1
    !------各种无量纲参数和准则数------!
    Re=20.0D0

    !------计算域尺度量------!
    LEFT=-15.0D0
    RIGH= 15.0D0
    BOTT=-10.0D0
    TOPP= 10.0D0

    !LEM1=-4.0D0
    !RIM1= 4.0D0
    !BOM1=-5.0D0
    !TOM1= 5.0D0
    !
    !LEM2=-2.0D0
    !RIM2= 2.0D0
    !BOM2=-4.0D0
    !TOM2= 4.0D0
    !
    !LEM3=-1.6D0
    !RIM3= 1.6D0
    !BOM3=-3.0D0
    !TOM3= 3.0D0

    LEIN=-1.0D0!-0.7D0!-1.2D0
    RIIN= 1.0D0! 0.7D0! 1.2D0
    BOIN=-1.0D0!-1.0D0!-2.5D0
    TOIN= 1.0D0! 1.0D0! 2.5D0

    !------网格密度量------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!外层
        DX21=1.0D0/10.0D0!中外层
        DX22=1.0D0/30.0D0!中中层
        DX23=1.0D0/60.0D0!中内层
        DX3 =1.0D0/50.0D0!内层
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!外层
        !DX21=0.0D0!中外层
        !DX22=0.0D0!中中层
        !DX23=0.0D0!中内层
        DX3 =1.0D0/50.0D0!内层
    END IF

    !------迭代控制------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=100!1500!100的倍数
    ELSE IF(TASK_TYPE==0)THEN
        NCYCLE=250
    ELSE IF(TASK_TYPE==2)THEN
        NCYCLE=500
    END IF
    IF(TASK_TYPE==1)THEN
        NDURATION=500*NCYCLE
    ELSE IF(TASK_TYPE==0)THEN
        NDURATION=NCYCLE
    ELSE IF(TASK_TYPE==2)THEN
        NDURATION=10*NCYCLE
    END IF

    !------确定时间步------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=1.0D0/DBLE(NCYCLE)

    !------调用输出次数控制------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=1
    NIB=100

    !------激活边界------!1-存在，0-不存在
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------扑翼坐标系旋转相关------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------背景速度场无量纲速度及背景速度场------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------布置探针（最多四个）------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------边界条件系数------!
    !1-进口,2-出口,3-固壁（无滑移）,4-无剪切（自由滑移）,5-出口（对流边条）
    !（左为进口，上下右为出口）
    BCTYPE_L=1
    BCTYPE_R=2
    BCTYPE_B=3
    BCTYPE_T=3

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN 专为对流条件而设
    !速度左进口：BC_A= 0，BC_B=U_FREESTREAM  ，BC_C=0；无剪切（自由滑移）：上下BC_A= 1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：上下BC_A=-1，BC_B=0，BC_C=0/左右BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN 专为对流条件而设
    !速度左进口：BC_A=-1，BC_B=2*V_FREESTREAM，BC_C=0；无剪切（自由滑移）：左右BC_A= 1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右出口：BC_A= 1，BC_B=0             ，BC_C=0；固壁（无滑移）    ：左右BC_A=-1，BC_B=0，BC_C=0/上下BC_A=0，BC_B=0，BC_C=0
    !速度右对流出口：BC_A=0，BC_B=0          ，BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !压力进口/出口/无剪切/固壁:BC_A=1，BC_B=0
    IF(BCTYPE_L==1)THEN!左进口
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!左出口
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!左固壁
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!左无剪切
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!左出口（对流）
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!右进口
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!右出口
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!右固壁
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!右无剪切
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!右出口（对流）
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!下进口
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!下出口
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!下固壁
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!下无剪切
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!下出口（对流）
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!上进口
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!上出口
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!上固壁
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!上无剪切
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!上出口（对流）
        BCU_AT=0.0D0
        BCU_BT=0.0D0
        BCU_CT=-DT*V_FREESTREAM
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=-DT*V_FREESTREAM
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    END IF


    RETURN
    END SUBROUTINE