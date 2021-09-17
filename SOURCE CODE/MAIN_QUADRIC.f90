    !######################################################################!
    !#                                                                    #!
    !#                            主程序                                  #!
    !#       本程序适用于很薄，特定的初始放置方式，翅膀形状的平板         #!
    !#                                                                    #!
    !######################################################################!



    PROGRAM IBM_QUADRIC_MOVING
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !------确定物体运动规律------!
    !-4-顶盖驱动流包含圆柱
    !-3-顶盖驱动流
    !-2-解析解1
    !-1-静止圆柱
    !0-静止边界，1-周期运动（扑翼），2-间歇性飞行（扑翼），11-WANG校核
    !12-极限负载飞行，13-极限负载飞行（水平来流），14-参考前飞-SUN，15-参考前飞-自设
    !3-圆柱突然启动
    !4-（暂时用于lock-in研究）41-静置圆柱初场，42-HEAVING&PLUNGING,
    !5-CAVITY_OSCILLATING
    !6-X_OSCILLATING,7-PURE_ROTATING(UNSTEADY),8-PURE_ROTATING(STEADY)
    IB_LOCOMOTION=13
    !------根据运动规律类型主要算例设置------!
    IF(IB_LOCOMOTION==0)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING
    ELSE IF(IB_LOCOMOTION==-3 .OR. IB_LOCOMOTION==-4)THEN
        CALL MAJOR_CONFIGURATION_DRIVEN_CAVITY
    ELSE IF(IB_LOCOMOTION==-2)THEN
        CALL MAJOR_CONFIGURATION_ANALYTICAL_CASES
    ELSE IF(IB_LOCOMOTION==-1)THEN
        CALL MAJOR_CONFIGURATION_STATIC_CYLINDER
        !CALL MAJOR_CONFIGURATION_STATIC_CYLINDER_IN_CHANNEL
    ELSE IF(IB_LOCOMOTION==1 .OR. IB_LOCOMOTION==2)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING
    ELSE IF(IB_LOCOMOTION==11)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING_WANG
    ELSE IF(IB_LOCOMOTION==12)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING_MAXIMUM
    ELSE IF(IB_LOCOMOTION==13)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING_MAXIMUM_HORI
    ELSE IF(IB_LOCOMOTION==14)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING_FORWARD_SUN
    ELSE IF(IB_LOCOMOTION==15)THEN
        CALL MAJOR_CONFIGURATION_WING_FLAPPING_MAXIMUM_HORI
    ELSE IF(IB_LOCOMOTION==3)THEN
        CALL MAJOR_CONFIGURATION_IMPULSIVE_STARTED_CIRCULAR_CYLINDER
    ELSE IF(IB_LOCOMOTION==41 .OR. IB_LOCOMOTION==42)THEN
        CALL MAJOR_CONFIGURATION_HEAVING_PLUNGING
    ELSE IF(IB_LOCOMOTION==5 .OR. IB_LOCOMOTION==7 .OR. IB_LOCOMOTION==8)THEN
        CALL MAJOR_CONFIGURATION_CAVITY_CIRCULAR_CYLINDER
    ELSE IF(IB_LOCOMOTION==6)THEN
        CALL MAJOR_CONFIGURATION_X_OSCILLATING_CYLINDER
    END IF

    !------确定网格规模------!
    CALL CAL_GRID_SIZE
    !------分配变量数组和文件------!
    CALL ALLOCATION_VARIABLE_FILE

    !------划分固定网格点------!
    CALL MESHING_FIXED
    !------初始化流场和时间层跨度------!
    CALL INITIATION

    !------确定探针信息------!
    CALL CAL_PROBE

    !------输出算例设置说明文件------!
    CALL OUTPUT_CONFIGURATION_LOG

    !---------------------处理流场数据------------------------!
    !转换相对流场等
    IF(TASK_TYPE==3)THEN
        DO NSTEP=0,3100,50
            CALL READFROMFILE_2D_STAGGERED(FILENAME_RESTART)
            CALL OUTPUT_PLT_1_STAGGERED_RELATIVE
        END DO
        STOP
    END IF
    !求解结果误差
    IF(TASK_TYPE==4)THEN
        CALL CAL_ERROR
        STOP
    END IF
    !输出一列速度分布
    IF(TASK_TYPE==5)THEN
        CALL OUTPUT_VELOCITY_PROFILE
        STOP
    END IF

    !---------------------开始迭代求解------------------------!
    DO NSTEP=NSTART,NMAX,1
        !------确定物理时间------!
        T=DT*DBLE(NSTEP)
        !------求解动边界位置和运动信息------!
        CALL CAL_QUADRIC_2D
        !------若只是输出运动规律------!
        IF( TASK_TYPE==2 )CYCLE

        !------重置并更新浸入式边界影响区域------!
        CALL IBM_INITIATION

        !------边界1------!
        IF( BOUNDARY_EXISTENCE_1==1 ) CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_1,QUADRIC_KINETIC_1)
        !------边界2------!
        IF( BOUNDARY_EXISTENCE_2==1 ) CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_2,QUADRIC_KINETIC_2)

        !------权宜之计：去除相邻两向外点------!
        CALL IBM_TYPE_FILTER
        !------求解目前时间层边界交点上一时间层速度------!
        !IF(IB_LOCOMOTION/=0)THEN
        !CALL IBM_CAL_VELON_FOR_BN1
        !END IF
        !!------IB方法原始信息转换为所需信息------!
        !CALL IBM_PRIMITIVE2DERIVATIVE
        !CALL IBM_PRIMITIVE2DERIVATIVE_2
        !CALL IBM_PRIMITIVE2DERIVATIVE_SECOND_ORDER

        !------输出相交点生成视频（如需要）------!
        IF( TASK_TYPE==0 )THEN
            CALL OUTPUT_INTERSECTION
        END IF
        IF( TASK_TYPE==0 )CYCLE
        !------输出浸入式边界参数分布（如需要）------!
        !IF( NCYCLE>=100 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NIB) ) )==0 ) CALL OUTPUT_IB_STAGGERED

        !------对N-S方程进行时间推进（动量方程RK-CN离散）------!
        CALL TIMEADVANCE_INSE_FSM_RK3CN

        !!------对AB-CN离散的动量方程进行时间推进得到UHAT------!
        !CALL TIMEADVANCE1_ME_AB_CN
        !!CALL TIMEADVANCE1_ME_AB_CN_CONVECT
        !!CALL TIMEADVANCE1_ME_NOCONVECT_CN
        !!CALL TIMEADVANCE1_ME_NOCONVECT_EXPLICIT
        !!------求解压力泊松方程得到PHI------!
        !CALL TIMEADVANCE2_PPE
        !!------使用PHI更新得到速度场和压力场并计算变化量------!
        !CALL TIMEADVANCE3_UPDATEUP_AB

        !!------求解气动参数------!
        !!------边界1------!
        !IF( BOUNDARY_EXISTENCE_1==1 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        !    IF(IB_SHAPE==1)CALL CAL_CLCT_CIRCULAR
        !    IF(IB_SHAPE==2)CALL CAL_CLCT_ELLIPTIC(QUADRIC_GEOMETRICAL_1,QUADRIC_KINETIC_1,1)
        !END IF
        !!------边界2------!
        !IF( BOUNDARY_EXISTENCE_2==1 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        !    IF(IB_SHAPE==1)CALL CAL_CLCT_CIRCULAR
        !    IF(IB_SHAPE==2)CALL CAL_CLCT_ELLIPTIC(QUADRIC_GEOMETRICAL_2,QUADRIC_KINETIC_2,2)
        !END IF

        !------求解输出流场信息------!
        IF(MOD( NSTEP,IDNINT( DBLE(NCYCLE)/NPLT ) )==0 .OR. NSTEP==NSTART )THEN! .OR. NSTEP>=2000
            CALL OUTPUT_PLT_1_STAGGERED
            !CALL OUTPUT_FULL_STAGGERED
            !CALL TEST_OUTPUT
            IF(IB_LOCOMOTION==-2)THEN
                !CALL ANALYTICAL_CASE_ERROR_ANALYSIS
                CALL ANALYTICAL_CASE_ERROR_OUTPUT
            END IF
            IF(IB_LOCOMOTION==3)THEN
                CALL OUTPUT_PLT_1_STAGGERED_RELATIVE
            END IF
        END IF
        !------读取探针信息------!
        IF( NCYCLE>=100 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPROBE) ) )==0 ) CALL OUTPUT_PROBE_STAGGERD_CONTINUOUS
        !------屏幕输出------!
        CALL POSTPROCESSING
        !------屏幕输出------!
        CALL OUTPUT_SCREEN
        !------跳出设置------!
        CALL CONVERGENCE_JUDGE
        IF( CONVERGENCE==1 )EXIT

    END DO
    !---------------------迭代求解结束------------------------!

    STOP
    END




















