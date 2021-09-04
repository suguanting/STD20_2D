    !######################################################################!
    !#                                                                    #!
    !#                            ������                                  #!
    !#       �����������ںܱ����ض��ĳ�ʼ���÷�ʽ�������״��ƽ��         #!
    !#                                                                    #!
    !######################################################################!



    PROGRAM IBM_QUADRIC_MOVING
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !------ȷ�������˶�����------!
    !-4-��������������Բ��
    !-3-����������
    !-2-������1
    !-1-��ֹԲ��
    !0-��ֹ�߽磬1-�����˶���������2-��Ъ�Է��У�������11-WANGУ��
    !3-Բ��ͻȻ����
    !4-����ʱ����lock-in�о���41-����Բ��������42-HEAVING&PLUNGING,
    !5-CAVITY_OSCILLATING
    !6-X_OSCILLATING,7-PURE_ROTATING(UNSTEADY),8-PURE_ROTATING(STEADY)
    IB_LOCOMOTION=1
    !------�����˶�����������Ҫ��������------!
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
    ELSE IF(IB_LOCOMOTION==3)THEN
        CALL MAJOR_CONFIGURATION_IMPULSIVE_STARTED_CIRCULAR_CYLINDER
    ELSE IF(IB_LOCOMOTION==41 .OR. IB_LOCOMOTION==42)THEN
        CALL MAJOR_CONFIGURATION_HEAVING_PLUNGING
    ELSE IF(IB_LOCOMOTION==5 .OR. IB_LOCOMOTION==7 .OR. IB_LOCOMOTION==8)THEN
        CALL MAJOR_CONFIGURATION_CAVITY_CIRCULAR_CYLINDER
    ELSE IF(IB_LOCOMOTION==6)THEN
        CALL MAJOR_CONFIGURATION_X_OSCILLATING_CYLINDER
    END IF

    !------ȷ�������ģ------!
    CALL CAL_GRID_SIZE
    !------�������������ļ�------!
    CALL ALLOCATION_VARIABLE_FILE

    !------���̶ֹ������------!
    CALL MESHING_FIXED
    !------��ʼ��������ʱ�����------!
    CALL INITIATION

    !------ȷ��̽����Ϣ------!
    CALL CAL_PROBE

    !------�����������˵���ļ�------!
    CALL OUTPUT_CONFIGURATION_LOG

    !---------------------������������------------------------!
    !ת�����������
    IF(TASK_TYPE==3)THEN
        DO NSTEP=0,3100,50
            CALL READFROMFILE_2D_STAGGERED(FILENAME_RESTART)
            CALL OUTPUT_PLT_1_STAGGERED_RELATIVE
        END DO
        STOP
    END IF
    !��������
    IF(TASK_TYPE==4)THEN
        CALL CAL_ERROR
        STOP
    END IF
    !���һ���ٶȷֲ�
    IF(TASK_TYPE==5)THEN
        CALL OUTPUT_VELOCITY_PROFILE
        STOP
    END IF

    !---------------------��ʼ�������------------------------!
    DO NSTEP=NSTART,NMAX,1
        !------ȷ������ʱ��------!
        T=DT*DBLE(NSTEP)
        !------��⶯�߽�λ�ú��˶���Ϣ------!
        CALL CAL_QUADRIC_2D
        !------��ֻ������˶�����------!
        IF( TASK_TYPE==2 )CYCLE

        !------���ò����½���ʽ�߽�Ӱ������------!
        CALL IBM_INITIATION

        !------�߽�1------!
        IF( BOUNDARY_EXISTENCE_1==1 ) CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_1,QUADRIC_KINETIC_1)
        !------�߽�2------!
        IF( BOUNDARY_EXISTENCE_2==1 ) CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_2,QUADRIC_KINETIC_2)

        !------Ȩ��֮�ƣ�ȥ�������������------!
        CALL IBM_TYPE_FILTER
        !------���Ŀǰʱ���߽罻����һʱ����ٶ�------!
        !IF(IB_LOCOMOTION/=0)THEN
        !CALL IBM_CAL_VELON_FOR_BN1
        !END IF
        !!------IB����ԭʼ��Ϣת��Ϊ������Ϣ------!
        !CALL IBM_PRIMITIVE2DERIVATIVE
        !CALL IBM_PRIMITIVE2DERIVATIVE_2
        !CALL IBM_PRIMITIVE2DERIVATIVE_SECOND_ORDER

        !------����ཻ��������Ƶ������Ҫ��------!
        IF( TASK_TYPE==0 )THEN
            CALL OUTPUT_INTERSECTION
        END IF
        IF( TASK_TYPE==0 )CYCLE
        !------�������ʽ�߽�����ֲ�������Ҫ��------!
        !IF( NCYCLE>=100 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NIB) ) )==0 ) CALL OUTPUT_IB_STAGGERED

        !------��N-S���̽���ʱ���ƽ�����������RK-CN��ɢ��------!
        CALL TIMEADVANCE_INSE_FSM_RK3CN

        !!------��AB-CN��ɢ�Ķ������̽���ʱ���ƽ��õ�UHAT------!
        !CALL TIMEADVANCE1_ME_AB_CN
        !!CALL TIMEADVANCE1_ME_AB_CN_CONVECT
        !!CALL TIMEADVANCE1_ME_NOCONVECT_CN
        !!CALL TIMEADVANCE1_ME_NOCONVECT_EXPLICIT
        !!------���ѹ�����ɷ��̵õ�PHI------!
        !CALL TIMEADVANCE2_PPE
        !!------ʹ��PHI���µõ��ٶȳ���ѹ����������仯��------!
        !CALL TIMEADVANCE3_UPDATEUP_AB

        !!------�����������------!
        !!------�߽�1------!
        !IF( BOUNDARY_EXISTENCE_1==1 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        !    IF(IB_SHAPE==1)CALL CAL_CLCT_CIRCULAR
        !    IF(IB_SHAPE==2)CALL CAL_CLCT_ELLIPTIC(QUADRIC_GEOMETRICAL_1,QUADRIC_KINETIC_1,1)
        !END IF
        !!------�߽�2------!
        !IF( BOUNDARY_EXISTENCE_2==1 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        !    IF(IB_SHAPE==1)CALL CAL_CLCT_CIRCULAR
        !    IF(IB_SHAPE==2)CALL CAL_CLCT_ELLIPTIC(QUADRIC_GEOMETRICAL_2,QUADRIC_KINETIC_2,2)
        !END IF

        !------������������Ϣ------!
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
        !------��ȡ̽����Ϣ------!
        IF( NCYCLE>=100 .AND. MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPROBE) ) )==0 ) CALL OUTPUT_PROBE_STAGGERD_CONTINUOUS
        !------��Ļ���------!
        CALL POSTPROCESSING
        !------��Ļ���------!
        CALL OUTPUT_SCREEN
        !------��������------!
        CALL CONVERGENCE_JUDGE
        IF( CONVERGENCE==1 )EXIT

    END DO
    !---------------------����������------------------------!

    STOP
    END




















