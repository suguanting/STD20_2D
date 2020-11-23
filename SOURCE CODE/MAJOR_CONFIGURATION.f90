    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !**********************************************��Ҫ�������ã�ĸ�棩**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.02D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------�ڱ߽��˶����------!0-��ֹ�߽磬1-�����˶���2-��Ъ�Է��У�3-Բ��ͻȻ����,4-HEAVING&PLUNGING,5-CAVITY_OSCILLATING
    IF(IB_LOCOMOTION>=4)THEN
        KH=0.8D0
        K=8.0D0
        H=KH/K
    END IF
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=500.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/100.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=1600!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=2.0D0*PI/K/DBLE(NCYCLE)!1.0D0/DBLE(NCYCLE)!

    !------���������������------!
    NPROBE=200
    NCLCT=200
    NPLT=20
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0!23.0D0/180.0D0*PI!��Ŀǰ���㷨ֻ��Ϊ0
    ABSX_UPSTROKE_ANGLE=67.0D0/180.0D0*PI!0.0D0!90.0D0/180.0D0*PI!
    TRUX_FLIGHT_ANGLE=113.0D0/180.0D0*PI!90.0D0/180.0D0*PI!23.0D0/180.0D0*PI!
    ABSX_TRUX_ANGLE=67.0D0/180.0D0*PI!

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------�߽�����ϵ��------!
    !BC_A*U+BC_B
    !BC_A*DU
    !�ٶ�����ڣ�BC_A=0��BC_B=U_FREESTREAM���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0/����BC_A=0��BC_B=0
    !�ٶ��ҳ��ڣ�BC_A=1��BC_B=0           ���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0/����BC_A=0��BC_B=0
    !����������ΪԶ��������
    BCU_AL=1.0D0
    BCU_BL=0.0D0
    BCU_AR=1.0D0
    BCU_BR=0.0D0
    BCU_AB=1.0D0
    BCU_BB=0.0D0
    BCU_AT=1.0D0
    BCU_BT=0.0D0
    !������ڣ������޼��������ҳ��ڣ�
    !BCU_AL=0.0D0
    !BCU_BL=U_FREESTREAM
    !BCU_AR=1.0D0
    !BCU_BR=0.0D0
    !BCU_AB=1.0D0
    !BCU_BB=0.0D0
    !BCU_AT=1.0D0
    !BCU_BT=0.0D0
    !����������Ϊ�޻��ƹ̱ڱ�����
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
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0/����BC_A=0��BC_B=V_FREESTREAM/0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0/����BC_A=0��BC_B=0
    !����������ΪԶ��������
    BCV_AL=1.0D0
    BCV_BL=0.0D0
    BCV_AR=1.0D0
    BCV_BR=0.0D0
    BCV_AB=1.0D0
    BCV_BB=0.0D0
    BCV_AT=1.0D0
    BCV_BT=0.0D0
    !������ڣ������޼��������ҳ��ڣ�
    !BCV_AL=-1.0D0
    !BCV_BL=2.0D0*V_FREESTREAM
    !BCV_AR=1.0D0
    !BCV_BR=0.0D0
    !BCV_AB=0.0D0
    !BCV_BB=V_FREESTREAM
    !BCV_AT=0.0D0
    !BCV_BT=V_FREESTREAM
    !����������Ϊ�޻��ƹ̱ڱ�����
    !BCV_AL=-1.0D0
    !BCV_BL=0.0D0
    !BCV_AR=-1.0D0
    !BCV_BR=0.0D0
    !BCV_AB=0.0D0
    !BCV_BB=0.0D0
    !BCV_AT=0.0D0
    !BCV_BT=0.0D0

    !BC_A*PHI+BC_B
    !ѹ������/�޼���/�̱�:BC_A=1��BC_B=0��ѹ������:BC_A=0��BC_B=0
    !����������ΪԶ��������
    BCPHI_AL=0.0D0
    BCPHI_BL=0.0D0
    BCPHI_AR=0.0D0
    BCPHI_BR=0.0D0
    BCPHI_AB=0.0D0
    BCPHI_BB=0.0D0
    BCPHI_AT=0.0D0
    BCPHI_BT=0.0D0
    !������ڣ������޼��������ҳ��ڣ�
    !BCPHI_AL=1.0D0
    !BCPHI_BL=0.0D0
    !BCPHI_AR=0.0D0
    !BCPHI_BR=0.0D0
    !BCPHI_AB=1.0D0
    !BCPHI_BB=0.0D0
    !BCPHI_AT=1.0D0
    !BCPHI_BT=0.0D0
    !����������Ϊ�޻��ƹ̱ڱ�����
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

    !********************************************��Ҫ��������(3-Բ��ͻȻ����)****************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_IMPULSIVE_STARTED_CIRCULAR_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=2
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4�������5���һ���ٶȷֲ�
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=1
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=3000.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/120.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1 .OR. TASK_TYPE==5)THEN
        NCYCLE=600!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=1.0D0/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00550N000100.PLT"

    !------���������������------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=20
    NIB=200

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=2
    BCTYPE_R=2
    BCTYPE_B=2
    BCTYPE_T=2

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !**********************��Ҫ��������(5-CAVITY_OSCILLATING,7-PURE_ROTATING(UNSTEADY),8-PURE_ROTATING(STEADY))*******************!
    SUBROUTINE MAJOR_CONFIGURATION_CAVITY_CIRCULAR_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=0
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.01D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=2
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=4
    !------�ڱ߽��˶����------!5-CAVITY_OSCILLATING,7-PURE_ROTATING
    IF(IB_LOCOMOTION==5)THEN
        KH=1.0D0
        H=0.125D0
        K=KH/H
        !K=4.0D0/PI
    END IF
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=100.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/225.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/100.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=3600!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    IF(IB_LOCOMOTION==5) DT=2.0D0*PI/K/DBLE(NCYCLE)
    IF(IB_LOCOMOTION==7) DT=PI**2.0D0/DBLE(NCYCLE)
    IF(IB_LOCOMOTION==8) DT=PI/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00100N008000.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    IF(IB_LOCOMOTION==7)THEN
        NPLT=4
    ELSE IF(IB_LOCOMOTION==8)THEN
        NPLT=5
    END IF
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����������Ϊ�̱ڣ�
    BCTYPE_L=3
    BCTYPE_R=3
    BCTYPE_B=3
    BCTYPE_T=3

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !*********************************��Ҫ��������(-3-DRIVEN_CAVITY)*****************************************!
    SUBROUTINE MAJOR_CONFIGURATION_DRIVEN_CAVITY
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=0
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=2
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4�������5���һ���ٶȷֲ�
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=20.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/25.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=1000!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=1.0D0/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00400N012608.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=10
    IF(IB_LOCOMOTION==-4) NPLT=10
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=0
    BOUNDARY_EXISTENCE_2=0
    IF(IB_LOCOMOTION==-4) BOUNDARY_EXISTENCE_1=1

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !��������Ϊ�޻��ƹ̱ڱ���,��Ϊ�����߽�������
    BCTYPE_L=3
    BCTYPE_R=3
    BCTYPE_B=3
    BCTYPE_T=3

    !!!!!!!!!!!!!!!BCU_AT=-1.0D0
    !!!!!!!!!!!!!!!BCU_BT=2.0D0

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !*********************************��Ҫ��������(-2-ANALYTICAL_CASES)*****************************************!
    SUBROUTINE MAJOR_CONFIGURATION_ANALYTICAL_CASES
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=0
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=1
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=256.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    DX1 =1.0D0/ 3.0D0!���
    DX21=1.0D0/10.0D0!�����
    DX22=1.0D0/30.0D0!���в�
    DX23=1.0D0/60.0D0!���ڲ�
    DX3 =PI/16.0D0!�ڲ�

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=64!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=PI**2.0D0/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00100N008000.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=4!NCYCLE
    NIB=NCYCLE

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=0
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=0.75D0
    PROBE_Y1=0.0D0

    PROBE_X2=-0.75D0
    PROBE_Y2=0.0D0

    PROBE_X3=0.0D0
    PROBE_Y3=0.75D0

    PROBE_X4=0.0D0
    PROBE_Y4=-0.75D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !�����������޼��У�
    BCTYPE_L=4
    BCTYPE_R=4
    BCTYPE_B=4
    BCTYPE_T=4

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !***********************************************��Ҫ��������(HEAVING&PLUNGING)************************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_HEAVING_PLUNGING
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.02D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    IF(IB_LOCOMOTION==41)THEN
        CASE_TYPE=1
    ELSE IF(IB_LOCOMOTION==42)THEN
        CASE_TYPE=2
    END IF
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4�������5���һ���ٶȷֲ�
    TASK_TYPE=1
    !------�ڱ߽��˶����------!
    !KH=1.0D0!�����ٶ����ֵ
    !H=1.0D0!5.0D0/(2.0D0*PI)!���������ٻ����KH/K
    !K=KH/H!���������ٻ�Ƶ��
    H=1.0D0!5.0D0/(2.0D0*PI)!���������ٻ����KH/K
    K=1.0D0!���������ٻ�Ƶ��
    KH=K*H!�����ٶ����ֵ
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=1
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=500.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    END IF

    !------��������------!
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

    !------ȷ��ʱ�䲽------!
    DT=2.0D0*PI/K/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00500N029775.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=20
    NIB=NCYCLE

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=1
    BCTYPE_R=2
    BCTYPE_B=2
    BCTYPE_T=2

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !***********************************************��Ҫ��������************************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_X_OSCILLATING_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.02D0
        BR=1.02D0
        BB=1.02D0
        BT=1.02D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=2
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4�������5���һ���ٶȷֲ�
    TASK_TYPE=1
    !------�ڱ߽��˶����------!
    KH=1.0D0!0.8D0
    H=5.0D0/(2.0D0*PI)!KH/K
    K=KH/H!0.4D0*PI!8.0D0
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=1
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=100.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1 .OR. TASK_TYPE==5)THEN
        NCYCLE=600!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=2.0D0*PI/K/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00100N015300.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=60
    NIB=NCYCLE

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0!1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=0.4D0

    PROBE_X2=1.2D0
    PROBE_Y2=-0.4D0

    PROBE_X3=2.5D0
    PROBE_Y3=0.4D0

    PROBE_X4=2.5D0
    PROBE_Y4=-0.4D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !**********************************************��Ҫ��������(1-��ֹ������2,8-��Ъ�Է��У�����)**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_WING_FLAPPING
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.03D0
        BT=1.03D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=2
    !------���������ٲ�����׼����------!
    Re=150.0D0!1750.938894D0!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------������߶���------!
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

    LEIN=-2.0D0!-0.7D0!-1.2D0
    RIIN= 2.5D0! 0.7D0! 1.2D0
    BOIN=-3.0D0!-1.0D0!-2.5D0
    TOIN= 3.0D0! 1.0D0! 2.5D0

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=2000!1600!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=8.623671834D0/DBLE(NCYCLE)!7.725054831D0/DBLE(NCYCLE)!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe01580N040000.PLT"

    !------���������������------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=1
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    BOUNDARY_EXISTENCE_2=1!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0!ֻ��Ϊ0
    ABSX_UPSTROKE_ANGLE=60.0D0/180.0D0*PI!67.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    TRUX_FLIGHT_ANGLE=180.0D0/180.0D0*PI!113.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    ABSX_TRUX_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0!0.399262959D0!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !**********************************************��Ҫ��������(1-��ֹ������2,8-��Ъ�Է��У�����)**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_WING_FLAPPING_WANG
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.03D0
        BT=1.03D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=2
    !------���������ٲ�����׼����------!
    Re=75.0D0!1954.616861D0!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/80.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=800!100�ı���
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

    !------ȷ��ʱ�䲽------!
    DT=2.8D0*PI/DBLE(NCYCLE)!7.725054831D0/DBLE(NCYCLE)!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe01580N040000.PLT"

    !------���������������------!
    NPROBE=200
    NCLCT=NCYCLE
    NPLT=4
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    BOUNDARY_EXISTENCE_2=0!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0!ֻ��Ϊ0
    ABSX_UPSTROKE_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    TRUX_FLIGHT_ANGLE=180.0D0/180.0D0*PI!113.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��
    ABSX_TRUX_ANGLE=0.0D0/180.0D0*PI!67.0D0/180.0D0*PI!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=0.0D0!0.399262959D0!ģ��1��ȷ�Ϸ���ģ��Ŀ��

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=5
    BCTYPE_R=5
    BCTYPE_B=5
    BCTYPE_T=5

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !**********************************************��Ҫ�������ã���ֹԲ����**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_STATIC_CYLINDER
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.01D0
        BR=1.01D0
        BB=1.01D0
        BT=1.01D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=25.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=100!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=1.0D0/DBLE(NCYCLE)

    !------�����ļ���------!
    FILENAME_RESTART="2DXYRe00020N005607.PLT"

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=1.0D0/20.0D0
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=1
    BCTYPE_R=5
    BCTYPE_B=4
    BCTYPE_T=4

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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

    !**********************************************��Ҫ�������ã��ܵ���ֹԲ����**********************************************************!
    SUBROUTINE MAJOR_CONFIGURATION_STATIC_CYLINDER_IN_CHANNEL
    USE DECLARATION
    USE QUADRIC_PARAMETER
    USE CAL_QUADRIC_DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::CFLC=1.0D0,CFLV=1.0D0
    REAL(KIND=8)::DTC,DTV

    !------��������------!1-��������0-�ּ�����ĳЩֻʹ�þ���������������0������LEM1����ֵ��
    CONTINUOUS_MESH=1
    IF(CONTINUOUS_MESH==1)THEN
        BL=1.03D0
        BR=1.03D0
        BB=1.03D0
        BT=1.03D0
    END IF
    !------��������------!1����Ϊ����������������2����Ϊ��ʵ����������
    CASE_TYPE=1
    !------��������------!1�������㣬0����intersection��ֲ���Ƶ��2����˶����ɣ�3ת�����������4��������
    TASK_TYPE=1
    !------ճ������㷽��------!1������2������һʱ�䲽��������ɢ��ʽ
    VISCOUS_TERM_METHOD=2
    !------���߽���״------!1Բ��2��Բ
    IB_SHAPE=1
    !------���������ٲ�����׼����------!
    Re=20.0D0

    !------������߶���------!
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

    !------�����ܶ���------!
    IF(CONTINUOUS_MESH==0)THEN
        DX1 =1.0D0/ 3.0D0!���
        DX21=1.0D0/10.0D0!�����
        DX22=1.0D0/30.0D0!���в�
        DX23=1.0D0/60.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !DX1 =0.0D0!���
        !DX21=0.0D0!�����
        !DX22=0.0D0!���в�
        !DX23=0.0D0!���ڲ�
        DX3 =1.0D0/50.0D0!�ڲ�
    END IF

    !------��������------!
    IF(TASK_TYPE==1)THEN
        NCYCLE=100!1500!100�ı���
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

    !------ȷ��ʱ�䲽------!
    !DTC=CFLC*DX3/( MAXVAL(DABS(U)) + MAXVAL(DABS(V)) )
    !DTV=CFLV*Re*DX3/2.0D0
    !DT=DMIN1(DTV,DTC)
    DT=1.0D0/DBLE(NCYCLE)

    !------���������������------!
    NPROBE=NCYCLE
    NCLCT=NCYCLE
    NPLT=1
    NIB=100

    !------����߽�------!1-���ڣ�0-������
    BOUNDARY_EXISTENCE_1=1
    BOUNDARY_EXISTENCE_2=0

    !------��������ϵ��ת���------!
    FREESTREAM_TILT=0.0D0
    ABSX_UPSTROKE_ANGLE=0.0D0
    TRUX_FLIGHT_ANGLE=0.0D0
    ABSX_TRUX_ANGLE=0.0D0

    !------�����ٶȳ��������ٶȼ������ٶȳ�------!
    VELO_RATIO=1.0D0!

    U_FREESTREAM=1.0D0*VELO_RATIO*DCOS(FREESTREAM_TILT)
    V_FREESTREAM=1.0D0*VELO_RATIO*DSIN(FREESTREAM_TILT)

    !------����̽�루����ĸ���------!
    PROBE_X1=1.2D0
    PROBE_Y1=1.0D0

    PROBE_X2=1.2D0
    PROBE_Y2=-1.0D0

    PROBE_X3=2.5D0
    PROBE_Y3=1.0D0

    PROBE_X4=2.5D0
    PROBE_Y4=-1.0D0

    !------�߽�����ϵ��------!
    !1-����,2-����,3-�̱ڣ��޻��ƣ�,4-�޼��У����ɻ��ƣ�,5-���ڣ�����������
    !����Ϊ���ڣ�������Ϊ���ڣ�
    BCTYPE_L=1
    BCTYPE_R=2
    BCTYPE_B=3
    BCTYPE_T=3

    !BC_A*U+BC_B
    !BC_A*DU
    !BC_A*DU+BC_C*DUN רΪ������������
    !�ٶ�����ڣ�BC_A= 0��BC_B=U_FREESTREAM  ��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*V+BC_B
    !BC_A*DV
    !BC_A*DV+BC_C*DVN רΪ������������
    !�ٶ�����ڣ�BC_A=-1��BC_B=2*V_FREESTREAM��BC_C=0���޼��У����ɻ��ƣ�������BC_A= 1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��ҳ��ڣ�BC_A= 1��BC_B=0             ��BC_C=0���̱ڣ��޻��ƣ�    ������BC_A=-1��BC_B=0��BC_C=0/����BC_A=0��BC_B=0��BC_C=0
    !�ٶ��Ҷ������ڣ�BC_A=0��BC_B=0          ��BC_C=-U0*DT/DX
    !BC_A*PHI+BC_B
    !ѹ������/����/�޼���/�̱�:BC_A=1��BC_B=0
    IF(BCTYPE_L==1)THEN!�����
        BCU_AL=0.0D0
        BCU_BL=U_FREESTREAM
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=2.0D0*V_FREESTREAM
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==2)THEN!�����
        BCU_AL=1.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==3)THEN!��̱�
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=-1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==4)THEN!���޼���
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=0.0D0
        BCV_AL=1.0D0
        BCV_BL=0.0D0
        BCV_CL=0.0D0
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    ELSE IF(BCTYPE_L==5)THEN!����ڣ�������
        BCU_AL=0.0D0
        BCU_BL=0.0D0
        BCU_CL=-DT*U_FREESTREAM
        BCV_AL=0.0D0
        BCV_BL=0.0D0
        BCV_CL=-DT*U_FREESTREAM
        BCPHI_AL=1.0D0
        BCPHI_BL=0.0D0
    END IF
    IF(BCTYPE_R==1)THEN!�ҽ���
        BCU_AR=0.0D0
        BCU_BR=U_FREESTREAM
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=2.0D0*V_FREESTREAM
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==2)THEN!�ҳ���
        BCU_AR=1.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==3)THEN!�ҹ̱�
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=-1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==4)THEN!���޼���
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=0.0D0
        BCV_AR=1.0D0
        BCV_BR=0.0D0
        BCV_CR=0.0D0
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    ELSE IF(BCTYPE_R==5)THEN!�ҳ��ڣ�������
        BCU_AR=0.0D0
        BCU_BR=0.0D0
        BCU_CR=-DT*U_FREESTREAM
        BCV_AR=0.0D0
        BCV_BR=0.0D0
        BCV_CR=-DT*U_FREESTREAM
        BCPHI_AR=1.0D0
        BCPHI_BR=0.0D0
    END IF
    IF(BCTYPE_B==1)THEN!�½���
        BCU_AB=-1.0D0
        BCU_BB=2.0D0*U_FREESTREAM
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=V_FREESTREAM
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==2)THEN!�³���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=1.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==3)THEN!�¹̱�
        BCU_AB=-1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==4)THEN!���޼���
        BCU_AB=1.0D0
        BCU_BB=0.0D0
        BCU_CB=0.0D0
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=0.0D0
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    ELSE IF(BCTYPE_B==5)THEN!�³��ڣ�������
        BCU_AB=0.0D0
        BCU_BB=0.0D0
        BCU_CB=-DT*V_FREESTREAM
        BCV_AB=0.0D0
        BCV_BB=0.0D0
        BCV_CB=-DT*V_FREESTREAM
        BCPHI_AB=1.0D0
        BCPHI_BB=0.0D0
    END IF
    IF(BCTYPE_T==1)THEN!�Ͻ���
        BCU_AT=-1.0D0
        BCU_BT=2.0D0*U_FREESTREAM
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=V_FREESTREAM
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==2)THEN!�ϳ���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=1.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==3)THEN!�Ϲ̱�
        BCU_AT=-1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==4)THEN!���޼���
        BCU_AT=1.0D0
        BCU_BT=0.0D0
        BCU_CT=0.0D0
        BCV_AT=0.0D0
        BCV_BT=0.0D0
        BCV_CT=0.0D0
        BCPHI_AT=1.0D0
        BCPHI_BT=0.0D0
    ELSE IF(BCTYPE_T==5)THEN!�ϳ��ڣ�������
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