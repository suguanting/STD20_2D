    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************�����Բ����������Ч��******************************************************!
    SUBROUTINE CAL_CLCT_ELLIPTIC(ELLIPSE_GEOMETRICAL,ELLIPSE_KINETIC,BOUNDARY_ID)
    USE QUADRIC_PARAMETER
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::ELLIPSE_GEOMETRICAL(36),ELLIPSE_KINETIC(9)
    INTEGER::BOUNDARY_ID
    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    CHARACTER(LEN=2)CHAR_ID
    INTEGER::REYNOLDS

    !��������
    REAL(KIND=8)::CEN_ELP(2)
    REAL(KIND=8)::CEN_DEVIATION(2)
    REAL(KIND=8)::LAXIS,SAXIS
    !����ת������
    REAL(KIND=8)::MAT_ABS2REL(2,2),MAT_REL2ABS(2,2),MAT_ABS2TRU(2,2),MAT_ABS2FLP(2,2)

    !��������ϵ�¶����������ѧ���ʽϵ��
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !�����ܶ�
    REAL(KIND=8)::RDA,RDN!Բ�Ľǣ�����
    !����ֵ
    REAL(KIND=8),ALLOCATABLE::XR(:,:),YR(:,:)
    REAL(KIND=8),ALLOCATABLE::XA(:,:),YA(:,:)
    REAL(KIND=8),ALLOCATABLE::NXR(:),NYR(:)
    REAL(KIND=8),ALLOCATABLE::NXA(:),NYA(:)
    REAL(KIND=8),ALLOCATABLE::ARC(:)
    REAL(KIND=8)::XATEMP,YATEMP
    !������
    INTEGER::RIM
    INTEGER::RJM
    INTEGER::RI,RJ,RJIM
    !---------��ֵ�漰�ľ����������---------!
    !�ǵ���Ϣ
    INTEGER::IAU,JAU!���½ǽű�ֵ
    INTEGER::IAV,JAV!���½ǽű�ֵ
    INTEGER::IAP,JAP!���½ǽű�ֵ

    !��ʵ��������ϵ���������ϵ���³����һ����ٶ�
    REAL(KIND=8)::UT,VT
    !�ٶ�ѹ����
    REAL(KIND=8),ALLOCATABLE::UR(:,:),VR(:,:),PR(:,:)
    REAL(KIND=8),ALLOCATABLE::UA(:,:),VA(:,:)
    REAL(KIND=8),ALLOCATABLE::VELON(:,:),VELOT(:,:)

    !---------�����������(������֮һrou*u�������ٻ�---------!
    REAL(KIND=8),ALLOCATABLE::CU(:),CP(:)
    REAL(KIND=8),ALLOCATABLE::CUX(:),CPX(:)
    REAL(KIND=8),ALLOCATABLE::CUY(:),CPY(:)
    REAL(KIND=8)::CUX_TOTAL,CPX_TOTAL
    REAL(KIND=8)::CUY_TOTAL,CPY_TOTAL
    REAL(KIND=8)::CU_TOTAL,CP_TOTAL
    !------��ϵ��------!
    REAL(KIND=8),ALLOCATABLE::CX(:),CY(:)
    REAL(KIND=8),ALLOCATABLE::CXT(:),CYT(:)
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!�������������ϵ�µ���
    REAL(KIND=8)::CXT_TOTAL,CYT_TOTAL!��ʵ��������ϵ�µ���
    REAL(KIND=8)::CXF_TOTAL,CYF_TOTAL!��������ϵ�µ���
    !------��ϵ��------!
    REAL(KIND=8),ALLOCATABLE::CP_PROP(:),CP_ELEV(:),CP_EFFE(:),CP_AERO(:)
    REAL(KIND=8),ALLOCATABLE::EFFICIENCY(:)
    REAL(KIND=8)::CP_PROP_TOTAL,CP_ELEV_TOTAL,CP_EFFE_TOTAL,CP_AERO_TOTAL,ELEV_EFFICIENCY_TOTAL,PROP_EFFICIENCY_TOTAL,EFFICIENCY_TOTAL


    !------��ʼ��------!
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
    RJIM=1!�����غϵ��������

    RDA=2.0D0*PI/DBLE(RIM)
    RDN=DSQRT(2.0D0)*DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��

    ALLOCATE( XR(RIM,RJM),YR(RIM,RJM) )
    ALLOCATE( XA(RIM,RJM),YA(RIM,RJM) )
    ALLOCATE( UR(RIM,RJM),VR(RIM,RJM),PR(RIM,RJM) )
    ALLOCATE( UA(RIM,RJM),VA(RIM,RJM) )
    ALLOCATE( VELON(RIM,RJM),VELOT(RIM,RJM) )
    ALLOCATE( NXR(RIM),NYR(RIM) )
    ALLOCATE( NXA(RIM),NYA(RIM) )
    ALLOCATE( ARC(RIM) )

    !------������������------!
    !�����㼰������
    RJ=RJIM
    DO RI=1,RIM,1
        !��X����������ʱ����ת
        XR(RI,RJ)=LAXIS*DCOS( RDA *(DBLE(RI)-0.5D0) )+CEN_DEVIATION(1)
        YR(RI,RJ)=SAXIS*DSIN( RDA *(DBLE(RI)-0.5D0) )+CEN_DEVIATION(2)
        CALL CRDNT_TRANSFORM(MAT_REL2ABS,XR(RI,RJ),YR(RI,RJ),XATEMP,YATEMP)
        XA(RI,RJ)=XATEMP+CEN_ELP(1)
        YA(RI,RJ)=YATEMP+CEN_ELP(2)
        CALL NORMALVECTOR_2D(ELLIPSE_GEOMETRICAL,XA(RI,RJ),YA(RI,RJ),NXA(RI),NYA(RI))
        CALL CRDNT_TRANSFORM(MAT_ABS2REL,NXA(RI),NYA(RI),NXR(RI),NYR(RI))
    END DO
    !����
    DO RI=1,RIM,1
        XATEMP=LAXIS*DABS( DCOS( RDA * DBLE(RI) )-DCOS( RDA * DBLE(RI-1) ) )
        YATEMP=SAXIS*DABS( DSIN( RDA * DBLE(RI) )-DSIN( RDA * DBLE(RI-1) ) )
        ARC(RI)=DSQRT(XATEMP**2.0D0+YATEMP**2.0D0)
    END DO
    !�����
    RJ=RJM
    DO RI=1,RIM,1
        XA(RI,RJ)=XA(RI,RJIM)+RDN*NXA(RI)
        YA(RI,RJ)=YA(RI,RJIM)+RDN*NYA(RI)
        XR(RI,RJ)=XR(RI,RJIM)+RDN*NXR(RI)
        YR(RI,RJ)=YR(RI,RJIM)+RDN*NYR(RI)
    END DO

    !------����������������ʼ��------!
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

    !------����ֱ�Ӳ�ֵ������з����ٶ�------!
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

    !------���㸳�趯�߽��˶��ٶȲ�����з����ٶ�------!
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

    !------������ϵ���ٶ�------!
    DO RJ=1,RJM,1
        DO RI=1,RIM,1
            CALL CRDNT_TRANSFORM(MAT_ABS2REL,UA(RI,RJ),VA(RI,RJ),UR(RI,RJ),VR(RI,RJ))
        END DO
    END DO


    !METHOD==1
    !------����Ч�ʵļ���------!
    RJ=RJIM
    DO RI=1,RIM,1
        CU(RI)=2.0D0/Re*( VELOT(RI,RJ+1)-VELOT(RI,RJ) )/RDN
        CP(RI)=2.0D0*PR(RI,RJ+1)
        !CU�غ�������������(NYA,-NXA),CP�ط�����������(-NXA,-NYA)
        IF(NYR(RI)>=CRITERIA)THEN
            CUX(RI)= CU(RI)*NYA(RI)
            CUY(RI)=-CU(RI)*NXA(RI)
        ELSE!CU�غ�������������(-NYA,NXA),CP�ط�����������(-NXA,-NYA)
            CUX(RI)=-CU(RI)*NYA(RI)
            CUY(RI)= CU(RI)*NXA(RI)
        END IF
        CPX(RI)=-CP(RI)*NXA(RI)
        CPY(RI)=-CP(RI)*NYA(RI)
        !���ѹ������ճ����
        CUX_TOTAL=CUX_TOTAL+CUX(RI)*ARC(RI)/(2.0D0*LAXIS)
        CPX_TOTAL=CPX_TOTAL+CPX(RI)*ARC(RI)/(2.0D0*LAXIS)

        CUY_TOTAL=CUY_TOTAL+CUY(RI)*ARC(RI)/(2.0D0*LAXIS)
        CPY_TOTAL=CPY_TOTAL+CPY(RI)*ARC(RI)/(2.0D0*LAXIS)
        CU_TOTAL=CU_TOTAL+CU(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_TOTAL=CP_TOTAL+CP(RI)*ARC(RI)/(2.0D0*LAXIS)

        !���������������ϵ��������
        CX(RI)=CUX(RI)+CPX(RI)
        CY(RI)=CUY(RI)+CPY(RI)
        CXC_TOTAL=CXC_TOTAL+CX(RI)*ARC(RI)/(2.0D0*LAXIS)
        CYC_TOTAL=CYC_TOTAL+CY(RI)*ARC(RI)/(2.0D0*LAXIS)

        !�����ʵ����ϵ��������������
        CALL CRDNT_TRANSFORM(MAT_ABS2TRU,CX(RI),CY(RI),CXT(RI),CYT(RI))
        CXT_TOTAL=CXT_TOTAL+CXT(RI)*ARC(RI)/(2.0D0*LAXIS)
        CYT_TOTAL=CYT_TOTAL+CYT(RI)*ARC(RI)/(2.0D0*LAXIS)

        !�����ʵ��������ϵ�����ù�����
        CP_PROP(RI)=CXT(RI)*VELO_RATIO*DCOS(TRUX_FLIGHT_ANGLE)
        CP_ELEV(RI)=CYT(RI)*VELO_RATIO*DSIN(TRUX_FLIGHT_ANGLE)
        CP_EFFE(RI)=CP_PROP(RI)+CP_ELEV(RI)
        !�����ʵ��������ϵ������������
        CALL CRDNT_TRANSFORM(MAT_ABS2TRU,UA(RI,RJ),VA(RI,RJ),UT,VT)
        CP_AERO(RI)=-CXT(RI)*(VELO_RATIO*DCOS(TRUX_FLIGHT_ANGLE)+UT)-CYT(RI)*(VELO_RATIO*DSIN(TRUX_FLIGHT_ANGLE)+VT)
        !Ч�ʷ���
        EFFICIENCY(RI)=CP_EFFE(RI)/CP_AERO(RI)
        !�����ʵ��������ϵ������������
        CP_PROP_TOTAL=CP_PROP_TOTAL+CP_PROP(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_ELEV_TOTAL=CP_ELEV_TOTAL+CP_ELEV(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_EFFE_TOTAL=CP_EFFE_TOTAL+CP_EFFE(RI)*ARC(RI)/(2.0D0*LAXIS)
        CP_AERO_TOTAL=CP_AERO_TOTAL+CP_AERO(RI)*ARC(RI)/(2.0D0*LAXIS)

    END DO

    !�����������ϵ��������
    CALL CRDNT_TRANSFORM(MAT_ABS2FLP,CXC_TOTAL,CYC_TOTAL,CXF_TOTAL,CYF_TOTAL)
    !�����ʵ����ϵ����Ч��
    PROP_EFFICIENCY_TOTAL=CP_PROP_TOTAL/CP_AERO_TOTAL
    ELEV_EFFICIENCY_TOTAL=CP_ELEV_TOTAL/CP_AERO_TOTAL
    EFFICIENCY_TOTAL=CP_EFFE_TOTAL/CP_AERO_TOTAL

    !*************************���ϵ��********************************!
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
    !*************************����������********************************!
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
    !*************************�����������********************************!
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
    !!*************************���ѹ��ϵ��********************************!
    !!IF( MOD(NSTEP,200)==0 )THEN
    !    OPEN(UNIT=10,FILE='P_COEFFICIENT'//TRIM(CHAR_ID)//'_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.DAT')
    !    RJ=RJIM
    !    DO RI=1,RIM,1
    !        WRITE(10,*) DBLE(RI)/DBLE(RIM),PR(RI,RJ+1)
    !    END DO
    !!END IF

    !------У�鷨������------!
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
    !------У���з����ٶȵ�------!
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
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************����ת��******************************************************!
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

    !***************************************************�����Բ�ֵ����U,V,W,P��ĳ����ֵ���������******************************************************!
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


    !***************************************************�������������Ч��******************************************************!
    SUBROUTINE CAL_CLCT_CIRCULAR
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !��������
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CYLINDER_KINETIC(9)
    !��������ϵ�¶����������ѧ���ʽϵ��
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !�����ܶ�
    REAL(KIND=8)::RDC,RDR!���򣬾���
    !����ֵ
    REAL(KIND=8),ALLOCATABLE::XR(:,:),YR(:,:)
    REAL(KIND=8),ALLOCATABLE::NXR(:),NYR(:)
    REAL(KIND=8)::XA,YA
    !������
    INTEGER::RIM
    INTEGER::RJM
    INTEGER::RI,RJ,RJIM
    !---------��ֵ�漰�ľ����������---------!
    !�ǵ���Ϣ
    INTEGER::IAU,JAU!���½ǽű�ֵ
    INTEGER::IAV,JAV!���½ǽű�ֵ
    INTEGER::IAP,JAP!���½ǽű�ֵ

    !�ٶ�ѹ����
    REAL(KIND=8),ALLOCATABLE::UR(:,:),VR(:,:),PR(:,:)
    !REAL(KIND=8)::UA,VA
    REAL(KIND=8),ALLOCATABLE::UA(:,:),VA(:,:)

    !---------�����(������֮һrou*u�������ٻ�---------!
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
    !INTEGER::METHOD!ճ�������㷽����1Ϊ��������ֵ��2Ϊ�����������ֵ

    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------�߽�1-----------------------------------------!
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
        !��ΪCOXY==0��COX2==COY2
        RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

        RJM=3
        RIM=IDNINT( 2.0D0*RADIUS/DX3 )
        RJIM=1!�����غϵ��������

        RDC=2.0D0*RADIUS*PI/DBLE(RIM)
        RDR=DSQRT(2.0D0)*DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��

        !METHOD=2

        WRITE(*,*) "������ϵ�������1��",RIM,RJM

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

        !------������������------!
        DO RJ=1,RJM,1
            DO RI=1,RIM,1
                !��X����ʱ����ת
                XR(RI,RJ)=CEN_CRC(1)+(RADIUS+DBLE(RJ-1)*RDR)*DCOS(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                YR(RI,RJ)=CEN_CRC(2)+(RADIUS+DBLE(RJ-1)*RDR)*DSIN(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                NXR(RI)=DCOS(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
                NYR(RI)=DSIN(2.0D0*PI*DBLE(RI-1)/DBLE(RIM))
            END DO
        END DO

        !------��ʼ��------!
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


        !------����ֱ�Ӳ�ֵ------!
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

        !------���㸳�趯�߽��˶��ٶ�------!
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
        !------����Ч�ʵļ���------!
        RJ=RJIM
        DO RI=1,RIM,1
            !IF(XR(RI)<=RXTRL .AND. XR(RI)>=RXLED)THEN
            CU(RI)=2.0D0/Re*( UR(RI,RJ+2)-UR(RI,RJ+1) )/RDR
            CP(RI)=2.0D0*PR(RI,RJ+1)
            !CU�غ�������������(NYR,-NXR),CP�ط�����������(-NXR,-NYR)
            IF(NYR(RI)>=CRITERIA)THEN
                CUX(RI)= CU(RI)*NYR(RI)
                CUY(RI)=-CU(RI)*NXR(RI)
            ELSE!CU�غ�������������(-NYR,NXR),CP�ط�����������(-NXR,-NYR)
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

        !*************************���ϵ��********************************!
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

        !*************************�������********************************!
        IF( MOD(NSTEP,1000)==0 )THEN
            OPEN(UNIT=10,FILE='SECTION2_Re'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')
            WRITE(10,*) 'TITLE="NONAME"'
            WRITE(10,*) 'VARIABLES="XR","YR","PR","UA","VA","CUX","CUY","CPX","CPY","CX","CY"'
            WRITE(10,*) 'ZONE T="NONAME", I=',RIM,', J=',RJM,', F=POINT'
            RJ=RJIM
            DO RI=1,RIM,1
                WRITE(10,*) XR(RI,RJ),YR(RI,RJ),PR(RI,RJ),UA(RI,RJ),VA(RI,RJ),CUX(RI),CUY(RI),CPX(RI),CPY(RI),CX(RI),CY(RI)
                !!!!!!!!ע��˴�UA UR������
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

    !***************************************************�������������Ч��******************************************************!
    SUBROUTINE CAL_CLCT_2DCURVE
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !---------�����㷨�����ԵĲ���---------!
    !��������
    REAL(KIND=8)::CEN_CRC(2)
    REAL(KIND=8)::RADIUS
    REAL(KIND=8)::CYLINDER_KINETIC(9)
    !��������ϵ�¶����������ѧ���ʽϵ��
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    !---------�����㷨�����ԵĲ���---------!

    !�����ܶ�
    REAL(KIND=8)::DN1!
    REAL(KIND=8),ALLOCATABLE::DN2(:)!
    !����ֵ
    REAL(KIND=8),ALLOCATABLE::XBF(:,:),YBF(:,:)
    REAL(KIND=8),ALLOCATABLE::N1X(:),N1Y(:),N2X(:),N2Y(:)
    REAL(KIND=8)::XA,YA
    !������
    INTEGER::BFIM
    INTEGER::BFJM
    INTEGER::BFI,BFJ
    !��ֵ�漰�ľ�������ǵ���Ϣ
    INTEGER::IAU,JAU!���½ǽű�ֵ
    INTEGER::IAV,JAV!���½ǽű�ֵ
    INTEGER::IAP,JAP!���½ǽű�ֵ

    !�����ٶ�ѹ����
    REAL(KIND=8),ALLOCATABLE::UBFN1(:,:),UBFN2(:,:),PBF(:,:)!��������ϵ�ٶ�
    REAL(KIND=8),ALLOCATABLE::UBF(:,:),VBF(:,:)!��������ϵ�ٶ�
    !����ת������
    REAL(KIND=8)::MAT_ABS2CUR(2,2),MAT_CUR2ABS(2,2)

    !---------��ϵ��---------!
    REAL(KIND=8),ALLOCATABLE::CFV1(:),CFV2(:),CFP(:)
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!�������������ϵ�µ���
    REAL(KIND=8)::CXT_TOTAL,CYT_TOTAL!��ʵ��������ϵ�µ���
    REAL(KIND=8)::CXF_TOTAL,CYF_TOTAL!��������ϵ�µ���

    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------�߽�1-----------------------------------------!
    CYLINDER_KINETIC=QUADRIC_KINETIC_1
    CEN_CRC(1)=QUADRIC_KINETIC_1(4)
    CEN_CRC(2)=QUADRIC_KINETIC_1(5)
    COX2=QUADRIC_GEOMETRICAL_1(7)
    COY2=QUADRIC_GEOMETRICAL_1(8)
    COXY=QUADRIC_GEOMETRICAL_1(9)
    COX=QUADRIC_GEOMETRICAL_1(10)
    COY=QUADRIC_GEOMETRICAL_1(11)
    COM=QUADRIC_GEOMETRICAL_1(13)
    !��ΪCOXY==0��COX2==COY2
    RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

    BFJM=3
    BFIM=IDNINT( 2.0D0*RADIUS/DX3 )

    ALLOCATE( DN2(BFJM) )

    DN1=DSQRT(2.0D0)*DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "������ϵ�������1��",BFIM,BFJM

    ALLOCATE( XBF(BFIM,BFJM),YBF(BFIM,BFJM) )
    ALLOCATE( N1X(BFIM),N1Y(BFIM),N2X(BFIM),N2Y(BFIM) )
    ALLOCATE( UBFN1(BFIM,BFJM),UBFN2(BFIM,BFJM),PBF(BFIM,BFJM) )
    ALLOCATE( UBF(BFIM,BFJM),VBF(BFIM,BFJM) )
    ALLOCATE( CFV1(BFIM),CFV2(BFIM),CFP(BFIM) )

    !------������������------!
    DO BFJ=1,BFJM,1
        DO BFI=1,BFIM,1
            !��X����ʱ����ת
            XBF(BFI,BFJ)=CEN_CRC(1)+(RADIUS+DBLE(BFJM-BFJ)*DN1)*DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            YBF(BFI,BFJ)=CEN_CRC(2)+(RADIUS+DBLE(BFJM-BFJ)*DN1)*DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N1X(BFI)=-DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N1Y(BFI)=-DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N2X(BFI)= DSIN(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
            N2Y(BFI)=-DCOS(-2.0D0*PI*DBLE(BFI-1)/DBLE(BFIM))
        END DO
    END DO

    !------��ʼ��------!
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

    !------����ֱ�Ӳ�ֵ------!
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

    !------���㸳�趯�߽��˶��ٶ�------!
    BFJ=BFJM
    DO BFI=1,BFIM,1
        XA=XBF(BFI,BFJ)
        YA=YBF(BFI,BFJ)
        CALL VELOCITY_LB(CYLINDER_KINETIC,XA,YA,UBF(BFI,BFJ),VBF(BFI,BFJ))
    END DO

    !------�ٶȳ�����任------!
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

    !------���������------!
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

    WRITE(30009,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE

    !***************************************************������������������巽��,�򻯣�******************************************************!
    SUBROUTINE CAL_CLCT_CONTROL_VOLUME_SIMPLIFIED
    USE DECLARATION
    IMPLICIT NONE

    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !����������
    REAL(KIND=8)::CV_LEFT,CV_RIGH,CV_BOTT,CV_TOPP

    !Ѱַ
    REAL(KIND=8),ALLOCATABLE::DISTANCE_X(:)
    INTEGER::INDEX_XL_U,INDEX_XR_U
    INTEGER::INDEX_XL_P,INDEX_XR_P

    !������Ϣ
    REAL(KIND=8),ALLOCATABLE::U_LEFT_BOUNDARY(:),U_RIGH_BOUNDARY(:)
    REAL(KIND=8),ALLOCATABLE::P_LEFT_BOUNDARY(:),P_RIGH_BOUNDARY(:)
    REAL(KIND=8),ALLOCATABLE::BOUNDARY_AREA(:)

    !------��ϵ��------!
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!�������������ϵ�µ���


    ALLOCATE( BOUNDARY_AREA(JM-1) )

    ALLOCATE( DISTANCE_X(IM) )
    ALLOCATE( U_LEFT_BOUNDARY(JM-1),U_RIGH_BOUNDARY(JM-1) )
    ALLOCATE( P_LEFT_BOUNDARY(JM-1),P_RIGH_BOUNDARY(JM-1) )


    CV_LEFT=-1.0D0
    CV_RIGH= 1.0D0
    CV_BOTT=BOTT
    CV_TOPP=TOPP


    !����������
    DISTANCE_X=X(:)-CV_LEFT
    INDEX_XL_U=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( X(INDEX_XL_U)>CV_LEFT ) INDEX_XL_U=INDEX_XL_U-1
    DISTANCE_X=X(:)-CV_RIGH
    INDEX_XR_U=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( X(INDEX_XR_U)>CV_RIGH ) INDEX_XR_U=INDEX_XR_U-1
    !����������
    DISTANCE_X=XPV(1:IM)-CV_LEFT
    INDEX_XL_P=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( XPV(INDEX_XL_P)>CV_LEFT ) INDEX_XL_P=INDEX_XL_P-1
    DISTANCE_X=XPV(1:IM)-CV_RIGH
    INDEX_XR_P=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( XPV(INDEX_XR_P)>CV_RIGH ) INDEX_XR_P=INDEX_XR_P-1

    !��ֵ�ٶȺ�ѹ��
    DO J=1,JM-1,1
        CALL LINEAR_INTERPOLATION(U(INDEX_XL_U,J),U_LEFT_BOUNDARY(J),U(INDEX_XL_U+1,J),X  (INDEX_XL_U),CV_LEFT,X  (INDEX_XL_U+1))
        CALL LINEAR_INTERPOLATION(U(INDEX_XR_U,J),U_RIGH_BOUNDARY(J),U(INDEX_XR_U+1,J),X  (INDEX_XR_U),CV_RIGH,X  (INDEX_XR_U+1))
        CALL LINEAR_INTERPOLATION(P(INDEX_XL_P,J),P_LEFT_BOUNDARY(J),P(INDEX_XL_P+1,J),XPV(INDEX_XL_P),CV_LEFT,XPV(INDEX_XL_P+1))
        CALL LINEAR_INTERPOLATION(P(INDEX_XR_P,J),P_RIGH_BOUNDARY(J),P(INDEX_XR_P+1,J),XPV(INDEX_XR_P),CV_RIGH,XPV(INDEX_XR_P+1))
    END DO

    !------��߽��С------!
    DO J=1,JM-1,1
        BOUNDARY_AREA(J)=Y(J+1)-Y(J)
    END DO

    CXC_TOTAL=-2.0D0*SUM((U_RIGH_BOUNDARY*U_RIGH_BOUNDARY+P_RIGH_BOUNDARY)*BOUNDARY_AREA-(U_LEFT_BOUNDARY*U_LEFT_BOUNDARY+P_LEFT_BOUNDARY)*BOUNDARY_AREA)


    WRITE(30001,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !WRITE(30001,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE


    !***************************************************������������������巽����******************************************************!
    SUBROUTINE CAL_CLCT_CONTROL_VOLUME
    USE DECLARATION
    IMPLICIT NONE

    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !����������
    REAL(KIND=8)::CV_LEFT,CV_RIGH,CV_BOTT,CV_TOPP

    !Ѱַ
    REAL(KIND=8),ALLOCATABLE::DISTANCE_X(:),DISTANCE_Y(:)
    INTEGER::INDEX_XL,INDEX_XR
    INTEGER::INDEX_YB,INDEX_YT

    !������Ϣ
    !REAL(KIND=8),ALLOCATABLE::U_LEFT_BOUNDARY(:),U_RIGH_BOUNDARY(:)
    !REAL(KIND=8),ALLOCATABLE::P_LEFT_BOUNDARY(:),P_RIGH_BOUNDARY(:)
    !REAL(KIND=8),ALLOCATABLE::BOUNDARY_AREA(:)

    !------��ϵ��------!
    !�������������ϵ�µ���
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL

    !�Ƕ���
    REAL(KIND=8)::CFTX,CFTY

    !����
    REAL(KIND=8)::CFCLX,CFCLY
    REAL(KIND=8)::CFCRX,CFCRY
    REAL(KIND=8)::CFCBX,CFCBY
    REAL(KIND=8)::CFCTX,CFCTY

    !ѹ��
    REAL(KIND=8)::CFPL
    REAL(KIND=8)::CFPR
    REAL(KIND=8)::CFPB
    REAL(KIND=8)::CFPT

    !ճ��Ӧ��
    REAL(KIND=8)::CFVLX,CFVLY
    REAL(KIND=8)::CFVRX,CFVRY
    REAL(KIND=8)::CFVBX,CFVBY
    REAL(KIND=8)::CFVTX,CFVTY

    !����
    REAL(KIND=8)::CFCX,CFCY
    REAL(KIND=8)::CFPX,CFPY
    REAL(KIND=8)::CFVX,CFVY


    ALLOCATE( DISTANCE_X(IM),DISTANCE_Y(JM) )

    !!2
    !CV_LEFT=-0.7D0
    !CV_RIGH= 0.7D0
    !CV_BOTT=-0.7D0
    !CV_TOPP= 0.7D0
    !
    !!����������
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!����ʼ��
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!����������
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!���Ƕ���
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!����
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!�������������ϵ�µ���
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30002,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL

    !!3
    !CV_LEFT=-0.6D0
    !CV_RIGH= 0.9D0
    !CV_BOTT=-0.6D0
    !CV_TOPP= 0.9D0
    !
    !!����������
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!����ʼ��
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!����������
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!���Ƕ���
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!����
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!�������������ϵ�µ���
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30003,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !
    !!4
    !CV_LEFT=-0.8D0
    !CV_RIGH= 0.8D0
    !CV_BOTT=-0.8D0
    !CV_TOPP= 0.8D0
    !
    !!����������
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!����ʼ��
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!����������
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!���Ƕ���
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!����
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!�������������ϵ�µ���
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30004,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !
    !!5
    !CV_LEFT=-0.9D0
    !CV_RIGH= 0.9D0
    !CV_BOTT=-0.9D0
    !CV_TOPP= 0.9D0
    !
    !!����������
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!����ʼ��
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!����������
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!���Ƕ���
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!����
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!�������������ϵ�µ���
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30005,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL

    !6
    IF(IB_LOCOMOTION==-1)THEN
        CV_LEFT=-0.6D0
        CV_RIGH= 0.6D0
        CV_BOTT=-0.6D0
        CV_TOPP= 0.6D0
    ELSE
        CV_LEFT=LEIN+0.1D0
        CV_RIGH=RIIN-0.1D0
        CV_BOTT=BOIN+0.1D0
        CV_TOPP=TOIN-0.1D0
    END IF

    !����������
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
    IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"

    !����ʼ��
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

    !����������
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
    !���Ƕ���
    DO J=INDEX_YB,INDEX_YT-1,1
        DO I=INDEX_XL,INDEX_XR-1,1
            CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
            CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
        END DO
    END DO

    !����
    CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    CFPX=-CFPL+CFPR
    CFPY=-CFPB+CFPT
    CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    CFVY=-CFVLY+CFVRY-CFVBY+CFVTY

    !�������������ϵ�µ���
    CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)

    WRITE(30006,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    WRITE(30007,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE