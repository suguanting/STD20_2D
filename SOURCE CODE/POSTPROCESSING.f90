    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                              ��������                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************����ģ��******************************************************!
    SUBROUTINE POSTPROCESSING
    USE DECLARATION
    IMPLICIT NONE
    INTEGER::NPDIS
    LOGICAL::PDIS_PROCESSED=.FALSE.

    NPDIS=NCYCLE
    IF(IB_LOCOMOTION==-1)THEN
        !------���Բ��ѹ���ֲ�------!
        IF( MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPDIS) ) )==0 ) CALL OUTPUT_PRESSURE_DISTRIBUTION(PDIS_PROCESSED)
    END IF

    RETURN
    END SUBROUTINE

    !***************************************************����Բ������ѹ��ϵ��******************************************************!
    SUBROUTINE OUTPUT_PRESSURE_DISTRIBUTION(PROCESSED)
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    LOGICAL::PROCESSED

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

    !---------ѹ��ϵ��---------!
    REAL(KIND=8),ALLOCATABLE::PRESSURE_COEFFICIENT(:)
    REAL(KIND=8)::PRESSURE_FREESTREAM

    !---------������---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------�߽�1-----------------------------------------!
    !CYLINDER_KINETIC=QUADRIC_KINETIC_1
    !CEN_CRC(1)=QUADRIC_KINETIC_1(4)
    !CEN_CRC(2)=QUADRIC_KINETIC_1(5)
    !COX2=QUADRIC_GEOMETRICAL_1(7)
    !COY2=QUADRIC_GEOMETRICAL_1(8)
    !COXY=QUADRIC_GEOMETRICAL_1(9)
    !COX=QUADRIC_GEOMETRICAL_1(10)
    !COY=QUADRIC_GEOMETRICAL_1(11)
    !COM=QUADRIC_GEOMETRICAL_1(13)
    !!��ΪCOXY==0��COX2==COY2
    !RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

    CEN_CRC(1)=0.0D0
    CEN_CRC(2)=0.0D0
    RADIUS=0.5D0

    BFJM=3
    BFIM=IDNINT( 2.0D0*RADIUS/DX3 )


    ALLOCATE( DN2(BFJM) )

    DN1=DSQRT(2.0D0)*DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "ѹ��ϵ�������1��",BFIM,BFJM

    ALLOCATE( XBF(BFIM,BFJM),YBF(BFIM,BFJM) )
    ALLOCATE( N1X(BFIM),N1Y(BFIM),N2X(BFIM),N2Y(BFIM) )
    ALLOCATE( UBFN1(BFIM,BFJM),UBFN2(BFIM,BFJM),PBF(BFIM,BFJM) )
    ALLOCATE( UBF(BFIM,BFJM),VBF(BFIM,BFJM) )
    ALLOCATE( PRESSURE_COEFFICIENT(BFIM) )

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

    !------���ѹ��ϵ��------!
    PRESSURE_FREESTREAM=0.0D0
    PRESSURE_COEFFICIENT=0.0D0
    BFJ=BFJM-1
    DO BFI=1,BFIM,1
        PRESSURE_COEFFICIENT(BFI) =2.0D0*(PBF(BFI,BFJ)-PRESSURE_FREESTREAM)
    END DO


    !DO BFI=1,BFIM,1
    !    WRITE(30,"( I6,2(1X,F9.5))")IDNINT(DBLE(BFI-1)/DBLE(BFIM)*360.0D0-180.0D0),PRESSURE_COEFFICIENT(BFI)
    !END DO
    !WRITE(30,"( I6,2(1X,F9.5))")IDNINT(180.0D0),PRESSURE_COEFFICIENT(1)

    IF(.NOT. PROCESSED)THEN
        DO BFI=1,BFIM,1
            WRITE(501,"(F10.5)") DBLE(BFI-1)/DBLE(BFIM)*360.0D0-180.0D0
        END DO
        WRITE(501,"(F10.5)") 180.0D0
        PROCESSED=.TRUE.
    END IF

    DO BFI=1,BFIM,1
        WRITE(500,"( 1X,F9.5)")PRESSURE_COEFFICIENT(BFI)
    END DO
    WRITE(500,"( 1X,F9.5)")PRESSURE_COEFFICIENT(1)

    !---------------------------------����2---------------------------------------------
    DN1=DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "ѹ��ϵ�������1��",BFIM,BFJM

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

    !------���ѹ��ϵ��------!
    PRESSURE_FREESTREAM=0.0D0
    PRESSURE_COEFFICIENT=0.0D0
    BFJ=BFJM-1
    DO BFI=1,BFIM,1
        PRESSURE_COEFFICIENT(BFI) =2.0D0*(PBF(BFI,BFJ)-PRESSURE_FREESTREAM)
    END DO


    !DO BFI=1,BFIM,1
    !    WRITE(30,"( I6,2(1X,F9.5))")IDNINT(DBLE(BFI-1)/DBLE(BFIM)*360.0D0-180.0D0),PRESSURE_COEFFICIENT(BFI)
    !END DO
    !WRITE(30,"( I6,2(1X,F9.5))")IDNINT(180.0D0),PRESSURE_COEFFICIENT(1)

    DO BFI=1,BFIM,1
        WRITE(502,"( 1X,F9.5)")PRESSURE_COEFFICIENT(BFI)
    END DO
    WRITE(502,"( 1X,F9.5)")PRESSURE_COEFFICIENT(1)

    !---------------------------------����3---------------------------------------------
    DN1=0.5D0*DX3!yֵ��Сֵ����֤���������Ķ����Բ�ֵֻ����������Ӱ��
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "ѹ��ϵ�������1��",BFIM,BFJM

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

    !------���ѹ��ϵ��------!
    PRESSURE_FREESTREAM=0.0D0
    PRESSURE_COEFFICIENT=0.0D0
    BFJ=BFJM-1
    DO BFI=1,BFIM,1
        PRESSURE_COEFFICIENT(BFI) =2.0D0*(PBF(BFI,BFJ)-PRESSURE_FREESTREAM)
    END DO


    !DO BFI=1,BFIM,1
    !    WRITE(30,"( I6,2(1X,F9.5))")IDNINT(DBLE(BFI-1)/DBLE(BFIM)*360.0D0-180.0D0),PRESSURE_COEFFICIENT(BFI)
    !END DO
    !WRITE(30,"( I6,2(1X,F9.5))")IDNINT(180.0D0),PRESSURE_COEFFICIENT(1)

    DO BFI=1,BFIM,1
        WRITE(503,"( 1X,F9.5)")PRESSURE_COEFFICIENT(BFI)
    END DO
    WRITE(503,"( 1X,F9.5)")PRESSURE_COEFFICIENT(1)


    RETURN
    END SUBROUTINE
