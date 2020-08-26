    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                              流场后处理                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************后处理模块******************************************************!
    SUBROUTINE POSTPROCESSING
    USE DECLARATION
    IMPLICIT NONE
    INTEGER::NPDIS
    LOGICAL::PDIS_PROCESSED=.FALSE.

    NPDIS=NCYCLE
    !------输出圆柱压力分布------!
    IF( MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPDIS) ) )==0 ) CALL OUTPUT_PRESSURE_DISTRIBUTION(PDIS_PROCESSED)


    RETURN
    END SUBROUTINE

    !***************************************************后处理圆柱表面压力系数******************************************************!
    SUBROUTINE OUTPUT_PRESSURE_DISTRIBUTION(PROCESSED)
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    LOGICAL::PROCESSED

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

    !---------压力系数---------!
    REAL(KIND=8),ALLOCATABLE::PRESSURE_COEFFICIENT(:)
    REAL(KIND=8)::PRESSURE_FREESTREAM

    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !----------------------------------边界1-----------------------------------------!
    !CYLINDER_KINETIC=QUADRIC_KINETIC_1
    !CEN_CRC(1)=QUADRIC_KINETIC_1(4)
    !CEN_CRC(2)=QUADRIC_KINETIC_1(5)
    !COX2=QUADRIC_GEOMETRICAL_1(7)
    !COY2=QUADRIC_GEOMETRICAL_1(8)
    !COXY=QUADRIC_GEOMETRICAL_1(9)
    !COX=QUADRIC_GEOMETRICAL_1(10)
    !COY=QUADRIC_GEOMETRICAL_1(11)
    !COM=QUADRIC_GEOMETRICAL_1(13)
    !!认为COXY==0，COX2==COY2
    !RADIUS=((0.25D0*COX**2.0D0/COX2+0.25D0*COY**2.0D0/COY2-COM)/COX2)**0.5D0

    CEN_CRC(1)=0.0D0
    CEN_CRC(2)=0.0D0
    RADIUS=0.5D0

    BFJM=3
    BFIM=IDNINT( 2.0D0*RADIUS/DX3 )


    ALLOCATE( DN2(BFJM) )

    DN1=DSQRT(2.0D0)*DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "压力系数格点数1：",BFIM,BFJM

    ALLOCATE( XBF(BFIM,BFJM),YBF(BFIM,BFJM) )
    ALLOCATE( N1X(BFIM),N1Y(BFIM),N2X(BFIM),N2Y(BFIM) )
    ALLOCATE( UBFN1(BFIM,BFJM),UBFN2(BFIM,BFJM),PBF(BFIM,BFJM) )
    ALLOCATE( UBF(BFIM,BFJM),VBF(BFIM,BFJM) )
    ALLOCATE( PRESSURE_COEFFICIENT(BFIM) )

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

    !------求解压力系数------!
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
    
    !---------------------------------试验2---------------------------------------------
    DN1=DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "压力系数格点数1：",BFIM,BFJM

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

    !------求解压力系数------!
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
    
    !---------------------------------试验3---------------------------------------------
    DN1=0.5D0*DX3!y值最小值，保证流场变量的二线性插值只受流场本侧影响
    DO BFJ=1,BFJM,1
        DN2(BFJ)=2.0D0*(RADIUS+DBLE(BFJM-BFJ)*DN1)*PI/DBLE(BFIM)
    END DO

    WRITE(*,*) "压力系数格点数1：",BFIM,BFJM

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

    !------求解压力系数------!
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
