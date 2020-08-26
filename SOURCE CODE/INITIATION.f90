    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !****************************************************初始化流场和时间层跨度*********************************************************!
    SUBROUTINE INITIATION
    USE DECLARATION
    !USE QUADRIC_PARAMETER
    IMPLICIT NONE
    CHARACTER(LEN=80)::FILENAME="插值数据-点.txt"
    INTEGER::TASK
    LOGICAL ALIVE
    INTEGER::ERROR

    WRITE(*,*)"输入数字确定初始化方式：1.赋初值 2.续算"
    READ(*,*)TASK

    IF(TASK==1)THEN
        !$!
        U=U_FREESTREAM
        V=V_FREESTREAM
        UN=U_FREESTREAM
        VN=V_FREESTREAM
        UN1=U_FREESTREAM
        VN1=V_FREESTREAM
        UHAT=U_FREESTREAM
        VHAT=V_FREESTREAM
        P=0.0D0
        !IF(IB_LOCOMOTION==-1)THEN
        !    CALL CYLINDER_POTENTIAL_FLOW_INITIATION
        !END IF
        IF(IB_LOCOMOTION==-2)THEN
            CALL ANALYTICAL_CASE_INITIATION
            !CALL OUTPUT_PLT_1_STAGGERED
            !CALL OUTPUT_FULL_STAGGERED
        END IF
        PHI=0.0D0
        WRITE(*,*)"输入起始时间步(解析解算例为1，其余为0)："
        READ(*,*)NSTART
        WRITE(*,*)"赋初值成功,起始时间步：",NSTART
    ELSE IF(TASK==2)THEN
        CALL READFROMFILE_2D_STAGGERED(FILENAME_RESTART)
        UN=U
        VN=V
        UN1=U
        VN1=V
        UHAT=U
        VHAT=V
        PHI=0.0D0
        WRITE(*,*)"输入起始时间步："
        READ(*,*)NSTART
        WRITE(*,*)"成功读入，开始续算，起始时间步：",NSTART
    END IF

    NMAX=NSTART+NDURATION

    RETURN
    END


    !****************************************************圆柱绕流势流初始化*********************************************************!
    SUBROUTINE CYLINDER_POTENTIAL_FLOW_INITIATION
    USE DECLARATION
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    REAL(KIND=8)::R2,SINT,COST,VR,VT,CEN_INI(2)

    CEN_INI(1)=QUADRIC_KINETIC_1(4)
    CEN_INI(2)=QUADRIC_KINETIC_1(5)

    !------U------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            R2=( X(I)-CEN_INI(1) )**2.0D0+( YPU(J)-CEN_INI(2) )**2.0D0
            IF( R2<0.25D0 )THEN
                U(I,J)=0.0D0
                UN(I,J)=0.0D0
                UN1(I,J)=0.0D0
                UHAT(I,J)=0.0D0
            ELSE
                COST=X(I)/DSQRT(R2)
                SINT=YPU(J)/DSQRT(R2)

                VR= U_FREESTREAM*COST*(1-0.25D0/R2)
                VT=-U_FREESTREAM*SINT*(1+0.25D0/R2)

                U(I,J)=VR*COST-VT*SINT
                UN(I,J)=U(I,J)
                UN1(I,J)=U(I,J)
                UHAT(I,J)=U(I,J)
            END IF
        END DO
    END DO

    !------施加边界条件------!
    U(1       ,1:JM-1:1)=BCU_AL*U(2       ,1:JM-1:1)+BCU_BL
    U(IM      ,1:JM-1:1)=BCU_AR*U(IM-1    ,1:JM-1:1)+BCU_BR
    U(:       ,0       )=BCU_AB*U(:       ,1       )+BCU_BB
    U(:       ,JM      )=BCU_AT*U(:       ,JM-1    )+BCU_BT

    !------V------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            R2=( XPV(I)-CEN_INI(1) )**2.0D0+( Y(J)-CEN_INI(1) )**2.0D0
            IF( R2<0.25D0 )THEN
                V(I,J)=0.0D0
                V(I,J)=0.0D0
                VN1(I,J)=0.0D0
                VHAT(I,J)=0.0D0
            ELSE
                COST=XPV(I)/DSQRT(R2)
                SINT=Y(J)/DSQRT(R2)

                VR= U_FREESTREAM*COST*(1-0.25D0/R2)
                VT=-U_FREESTREAM*SINT*(1+0.25D0/R2)

                V(I,J)=VR*SINT+VT*COST
                VN(I,J)=V(I,J)
                VN1(I,J)=V(I,J)
                VHAT(I,J)=V(I,J)
            END IF
        END DO
    END DO

    !------施加边界条件------!
    V(1:IM-1:1,1       )=BCV_AB*V(1:IM-1:1,2       )+BCV_BB
    V(1:IM-1:1,JM      )=BCV_AT*V(1:IM-1:1,JM-1    )+BCV_BT
    V(0       ,:       )=BCV_AL*V(1       ,:       )+BCV_BL
    V(IM      ,:       )=BCV_AR*V(IM-1    ,:       )+BCV_BR


    RETURN
    END