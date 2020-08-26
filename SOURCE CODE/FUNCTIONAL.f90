    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !******************************************求解二次曲面边界格点外法向量*****************************************************!
    SUBROUTINE NORMALVECTOR(QUADRIC_GEOMETRICAL,X,Y,Z,NX,NY,NZ)
    !在其他未知曲面方程的情况下使用九点拟合的二次曲面来求解法向量
    !NX=Fx(X,Y,Z),NY=Fy(X,Y,Z),NZ=Fz(X,Y,Z)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!注意标号
    IMPLICIT NONE

    REAL(KIND=8)::X,Y,Z,NX,NY,NZ,NORM
    REAL(KIND=8)::QUADRIC_GEOMETRICAL(36)
    REAL(KIND=8)::COX2,COY2,COZ2,COXY,COXZ,COYZ,COX,COY,COZ

    COX2=QUADRIC_GEOMETRICAL(10)
    COY2=QUADRIC_GEOMETRICAL(11)
    COZ2=QUADRIC_GEOMETRICAL(12)
    COXY=QUADRIC_GEOMETRICAL(13)
    COXZ=QUADRIC_GEOMETRICAL(14)
    COYZ=QUADRIC_GEOMETRICAL(15)
    COX=QUADRIC_GEOMETRICAL(16)
    COY=QUADRIC_GEOMETRICAL(17)
    COZ=QUADRIC_GEOMETRICAL(18)

    NX=2.0D0*COX2*X+COXY*Y+COXZ*Z+COX
    NY=2.0D0*COY2*Y+COXY*X+COYZ*Z+COY
    NZ=2.0D0*COZ2*Z+COXZ*X+COYZ*Y+COZ

    NORM=DSQRT(NX**2.0D0+NY**2.0D0+NZ**2.0D0)
    NX=NX/NORM
    NY=NY/NORM
    NZ=NZ/NORM

    RETURN
    END SUBROUTINE
    
    !******************************************求解二次曲线边界格点外法向量*****************************************************!
    SUBROUTINE NORMALVECTOR_2D(QUADRIC_GEOMETRICAL,X,Y,NX,NY)
    !在其他未知曲面方程的情况下使用九点拟合的二次曲面来求解法向量
    !NX=Fx(X,Y,Z),NY=Fy(X,Y,Z),NZ=Fz(X,Y,Z)
    IMPLICIT NONE

    REAL(KIND=8)::X,Y,NX,NY,NORM
    REAL(KIND=8)::QUADRIC_GEOMETRICAL(36)
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY

    COX2=QUADRIC_GEOMETRICAL(7)
    COY2=QUADRIC_GEOMETRICAL(8)
    COXY=QUADRIC_GEOMETRICAL(9)
    COX =QUADRIC_GEOMETRICAL(10)
    COY =QUADRIC_GEOMETRICAL(11)

    NX=2.0D0*COX2*X+COXY*Y+COX
    NY=2.0D0*COY2*Y+COXY*X+COY

    NORM=DSQRT(NX**2.0D0+NY**2.0D0)
    NX=NX/NORM
    NY=NY/NORM

    RETURN
    END SUBROUTINE

    !************************************************冒泡法排序*************************************************************!
    SUBROUTINE BUBBLESORT(A,M,N,NS)
    IMPLICIT NONE

    INTEGER::M,N,NS
    INTEGER I,J
    REAL(KIND=8)::A(M,N),TEMP(N)
    DO I=M-1,1,-1
        DO J=1,I
            IF(A(J,NS)>A(J+1,NS))THEN
                TEMP=A(J,:)
                A(J,:)=A(J+1,:)
                A(J+1,:)=TEMP
            END IF
        END DO
    END DO

    RETURN
    END SUBROUTINE

    !****************************************求解运动边界上格点的速度信息***************************************************!
    !------待改进------!
    SUBROUTINE VELOCITY_LB(BOUNDARY_KINETIC,X,Y,U,V)
    IMPLICIT NONE
    REAL(KIND=8)::CENB(2),TRANB(2),OMEGAB,R(2),TEMP(2)
    REAL(KIND=8)::X,Y,U,V
    REAL(KIND=8)::BOUNDARY_KINETIC(9)

    TRANB(1)=BOUNDARY_KINETIC(1)
    TRANB(2)=BOUNDARY_KINETIC(2)
    CENB(1)=BOUNDARY_KINETIC(4)
    CENB(2)=BOUNDARY_KINETIC(5)
    OMEGAB=BOUNDARY_KINETIC(7)

    R(1)=X-CENB(1)
    R(2)=Y-CENB(2)
    !------V=W叉乘R------!
    TEMP(1)=-OMEGAB*R(2)!OMEGAB(2)*R(3)
    TEMP(2)=OMEGAB*R(1)!-OMEGAB(1)*R(3)
    !TEMP(3)=OMEGAB(1)*R(2)-OMEGAB(2)*R(1)
    !------坐标系变换------!
    U=TEMP(1)+TRANB(1)
    V=TEMP(2)+TRANB(2)

    RETURN
    END SUBROUTINE

    !************************************************返回某个量的正负号*************************************************************!
    SUBROUTINE GETSIGN(X,Y)
    IMPLICIT NONE
    REAL(KIND=8)::X
    INTEGER::Y

    IF(X>=0)THEN
        Y=1.0D0
    ELSE
        Y=-1.0D0
    END IF

    RETURN
    END SUBROUTINE

    !************************************************给出转换矩阵*************************************************************!
    SUBROUTINE CAL_TRANMAT(ANGLE,MATRIX)
    IMPLICIT NONE
    REAL(KIND=8)::ANGLE
    INTEGER::AXIS
    REAL(KIND=8)::MATRIX(2,2)

    MATRIX=0.0D0
    MATRIX(1,1)=DCOS(ANGLE)
    MATRIX(1,2)=DSIN(ANGLE)
    MATRIX(2,1)=-DSIN(ANGLE)
    MATRIX(2,2)=DCOS(ANGLE)


    RETURN
    END SUBROUTINE

    !************************************************线性插值*************************************************************!
    SUBROUTINE LINEAR_INTERPOLATION(U2,U1,U0,X2,X1,X0)
    IMPLICIT NONE
    REAL(KIND=8)::U2,U1,U0,X2,X1,X0

    U1=(U2-U0)*(X1-X0)/(X2-X0)+U0

    RETURN
    END SUBROUTINE

    !************************************************确定一点从属的动边界*************************************************************!
    SUBROUTINE DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,VALUE_BOUNDARY_ID)
    USE QUADRIC_PARAMETER
    IMPLICIT NONE
    REAL(KIND=8)::COX2_1,COY2_1,COXY_1,COX_1,COY_1,COM_1!,COZ2,COXZ,COYZ,COZ
    REAL(KIND=8)::COX2_2,COY2_2,COXY_2,COX_2,COY_2,COM_2!,COZ2,COXZ,COYZ,COZ
    REAL(KIND=8)::QUADRIC_VALUE_1,QUADRIC_VALUE_2
    REAL(KIND=8)::VALUE_X,VALUE_Y
    INTEGER::VALUE_BOUNDARY_ID

    COX2_1=QUADRIC_GEOMETRICAL_N1(7)
    COY2_1=QUADRIC_GEOMETRICAL_N1(8)
    COXY_1=QUADRIC_GEOMETRICAL_N1(9)
    COX_1=QUADRIC_GEOMETRICAL_N1(10)
    COY_1=QUADRIC_GEOMETRICAL_N1(11)
    COM_1=QUADRIC_GEOMETRICAL_N1(13)

    COX2_2=QUADRIC_GEOMETRICAL_N2(7)
    COY2_2=QUADRIC_GEOMETRICAL_N2(8)
    COXY_2=QUADRIC_GEOMETRICAL_N2(9)
    COX_2=QUADRIC_GEOMETRICAL_N2(10)
    COY_2=QUADRIC_GEOMETRICAL_N2(11)
    COM_2=QUADRIC_GEOMETRICAL_N2(13)

    QUADRIC_VALUE_1=COX2_1*VALUE_X**2.0D0+COY2_1*VALUE_Y**2.0D0+COXY_1*VALUE_X*VALUE_Y+COX_1*VALUE_X+COY_1*VALUE_Y+COM_1
    QUADRIC_VALUE_2=COX2_2*VALUE_X**2.0D0+COY2_2*VALUE_Y**2.0D0+COXY_2*VALUE_X*VALUE_Y+COX_2*VALUE_X+COY_2*VALUE_Y+COM_2

    IF( QUADRIC_VALUE_1<0.0D0 .AND. QUADRIC_VALUE_2<0.0D0 )THEN
        WRITE(*,*)"ERROR:",VALUE_X,VALUE_Y
    ELSE IF( BOUNDARY_EXISTENCE_1==1 .AND. QUADRIC_VALUE_1<0.0D0 )THEN
        VALUE_BOUNDARY_ID=1
    ELSE IF( BOUNDARY_EXISTENCE_2==1 .AND. QUADRIC_VALUE_2<0.0D0 )THEN
        VALUE_BOUNDARY_ID=2
    ELSE
        VALUE_BOUNDARY_ID=0
    END IF

    RETURN
    END SUBROUTINE