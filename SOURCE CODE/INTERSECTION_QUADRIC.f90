    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !*************************************求解二次曲面交点坐标,外法向量,速度,面积信息******************************************!
    SUBROUTINE INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL,QUADRIC_KINETIC)
    USE QUADRIC_PARAMETER
    USE IMMERSED_BOUNDARY
    USE DECLARATION
    IMPLICIT NONE

    !------几何参数------!
    REAL(KIND=8)::QUADRIC_GEOMETRICAL(36)
    REAL(KIND=8)::QUADRIC_KINETIC(9)

    !坐标转换矩阵
    REAL(KIND=8)::T11,T12,T21,T22

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ

    !------非无限线段端点绝对坐标系限定条件------!
    !上表面
    REAL(KIND=8)::XMAXT,XMINT
    REAL(KIND=8)::YMAXT,YMINT
    !下表面
    REAL(KIND=8)::XMAXB,XMINB
    REAL(KIND=8)::YMAXB,YMINB

    !------前后缘相对坐标系限定条件------!
    REAL(KIND=8)::RXTRL,RXLED
    REAL(KIND=8)::RYTOP,RYBOT

    !------转动中心绝对坐标------!
    REAL(KIND=8)::CEN_INT(2)

    !------法向量------!
    REAL(KIND=8)::NXT,NYT,NXB,NYB

    !------二次方程系数和解------!
    REAL(KIND=8)COEQ2,COEQ1,COEQ0,DETA2,SOLUTION1,SOLUTION2,SOLUTION

    !------临时数据定义及分配------!
    REAL(KIND=8)::UBTEMP,VBTEMP
    REAL(KIND=8),ALLOCATABLE::TEMPS(:,:,:),TEMPE(:,:,:)!储存相应类别交点的临时信息
    INTEGER::TEMP_CHORD

    ALLOCATE( TEMPS(3,IM*JM,3),TEMPE(2,IM*JM,3) )!FOTRAN: COLUMN MAJOR

    !------读取几何参数------!
    T11=QUADRIC_GEOMETRICAL(1)
    T12=QUADRIC_GEOMETRICAL(2)
    T21=QUADRIC_GEOMETRICAL(4)
    T22=QUADRIC_GEOMETRICAL(5)

    COX2=QUADRIC_GEOMETRICAL(7)
    COY2=QUADRIC_GEOMETRICAL(8)
    COXY=QUADRIC_GEOMETRICAL(9)
    COX=QUADRIC_GEOMETRICAL(10)
    COY=QUADRIC_GEOMETRICAL(11)
    COM=QUADRIC_GEOMETRICAL(13)

    XMAXT=QUADRIC_GEOMETRICAL(15)
    XMINT=QUADRIC_GEOMETRICAL(16)
    YMAXT=QUADRIC_GEOMETRICAL(17)
    YMINT=QUADRIC_GEOMETRICAL(18)

    XMAXB=QUADRIC_GEOMETRICAL(21)
    XMINB=QUADRIC_GEOMETRICAL(22)
    YMAXB=QUADRIC_GEOMETRICAL(23)
    YMINB=QUADRIC_GEOMETRICAL(24)

    RYTOP=QUADRIC_GEOMETRICAL(29)
    RYBOT=QUADRIC_GEOMETRICAL(30)
    RXTRL=QUADRIC_GEOMETRICAL(31)
    RXLED=QUADRIC_GEOMETRICAL(32)

    CEN_INT(1)=QUADRIC_KINETIC(4)
    CEN_INT(2)=QUADRIC_KINETIC(5)


    !!------确定法向量------!
    !NXT=COX
    !NYT=COY
    !NXB=-COX
    !NYB=-COY
    !IF( DABS(NXT**2.0D0+NYT**2.0D0-1.0D0)>CRITERIA  )THEN
    !    WRITE(*,*) NSTEP,"法向量非单位向量，坐标转换有问题"
    !    STOP
    !END IF

    !-------------------------------------开始求解交点-------------------------------------!

    !------X方向网格线Y/V------!
    DO J=1,JM,1
        COEQ2= COX2
        COEQ1= COXY*Y(J) +COX
        COEQ0= COY2*Y(J)**2.0D0 + COY*Y(J) + COM
        DETA2= COEQ1**2.0D0 - 4.0D0*COEQ2*COEQ0
        IF( DETA2>=CRITERIA )THEN
            SOLUTION1=( -COEQ1 + DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            SOLUTION2=( -COEQ1 - DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            !XREL=T11*( SOLUTION-CEN_INT(1) )+T12*( Y(J)-CEN_INT(2) )
            !YREL=T21*( SOLUTION-CEN_INT(1) )+T22*( Y(J)-CEN_INT(2) )

            !!判断相对坐标系中在平板范围内
            !IF( XREL<RXTRL+CRITERIA .AND. XREL>RXLED-CRITERIA )THEN
            !    IF( DABS(YREL-RYTOP)>CRITERIA )THEN
            !        WRITE(*,*) "交点不在平面上，ERRORJ:",J
            !    END IF

            IF( DABS( SOLUTION1 - DBLE(FLOOR(SOLUTION1/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION1-LEIN) / DX3 - 0.5D0)
                TYPEVX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VX(I,J)=XPV(I)
                IB_IPSVL_VX(I,J)=VBTEMP
                IB_IPSVL_VXU(I,J)=UBTEMP

            ELSE IF( DABS( SOLUTION1 - DBLE(FLOOR(SOLUTION1/DX3))*DX3 - 0.5D0*DX3 )>=CRITERIA )THEN!非网格交点
                I=IL+CEILING( (SOLUTION1-LEIN) / DX3 - 0.5D0)
                TYPEVX(I,J)=+1
                CALL VELOCITY_LB(QUADRIC_KINETIC,SOLUTION1,Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VX(I,J)=SOLUTION1
                IB_IPSVL_VX(I,J)=VBTEMP
                IB_IPSVL_VXU(I,J)=UBTEMP

            END IF

            IF( DABS( SOLUTION2 - DBLE(FLOOR(SOLUTION2/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION2-LEIN) / DX3 - 0.5D0)
                TYPEVX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VX(I,J)=XPV(I)
                IB_IPSVL_VX(I,J)=VBTEMP
                IB_IPSVL_VXU(I,J)=UBTEMP

            ELSE IF( DABS( SOLUTION2 - DBLE(FLOOR(SOLUTION2/DX3))*DX3 - 0.5D0*DX3 )>=CRITERIA )THEN!非网格交点                
                I=IL+FLOOR( (SOLUTION2-LEIN) / DX3 - 0.5D0)
                TYPEVX(I,J)=-1
                CALL VELOCITY_LB(QUADRIC_KINETIC,SOLUTION2,Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VX(I,J)=SOLUTION2
                IB_IPSVL_VX(I,J)=VBTEMP
                IB_IPSVL_VXU(I,J)=UBTEMP

            END IF

            !------开始求解内点------!
            DO N=I,IM,1
                IF(TYPEVX(N,J)==-10 .AND. TYPEVY(N,J)==-10)THEN
                    CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(N),Y(J),UBTEMP,VBTEMP)
                    IB_ITSCT_VX(N,J)=XPV(N)
                    IB_IPSVL_VX(N,J)=VBTEMP
                    IB_IPSVL_VXU(N,J)=UBTEMP
                ELSE IF(TYPEVX(N,J)==10 .AND. TYPEVY(N,J)==10)THEN
                    EXIT
                END IF
            END DO

        ELSE IF( DABS(DETA2)<CRITERIA )THEN
            SOLUTION=( -COEQ1 )/( 2.0D0*COEQ2 )
            IF( DABS( SOLUTION - DBLE(FLOOR(SOLUTION/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION-LEIN) / DX3 - 0.5D0)
                TYPEVX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VX(I,J)=XPV(I)
                IB_IPSVL_VX(I,J)=VBTEMP
                IB_IPSVL_VXU(I,J)=UBTEMP
            END IF

        END IF

    END DO

    !------Y方向网格线XPV/V------!
    DO I=1,IM,1
        COEQ2= COY2
        COEQ1= COXY*XPV(I) +COY
        COEQ0= COX2*XPV(I)**2.0D0 + COX*XPV(I) + COM
        DETA2= COEQ1**2.0D0 - 4.0D0*COEQ2*COEQ0
        IF( DETA2>=CRITERIA )THEN!存在两个交点
            SOLUTION1=( -COEQ1 + DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            SOLUTION2=( -COEQ1 - DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            !XREL=T11*( XPV(I)-CEN_INT(1) )+T12*( SOLUTION-CEN_INT(2) )
            !YREL=T21*( XPV(I)-CEN_INT(1) )+T22*( SOLUTION-CEN_INT(2) )
            !
            !!本例为相对坐标系下在平板范围内
            !IF( XREL<RXTRL+CRITERIA .AND. XREL>RXLED-CRITERIA )THEN
            !    IF( DABS(YREL-RYTOP)>CRITERIA )THEN
            !        WRITE(*,*) "交点不在平面上，ERRORI:",I
            !    END IF

            IF( DABS( SOLUTION1 - DBLE(IDNINT(SOLUTION1/DX3))*DX3 )<CRITERIA )THEN!网格交点 DABS(COX)<CRITERIA
                J=JB+IDNINT( (SOLUTION1-BOIN) / DX3 )
                TYPEVY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VY(I,J)=Y(J)
                IB_IPSVL_VY(I,J)=VBTEMP
                IB_IPSVL_VYU(I,J)=UBTEMP

            ELSE IF( DABS( SOLUTION1 - DBLE(IDNINT(SOLUTION1/DX3))*DX3 )>=CRITERIA )THEN!非网格交点
                J=JB+CEILING( (SOLUTION1-BOIN) / DX3 )
                TYPEVY(I,J)=+1
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),SOLUTION1,UBTEMP,VBTEMP)
                IB_ITSCT_VY(I,J)=SOLUTION1
                IB_IPSVL_VY(I,J)=VBTEMP
                IB_IPSVL_VYU(I,J)=UBTEMP
            END IF

            IF( DABS( SOLUTION2 - DBLE(IDNINT(SOLUTION2/DX3))*DX3 )<CRITERIA )THEN!网格交点 DABS(COX)<CRITERIA
                J=JB+IDNINT( (SOLUTION2-BOIN) / DX3 )
                TYPEVY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VY(I,J)=Y(J)
                IB_IPSVL_VY(I,J)=VBTEMP
                IB_IPSVL_VYU(I,J)=UBTEMP

            ELSE IF( DABS( SOLUTION2 - DBLE(IDNINT(SOLUTION2/DX3))*DX3 )>=CRITERIA )THEN!非网格交点
                J=JB+FLOOR( (SOLUTION2-BOIN) / DX3 )
                TYPEVY(I,J)=-1
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),SOLUTION2,UBTEMP,VBTEMP)
                IB_ITSCT_VY(I,J)=SOLUTION2
                IB_IPSVL_VY(I,J)=VBTEMP
                IB_IPSVL_VYU(I,J)=UBTEMP
            END IF
            
            !------开始求解内点------!
            DO N=J,JM,1
                IF(TYPEVX(I,N)==-10 .AND. TYPEVY(I,N)==-10)THEN
                    CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(N),UBTEMP,VBTEMP)
                    IB_ITSCT_VY(I,N)=Y(N)
                    IB_IPSVL_VY(I,N)=VBTEMP
                    IB_IPSVL_VYU(I,N)=UBTEMP
                ELSE IF(TYPEVX(I,N)==10 .AND. TYPEVY(I,N)==10)THEN
                    EXIT
                END IF
            END DO

        ELSE IF( DABS(DETA2)<CRITERIA )THEN!存在1个交点
            SOLUTION=( -COEQ1 )/( 2.0D0*COEQ2 )
            IF( DABS( SOLUTION - DBLE(IDNINT(SOLUTION/DX3))*DX3 )<CRITERIA )THEN!网格交点
                J=JB+IDNINT( (SOLUTION-BOIN) / DX3 )
                TYPEVY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,XPV(I),Y(J),UBTEMP,VBTEMP)
                IB_ITSCT_VY(I,J)=Y(J)
                IB_IPSVL_VY(I,J)=VBTEMP
                IB_IPSVL_VYU(I,J)=UBTEMP
            END IF

        END IF

    END DO

    !------X方向网格线YPU/U------!
    DO J=1,JM,1
        COEQ2= COX2
        COEQ1= COXY*YPU(J) +COX
        COEQ0= COY2*YPU(J)**2.0D0 + COY*YPU(J) + COM
        DETA2= COEQ1**2.0D0 - 4.0D0*COEQ2*COEQ0
        IF( DETA2>=CRITERIA )THEN!网格线与曲面存在两个交点
            SOLUTION1=( -COEQ1 + DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            SOLUTION2=( -COEQ1 - DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            !XREL=T11*( SOLUTION-CEN_INT(1) )+T12*( YPU(J)-CEN_INT(2) )
            !YREL=T21*( SOLUTION-CEN_INT(1) )+T22*( YPU(J)-CEN_INT(2) )
            !
            !!判断相对坐标系中在平板范围内
            !IF( XREL<RXTRL+CRITERIA .AND. XREL>RXLED-CRITERIA )THEN
            !    IF( DABS(YREL-RYTOP)>CRITERIA )THEN
            !        WRITE(*,*) "交点不在平面上，ERRORJ:",J
            !    END IF

            IF( DABS( SOLUTION1 - DBLE(IDNINT(SOLUTION1/DX3))*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION1-LEIN) / DX3 )
                TYPEUX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UX(I,J)=X(I)
                IB_IPSVL_UX(I,J)=UBTEMP
                IB_IPSVL_UXV(I,J)=VBTEMP

            ELSE IF( DABS( SOLUTION1 - DBLE(IDNINT(SOLUTION1/DX3))*DX3 )>=CRITERIA )THEN!非网格交点
                I=IL+CEILING( (SOLUTION1-LEIN) / DX3 )
                TYPEUX(I,J)=+1
                CALL VELOCITY_LB(QUADRIC_KINETIC,SOLUTION1,YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UX(I,J)=SOLUTION1
                IB_IPSVL_UX(I,J)=UBTEMP
                IB_IPSVL_UXV(I,J)=VBTEMP
            END IF

            IF( DABS( SOLUTION2 - DBLE(IDNINT(SOLUTION2/DX3))*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION2-LEIN) / DX3 )
                TYPEUX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UX(I,J)=X(I)
                IB_IPSVL_UX(I,J)=UBTEMP
                IB_IPSVL_UXV(I,J)=VBTEMP

            ELSE IF( DABS( SOLUTION2 - DBLE(IDNINT(SOLUTION2/DX3))*DX3 )>=CRITERIA )THEN!非网格交点
                I=IL+FLOOR( (SOLUTION2-LEIN) / DX3 )
                TYPEUX(I,J)=-1
                CALL VELOCITY_LB(QUADRIC_KINETIC,SOLUTION2,YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UX(I,J)=SOLUTION2
                IB_IPSVL_UX(I,J)=UBTEMP
                IB_IPSVL_UXV(I,J)=VBTEMP
            END IF
            
            !------开始求解内点------!
            DO N=I,IM,1
                IF(TYPEUX(N,J)==-10 .AND. TYPEUY(N,J)==-10)THEN
                    CALL VELOCITY_LB(QUADRIC_KINETIC,X(N),YPU(J),UBTEMP,VBTEMP)
                    IB_ITSCT_UX(N,J)=X(N)
                    IB_IPSVL_UX(N,J)=UBTEMP
                    IB_IPSVL_UXV(N,J)=VBTEMP
                ELSE IF(TYPEUX(N,J)==10 .AND. TYPEUY(N,J)==10)THEN
                    EXIT
                END IF
            END DO
            
        ELSE IF( DABS(DETA2)<CRITERIA )THEN!网格线与曲面存在1个交点
            SOLUTION=( -COEQ1 )/( 2.0D0*COEQ2 )
            IF( DABS( SOLUTION - DBLE(IDNINT(SOLUTION/DX3))*DX3 )<CRITERIA )THEN!网格交点
                I=IL+IDNINT( (SOLUTION-LEIN) / DX3 )
                TYPEUX(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UX(I,J)=X(I)
                IB_IPSVL_UX(I,J)=UBTEMP
                IB_IPSVL_UXV(I,J)=VBTEMP
            END IF
        END IF
    END DO

    !------Y方向网格线X/U------!
    DO I=1,IM,1
        COEQ2= COY2
        COEQ1= COXY*X(I) +COY
        COEQ0= COX2*X(I)**2.0D0 + COX*X(I) + COM
        DETA2= COEQ1**2.0D0 - 4.0D0*COEQ2*COEQ0
        IF( DETA2>=CRITERIA )THEN!存在两个交点
            SOLUTION1=( -COEQ1 + DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            SOLUTION2=( -COEQ1 - DSQRT(DETA2) )/( 2.0D0*COEQ2 )
            !XREL=T11*( X(I)-CEN_INT(1) )+T12*( SOLUTION-CEN_INT(2) )
            !YREL=T21*( X(I)-CEN_INT(1) )+T22*( SOLUTION-CEN_INT(2) )
            !
            !!本例为相对坐标系下在平板范围内
            !IF( XREL<RXTRL+CRITERIA .AND. XREL>RXLED-CRITERIA )THEN
            !    IF( DABS(YREL-RYTOP)>CRITERIA )THEN
            !        WRITE(*,*) "交点不在平面上，ERRORI:",I
            !    END IF

            IF( DABS( SOLUTION1 - DBLE(FLOOR(SOLUTION1/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点 DABS(COX)<CRITERIA
                J=JB+IDNINT( (SOLUTION1-BOIN) / DX3 - 0.5D0 )
                TYPEUY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UY(I,J)=YPU(J)
                IB_IPSVL_UY(I,J)=UBTEMP
                IB_IPSVL_UYV(I,J)=VBTEMP

            ELSE IF( DABS( SOLUTION1 - DBLE(FLOOR(SOLUTION1/DX3))*DX3 - 0.5D0*DX3 )>=CRITERIA )THEN!非网格交点
                J=JB+CEILING( (SOLUTION1-BOIN) / DX3 - 0.5D0 )
                TYPEUY(I,J)=+1
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),SOLUTION1,UBTEMP,VBTEMP)
                IB_ITSCT_UY(I,J)=SOLUTION1
                IB_IPSVL_UY(I,J)=UBTEMP
                IB_IPSVL_UYV(I,J)=VBTEMP

            END IF


            IF( DABS( SOLUTION2 - DBLE(FLOOR(SOLUTION2/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点 DABS(COX)<CRITERIA
                J=JB+IDNINT( (SOLUTION2-BOIN) / DX3 - 0.5D0 )
                TYPEUY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UY(I,J)=YPU(J)
                IB_IPSVL_UY(I,J)=UBTEMP
                IB_IPSVL_UYV(I,J)=VBTEMP

            ELSE IF( DABS( SOLUTION2 - DBLE(FLOOR(SOLUTION2/DX3))*DX3 - 0.5D0*DX3 )>=CRITERIA )THEN!非网格交点
                J=JB+FLOOR( (SOLUTION2-BOIN) / DX3 - 0.5D0 )
                TYPEUY(I,J)=-1
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),SOLUTION2,UBTEMP,VBTEMP)
                IB_ITSCT_UY(I,J)=SOLUTION2
                IB_IPSVL_UY(I,J)=UBTEMP
                IB_IPSVL_UYV(I,J)=VBTEMP

            END IF
            
            !------开始求解内点------!
            DO N=J,JM,1
                IF(TYPEUX(I,N)==-10 .AND. TYPEUY(I,N)==-10)THEN
                    CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(N),UBTEMP,VBTEMP)
                    IB_ITSCT_UY(I,N)=YPU(N)
                    IB_IPSVL_UY(I,N)=UBTEMP
                    IB_IPSVL_UYV(I,N)=VBTEMP
                ELSE IF(TYPEUX(I,N)==10 .AND. TYPEUY(I,N)==10)THEN
                    EXIT
                END IF
            END DO
            
        ELSE IF( DABS(DETA2)<CRITERIA )THEN!存在1个交点
            SOLUTION=( -COEQ1 )/( 2.0D0*COEQ2 )
            IF( DABS( SOLUTION - DBLE(FLOOR(SOLUTION/DX3))*DX3 - 0.5D0*DX3 )<CRITERIA )THEN!网格交点 DABS(COX)<CRITERIA
                J=JB+IDNINT( (SOLUTION-BOIN) / DX3 - 0.5D0 )
                TYPEUY(I,J)=0
                CALL VELOCITY_LB(QUADRIC_KINETIC,X(I),YPU(J),UBTEMP,VBTEMP)
                IB_ITSCT_UY(I,J)=YPU(J)
                IB_IPSVL_UY(I,J)=UBTEMP
                IB_IPSVL_UYV(I,J)=VBTEMP
            END IF

        END IF

    END DO



    RETURN
    END SUBROUTINE