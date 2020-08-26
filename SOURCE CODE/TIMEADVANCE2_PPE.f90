    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !**************************************************求解压力泊松方程*****************************************************!
    SUBROUTINE TIMEADVANCE2_PPE
    USE DECLARATION
    USE CSV_FILE
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::RHS(:,:),ERRORPHI(:,:),TEMP(:,:)
    REAL(KIND=8)::RELAX=1.9D0
    REAL(KIND=8),ALLOCATABLE::A1(:),A2(:),A3(:),B1(:),B2(:),B3(:)!C1,C2,C3
    REAL(KIND=8),ALLOCATABLE::LAMBDAX(:),LAMBDAY(:),DX_PMAT(:,:),DX_PINV(:,:),DY_QMAT(:,:),DY_QINV(:,:)
    INTEGER,PARAMETER::ITERATION_MAX=300000
    REAL(KIND=8)::ERRORPHI_HISTORY(ITERATION_MAX)
    REAL(KIND=8)::RHS_MAX,CORRECTION
    INTEGER::SOLVERTYPE=2 !2-矩阵分解，3-PSOR
    INTEGER::PSORTYPE=4

    SAVE LAMBDAX,LAMBDAY,DX_PMAT,DX_PINV,DY_QMAT,DY_QINV

    !------差分离散系数------!
    IF( ALLOCATED(A1) )THEN
        DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    END IF
    ALLOCATE( A1(IM-1),A2(IM-1),A3(IM-1),B1(JM-1),B2(JM-1),B3(JM-1) )
    A1=0.0D0
    A2=0.0D0
    A3=0.0D0
    B1=0.0D0
    B2=0.0D0
    B3=0.0D0
    DO I=1,IM-1,1
        A1(I)=2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
        A2(I)=2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
        A3(I)=2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    END DO
    DO J=1,JM-1,1
        B1(J)=2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
        B2(J)=2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
        B3(J)=2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    END DO


    ALLOCATE( RHS(IM-1,JM-1),ERRORPHI(0:IM,0:JM) )
    RHS=0.0D0
    PHI=0.0D0
    ERRORPHI=0.0D0
    !------计算右端项------!
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            !!A-B C-N
            !RHS(I,J)=( ( UHAT(I+1,J)-UHAT(I,J) )/( X(I+1)-X(I) )+( VHAT(I,J+1)-VHAT(I,J) )/( Y(J+1)-Y(J) ) )/DT
            !R-K C-N
            RHS(I,J)=( ( UHAT(I+1,J)-UHAT(I,J) )/( X(I+1)-X(I) )+( VHAT(I,J+1)-VHAT(I,J) )/( Y(J+1)-Y(J) ) )/DT/ALPHA(NSUBSTEP)
        END DO
    END DO
    RHS_MAX=MAXVAL(DABS(RHS))

    IF(SOLVERTYPE==2)THEN
        ALLOCATE( TEMP(1:IM-1,1:JM-1) )
        TEMP=0.0D0
        !第一步就准备系数矩阵
        IF(NSTEP==NSTART .AND. NSUBSTEP==1)THEN
            ALLOCATE( LAMBDAX(IM-1),LAMBDAY(JM-1),DX_PMAT(IM-1,IM-1),DX_PINV(IM-1,IM-1),DY_QMAT(JM-1,JM-1),DY_QINV(JM-1,JM-1) )

            LAMBDAX=0.0D0
            LAMBDAY=0.0D0
            DX_PMAT=0.0D0
            DX_PINV=0.0D0
            DY_QMAT=0.0D0
            DY_QINV=0.0D0

            !构建DX矩阵
            DO I=1,IM-2
                DX_PMAT(I+1,I)=A1(I+1)
                DX_PMAT(I,I+1)=A3(I)
            END DO
            DO I=2,IM-2
                DX_PMAT(I,I)=-A2(I)
            END DO
            DX_PMAT(1,1)=A1(1)*BCPHI_AL-A2(1)
            DX_PMAT(IM-1,IM-1)=A3(IM-1)*BCPHI_AR-A2(IM-1)

            !输出构建结果
            !OPEN(UNIT=10,FILE="MATRIX_4_EIGENDECOMPOSITION_X.CSV")
            !CALL CSV_WRITE(10, TRANSPOSE(DX_PMAT))
            !CLOSE(10)
            OPEN(UNIT=10,FILE="MATRIX_4_EIGENDECOMPOSITION.DAT")
            WRITE(10,*) DX_PMAT!((DX_PMAT(I,J),J=1,IM-1),I=1,IM-1)
            CLOSE(10)
            OPEN(UNIT=10,FILE="DIM_MATRIX.DAT")
            WRITE(10,*) SHAPE(DX_PMAT)
            CLOSE(10)

            !特征分解
            CALL SPECTRAL_DECOMPOSITION

            !读取分解结果

            DX_PMAT=0.0D0
            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTOR.DAT")
            READ(10,*) DX_PMAT
            CLOSE(10)

            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTORINV.DAT")
            READ(10,*) DX_PINV
            CLOSE(10)

            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVALUE.DAT")
            READ(10,*) LAMBDAX
            CLOSE(10)

            !OPEN(UNIT=10,FILE="EIGENVECTOR_X.CSV")
            !CALL CSV_WRITE(10, TRANSPOSE(DX_PMAT))
            !CLOSE(10)
            !
            !OPEN(UNIT=10,FILE="EIGENVALUE_X.CSV")
            !CALL CSV_WRITE(10, LAMBDAX )
            !CLOSE(10)


            !构建DY矩阵
            DO J=1,JM-2
                DY_QMAT(J+1,J)=B1(J+1)
                DY_QMAT(J,J+1)=B3(J)
            END DO
            DO J=2,JM-2
                DY_QMAT(J,J)=-B2(J)
            END DO
            DY_QMAT(1,1)=B1(1)*BCPHI_AB-B2(1)
            DY_QMAT(JM-1,JM-1)=B3(JM-1)*BCPHI_AT-B2(JM-1)

            !DY要转置一下，别给忘了
            DY_QMAT=TRANSPOSE(DY_QMAT)

            !输出构建结果
            !OPEN(UNIT=10,FILE="MATRIX_4_EIGENDECOMPOSITION_Y.CSV")
            !CALL CSV_WRITE(10, TRANSPOSE(DY_QMAT))
            !CLOSE(10)
            OPEN(UNIT=10,FILE="MATRIX_4_EIGENDECOMPOSITION.DAT")
            WRITE(10,*)  DY_QMAT!((DY_QMAT(I,J),J=1,JM-1),I=1,JM-1)
            CLOSE(10)
            OPEN(UNIT=10,FILE="DIM_MATRIX.DAT")
            WRITE(10,*) SHAPE(DY_QMAT)
            CLOSE(10)

            !特征分解
            CALL SPECTRAL_DECOMPOSITION

            !读取分解结果

            DY_QMAT=0.0D0
            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTOR.DAT")
            READ(10,*) DY_QMAT
            CLOSE(10)

            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVECTORINV.DAT")
            READ(10,*) DY_QINV
            CLOSE(10)

            OPEN(UNIT=10,FILE="DECOMPOSITION_EIGENVALUE.DAT")
            READ(10,*) LAMBDAY
            CLOSE(10)

        END IF

        !计算RHS,补充边界条件
        I=1
        DO J=1,JM-1,1
            RHS(I,J)=RHS(I,J)-BCPHI_BL*A1(I)
        END DO

        I=IM-1
        DO J=1,JM-1,1
            RHS(I,J)=RHS(I,J)-BCPHI_BR*A3(I)
        END DO

        J=1
        DO I=1,IM-1,1
            RHS(I,J)=RHS(I,J)-BCPHI_BB*B1(J)
        END DO

        J=JM-1
        DO I=1,IM-1,1
            RHS(I,J)=RHS(I,J)-BCPHI_BT*B3(J)
        END DO

        !计算TEMP
        RHS=MATMUL(DX_PINV,RHS)
        RHS=MATMUL(RHS,DY_QMAT)

        DO J=1,JM-1,1
            DO I=1,IM-1,1
                TEMP(I,J)=RHS(I,J)/( LAMBDAX(I)+LAMBDAY(J) )
            END DO
        END DO

        TEMP=MATMUL(DX_PMAT,TEMP)
        TEMP=MATMUL(TEMP,DY_QINV)

        PHI(1:IM-1,1:JM-1)=TEMP

        !------边界条件------!
        PHI(0       ,1:JM-1:1)=BCPHI_AL*PHI(1       ,1:JM-1:1)+BCPHI_BL
        PHI(IM      ,1:JM-1:1)=BCPHI_AR*PHI(IM-1    ,1:JM-1:1)+BCPHI_BR
        PHI(:       ,0       )=BCPHI_AB*PHI(:       ,1       )+BCPHI_BB
        PHI(:       ,JM      )=BCPHI_AT*PHI(:       ,JM-1    )+BCPHI_BT

        !------数值锚------!
        IF( DABS(BCPHI_AL*BCPHI_AR*BCPHI_AB*BCPHI_AT-1.0D0)<=CRITERIA )THEN
            CORRECTION=PHI(0,0)
            PHI=PHI-CORRECTION
            WRITE(1234567,*) "PHI修正：",NSTEP,CORRECTION
        END IF



    ELSE IF(SOLVERTYPE==3)THEN
        ALLOCATE( TEMP(0:IM,0:JM) )
        TEMP=0.0D0


        !------PSOR------!
        DO N=1,ITERATION_MAX,1
            IF(PSORTYPE==1)THEN
                !PSOR
                DO J=1,JM-1,1!JM-1,1,-1!
                    IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                        DO I=1,IM-1,1!IM-1,1,-1!
                            IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                            ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                            ELSE!中间区
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                            END IF
                        END DO
                    ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                        DO I=1,IM-1,1!IM-1,1,-1!
                            IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                            ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                            ELSE!中间区
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                            END IF
                        END DO
                    ELSE!中间区
                        DO I=1,IM-1,1!IM-1,1,-1!
                            IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                            ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                            ELSE!中间区
                                TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                            END IF
                        END DO
                    END IF
                END DO
            ELSE IF(PSORTYPE==2)THEN
                !两次来回扫PSOR
                IF( MOD( N,2 )==0 )THEN
                    DO J=1,JM-1,1!JM-1,1,-1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                ELSE IF( MOD( N,2 )==1 )THEN
                    DO J=JM-1,1,-1!1,JM-1,1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                END IF
            ELSE IF(PSORTYPE==4)THEN
                !四次来回扫PSOR
                IF( MOD( N,4 )==0 )THEN
                    DO J=1,JM-1,1!JM-1,1,-1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                ELSE IF( MOD( N,4 )==1 )THEN
                    DO J=JM-1,1,-1!1,JM-1,1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=1,IM-1,1!IM-1,1,-1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                ELSE IF( MOD( N,4 )==2 )THEN
                    DO J=1,JM-1,1!JM-1,1,-1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                ELSE IF( MOD( N,4 )==3 )THEN
                    DO J=JM-1,1,-1!1,JM-1,1!
                        IF(J==1)THEN!下边界TEMP(I,J-1)=BCPHI_AB*TEMP(I,J)+BCPHI_BB
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*BCPHI_BB   + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AB*B1(J))
                                END IF
                            END DO
                        ELSE IF(J==JM-1)THEN!上边界TEMP(I,J+1)=BCPHI_AT*TEMP(I,J)+BCPHI_BT
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*BCPHI_BT   - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AT*B3(J))
                                END IF
                            END DO
                        ELSE!中间区
                            DO I=IM-1,1,-1!1,IM-1,1!
                                IF(I==1)THEN        !左边界TEMP(I-1,J)=BCPHI_AL*TEMP(I,J)+BCPHI_BL
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*BCPHI_BL    + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AL*A1(I))
                                ELSE IF(I==IM-1)THEN!右边界TEMP(I+1,J)=BCPHI_AR*TEMP(I,J)+BCPHI_BR
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*BCPHI_BR    + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J)-BCPHI_AR*A3(I))
                                ELSE!中间区
                                    TEMP(I,J)=(1.0D0-RELAX)*TEMP(I,J)+RELAX*( A1(I)*TEMP(I-1,J) + A3(I)*TEMP(I+1,J) + B1(J)*TEMP(I,J-1) + B3(J)*TEMP(I,J+1) - RHS(I,J) )/(A2(I)+B2(J))
                                END IF
                            END DO
                        END IF
                    END DO
                END IF
            END IF

            !------边界条件------!
            TEMP(0       ,1:JM-1:1)=BCPHI_AL*TEMP(1       ,1:JM-1:1)+BCPHI_BL
            TEMP(IM      ,1:JM-1:1)=BCPHI_AR*TEMP(IM-1    ,1:JM-1:1)+BCPHI_BR
            TEMP(:       ,0       )=BCPHI_AB*TEMP(:       ,1       )+BCPHI_BB
            TEMP(:       ,JM      )=BCPHI_AT*TEMP(:       ,JM-1    )+BCPHI_BT
            !------数值锚------!
            IF( DABS(BCPHI_AL*BCPHI_AR*BCPHI_AB*BCPHI_AT-1.0D0)<=CRITERIA )THEN
                CORRECTION=TEMP(1,1)
                TEMP=TEMP-CORRECTION
                WRITE(1234567,*) NSTEP,CORRECTION
            END IF
            !------判断收敛------!
            ERRORPHI=DABS(TEMP-PHI)
            PHI=TEMP
            ERRORPHI_HISTORY(N)=MAXVAL(ERRORPHI)

            IF( ERRORPHI_HISTORY(N)<=CRITERIA )THEN
                IF(N>400)THEN
                    IF( ERRORPHI_HISTORY(N)<=100.0D0*CRITERIA*MAXVAL(DABS(PHI)) .OR. ( ERRORPHI_HISTORY(N-200)/ERRORPHI_HISTORY(N)<=1.05D0 .AND. ERRORPHI_HISTORY(N-400)/ERRORPHI_HISTORY(N-200)<=1.05D0 ) )THEN
                        WRITE(*,*) N,2
                        EXIT
                    END IF
                ELSE IF( ERRORPHI_HISTORY(N)<=100.0D0*CRITERIA*MAXVAL(DABS(PHI)) )THEN
                    WRITE(*,*) N,1
                    EXIT
                END IF

            END IF

            !第二跳出条件
            !IF(N>400)THEN
            !    IF( ERRORPHI_HISTORY(N-200)/ERRORPHI_HISTORY(N)<=1.05D0 .AND. ERRORPHI_HISTORY(N-400)/ERRORPHI_HISTORY(N-200)<=1.05D0 )THEN
            !        WRITE(*,*) N,2
            !        EXIT
            !    END IF
            !END IF

            IF(MOD(N,200)==1) WRITE(*,*) N,ERRORPHI_HISTORY(N)

        END DO

    END IF

    RETURN
    END SUBROUTINE