    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !*************************************���Ŀǰʱ���߽罻����һʱ����ٶ�******************************************!
    SUBROUTINE IBM_CAL_VELON_FOR_BN1
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    REAL(KIND=8)::VALUE_X,VALUE_Y
    INTEGER::BOUNDARY_ID

    REAL(KIND=8)::UBTEMP,VBTEMP

    !------U------!
    DO J=0,JM,1
        DO I=1,IM,1
            !------X------!
            IF     (TYPEUX(I,J)==0)THEN
                IBN1_IPSVL_UXN(I,J)=UN(I,J)

            ELSE IF(TYPEUX(I,J)==1)THEN
                IF     (TYPEUX(I-1,J)==-10 .AND. TYPEUXM(I,J)==1 .AND. TYPEUXM(I-1,J)==-10)THEN
                    IF     (IB_ITSCT_UXN(I,J)<IB_ITSCT_UX(I,J)-CRITERIA)THEN!+1+1+1+1+1
                        CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),IB_IPSVL_UXN(I,J),X(I),IB_ITSCT_UX(I,J),IB_ITSCT_UXN(I,J))
                    ELSE IF(IB_ITSCT_UXN(I,J)>IB_ITSCT_UX(I,J)+CRITERIA)THEN!+3+3+3+3+3
                        VALUE_X=IB_ITSCT_UX(I,J)
                        VALUE_Y=YPU(J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_UXN(I,J)-IB_ITSCT_UX(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_UXN(I,J)=IB_IPSVL_UXN(I,J)!0000000
                    END IF
                ELSE IF(TYPEUX(I-1,J)==-10 .AND. TYPEUXM(I-1,J)==1 .AND. TYPEUXM(I-2,J)==-10)THEN!+2+2+2+2+2
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I-1,J),X(I),IB_ITSCT_UX(I,J),X(I-1))
                ELSE IF(TYPEUX(I-1,J)==-10 .AND. TYPEUXM(I+1,J)==1 .AND. TYPEUXM(I  ,J)==-10)THEN!+4+4+4+4+4
                    VALUE_X=IB_ITSCT_UX(I,J)
                    VALUE_Y=YPU(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEUXM(I-1,J)==10 .AND. TYPEUXM(I,J)==10 .AND. TYPEUXM(I+1,J)==10)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I-1,J),X(I),IB_ITSCT_UX(I,J),X(I-1))
                    !��10���������һ��I/I-1�������
                ELSE IF(TYPEUXM(I,J)==0 .OR. TYPEUXM(I-1,J)==0)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I-1,J),X(I),IB_ITSCT_UX(I,J),X(I-1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEUXM(I-1,J)==-10 .AND. TYPEUXM(I,J)==-10 .AND. TYPEUXM(I+1,J)==-10)THEN
                    VALUE_X=IB_ITSCT_UX(I,J)
                    VALUE_Y=YPU(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  UX+1",I,J,TYPEUX(I-1,J),TYPEUX(I,J),TYPEUX(I+1,J),TYPEUXM(I-2,J),TYPEUXM(I-1,J),TYPEUXM(I,J),TYPEUXM(I+1,J),TYPEUXM(I+2,J)
                END IF

            ELSE IF(TYPEUX(I,J)==-1)THEN!(������һ�£�����ֻ���������)
                IF     (TYPEUX(I+1,J)==-10 .AND. TYPEUXM(I,J)==-1 .AND. TYPEUXM(I+1,J)==-10)THEN
                    IF     (IB_ITSCT_UX(I,J)<IB_ITSCT_UXN(I,J)-CRITERIA)THEN!-1-1-1-1-1
                        CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),IB_IPSVL_UXN(I,J),X(I),IB_ITSCT_UX(I,J),IB_ITSCT_UXN(I,J))
                    ELSE IF(IB_ITSCT_UX(I,J)>IB_ITSCT_UXN(I,J)+CRITERIA)THEN!-3-3-3-3-3
                        VALUE_X=IB_ITSCT_UX(I,J)
                        VALUE_Y=YPU(J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_UXN(I,J)-IB_ITSCT_UX(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_UXN(I,J)=IB_IPSVL_UXN(I,J)!0000000
                    END IF
                ELSE IF(TYPEUX(I+1,J)==-10 .AND. TYPEUXM(I+1,J)==-1 .AND. TYPEUXM(I+2,J)==-10)THEN!-2-2-2-2-2
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I+1,J),X(I),IB_ITSCT_UX(I,J),X(I+1))
                ELSE IF(TYPEUX(I+1,J)==-10 .AND. TYPEUXM(I-1,J)==-1 .AND. TYPEUXM(I  ,J)==-10)THEN!-4-4-4-4-4
                    VALUE_X=IB_ITSCT_UX(I,J)
                    VALUE_Y=YPU(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEUXM(I-1,J)==10 .AND. TYPEUXM(I,J)==10 .AND. TYPEUXM(I+1,J)==10)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I+1,J),X(I),IB_ITSCT_UX(I,J),X(I+1))
                    !��10���������һ��I/I+1�������
                ELSE IF(TYPEUXM(I,J)==0 .OR. TYPEUXM(I+1,J)==0)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UXN(I,J),UN(I+1,J),X(I),IB_ITSCT_UX(I,J),X(I+1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEUXM(I-1,J)==-10 .AND. TYPEUXM(I,J)==-10 .AND. TYPEUXM(I+1,J)==-10)THEN
                    VALUE_X=IB_ITSCT_UX(I,J)
                    VALUE_Y=YPU(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UXN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UXN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  UX-1",I,J,TYPEUX(I-1,J),TYPEUX(I,J),TYPEUX(I+1,J),TYPEUXM(I-2,J),TYPEUXM(I-1,J),TYPEUXM(I,J),TYPEUXM(I+1,J),TYPEUXM(I+2,J)
                END IF

            END IF

            !------Y------!
            IF     (TYPEUY(I,J)==0)THEN
                IBN1_IPSVL_UYN(I,J)=UN(I,J)

            ELSE IF(TYPEUY(I,J)==1)THEN
                IF     (TYPEUY(I,J-1)==-10 .AND. TYPEUYM(I,J)==1 .AND. TYPEUYM(I,J-1)==-10)THEN
                    IF     (IB_ITSCT_UYN(I,J)<IB_ITSCT_UY(I,J)-CRITERIA)THEN!+1+1+1+1+1
                        CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),IB_IPSVL_UYN(I,J),YPU(J),IB_ITSCT_UY(I,J),IB_ITSCT_UYN(I,J))
                    ELSE IF(IB_ITSCT_UYN(I,J)>IB_ITSCT_UY(I,J)+CRITERIA)THEN!+3+3+3+3+3
                        VALUE_X=X(I)
                        VALUE_Y=IB_ITSCT_UY(I,J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_UYN(I,J)-IB_ITSCT_UY(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_UYN(I,J)=IB_IPSVL_UYN(I,J)!0000000
                    END IF
                ELSE IF(TYPEUY(I,J-1)==-10 .AND. TYPEUYM(I,J-1)==1 .AND. TYPEUYM(I,J-2)==-10)THEN!+2+2+2+2+2
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J-1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                ELSE IF(TYPEUY(I,J-1)==-10 .AND. TYPEUYM(I,J+1)==1 .AND. TYPEUYM(I,J  )==-10)THEN!+4+4+4+4+4
                    VALUE_X=X(I)
                    VALUE_Y=IB_ITSCT_UY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEUYM(I,J-1)==10 .AND. TYPEUYM(I,J)==10 .AND. TYPEUYM(I,J+1)==10)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J-1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                    !��10���������һ��I/I-1�������
                ELSE IF(TYPEUYM(I,J)==0 .OR. TYPEUYM(I,J-1)==0)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J-1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEUYM(I,J-1)==-10 .AND. TYPEUYM(I,J)==-10 .AND. TYPEUYM(I,J+1)==-10)THEN
                    VALUE_X=X(I)
                    VALUE_Y=IB_ITSCT_UY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  UY+1",I,J,TYPEUY(I,J-1),TYPEUY(I,J),TYPEUY(I,J+1),TYPEUYM(I,J-2),TYPEUYM(I,J-1),TYPEUYM(I,J),TYPEUYM(I,J+1),TYPEUYM(I,J+2)
                END IF

            ELSE IF(TYPEUY(I,J)==-1)THEN!(������һ�£�����ֻ���������)
                IF     (TYPEUY(I,J+1)==-10 .AND. TYPEUYM(I,J)==-1 .AND. TYPEUYM(I,J+1)==-10)THEN
                    IF     (IB_ITSCT_UY(I,J)<IB_ITSCT_UYN(I,J)-CRITERIA)THEN!-1-1-1-1-1
                        CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),IB_IPSVL_UYN(I,J),YPU(J),IB_ITSCT_UY(I,J),IB_ITSCT_UYN(I,J))
                    ELSE IF(IB_ITSCT_UY(I,J)>IB_ITSCT_UYN(I,J)+CRITERIA)THEN!-3-3-3-3-3
                        VALUE_X=X(I)
                        VALUE_Y=IB_ITSCT_UY(I,J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_UYN(I,J)-IB_ITSCT_UY(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_UYN(I,J)=IB_IPSVL_UYN(I,J)!0000000
                    END IF
                ELSE IF(TYPEUY(I,J+1)==-10 .AND. TYPEUYM(I,J+1)==-1 .AND. TYPEUYM(I,J+2)==-10)THEN!-2-2-2-2-2
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J+1))
                ELSE IF(TYPEUY(I,J+1)==-10 .AND. TYPEUYM(I,J-1)==-1 .AND. TYPEUYM(I,J  )==-10)THEN!-4-4-4-4-4
                    VALUE_X=X(I)
                    VALUE_Y=IB_ITSCT_UY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEUYM(I,J-1)==10 .AND. TYPEUYM(I,J)==10 .AND. TYPEUYM(I,J+1)==10)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J+1))
                    !��10���������һ��I/I+1�������
                ELSE IF(TYPEUYM(I,J)==0 .OR. TYPEUYM(I,J+1)==0)THEN
                    CALL LINEAR_INTERPOLATION(UN(I,J),IBN1_IPSVL_UYN(I,J),UN(I,J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J+1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEUYM(I,J-1)==-10 .AND. TYPEUYM(I,J)==-10 .AND. TYPEUYM(I,J+1)==-10)THEN
                    VALUE_X=X(I)
                    VALUE_Y=IB_ITSCT_UY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,IBN1_IPSVL_UYN(I,J),VBTEMP)
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_UYN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  UY-1",I,J,TYPEUY(I,J-1),TYPEUY(I,J),TYPEUY(I,J+1),TYPEUYM(I,J-2),TYPEUYM(I,J-1),TYPEUYM(I,J),TYPEUYM(I,J+1),TYPEUYM(I,J+2)
                END IF

            END IF

        END DO
    END DO

    !------V------!
    DO J=1,JM,1
        DO I=0,IM,1
            !------X------!
            IF     (TYPEVX(I,J)==0)THEN
                IBN1_IPSVL_VXN(I,J)=VN(I,J)

            ELSE IF(TYPEVX(I,J)==1)THEN
                IF     (TYPEVX(I-1,J)==-10 .AND. TYPEVXM(I,J)==1 .AND. TYPEVXM(I-1,J)==-10)THEN
                    IF     (IB_ITSCT_VXN(I,J)<IB_ITSCT_VX(I,J)-CRITERIA)THEN!+1+1+1+1+1
                        CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),IB_IPSVL_VXN(I,J),XPV(I),IB_ITSCT_VX(I,J),IB_ITSCT_VXN(I,J))
                    ELSE IF(IB_ITSCT_VXN(I,J)>IB_ITSCT_VX(I,J)+CRITERIA)THEN!+3+3+3+3+3
                        VALUE_X=IB_ITSCT_VX(I,J)
                        VALUE_Y=Y(J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_VXN(I,J)-IB_ITSCT_VX(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_VXN(I,J)=IB_IPSVL_VXN(I,J)!0000000
                    END IF
                ELSE IF(TYPEVX(I-1,J)==-10 .AND. TYPEVXM(I-1,J)==1 .AND. TYPEVXM(I-2,J)==-10)THEN!+2+2+2+2+2
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I-1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                ELSE IF(TYPEVX(I-1,J)==-10 .AND. TYPEVXM(I+1,J)==1 .AND. TYPEVXM(I  ,J)==-10)THEN!+4+4+4+4+4
                    VALUE_X=IB_ITSCT_VX(I,J)
                    VALUE_Y=Y(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEVXM(I-1,J)==10 .AND. TYPEVXM(I,J)==10 .AND. TYPEVXM(I+1,J)==10)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I-1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                    !��10���������һ��I/I-1�������
                ELSE IF(TYPEVXM(I,J)==0 .OR. TYPEVXM(I-1,J)==0)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I-1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEVXM(I-1,J)==-10 .AND. TYPEVXM(I,J)==-10 .AND. TYPEVXM(I+1,J)==-10)THEN
                    VALUE_X=IB_ITSCT_VX(I,J)
                    VALUE_Y=Y(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  VX+1",I,J,TYPEVX(I-1,J),TYPEVX(I,J),TYPEVX(I+1,J),TYPEVXM(I-2,J),TYPEVXM(I-1,J),TYPEVXM(I,J),TYPEVXM(I+1,J),TYPEVXM(I+2,J)
                END IF

            ELSE IF(TYPEVX(I,J)==-1)THEN!(������һ�£�����ֻ���������)
                IF     (TYPEVX(I+1,J)==-10 .AND. TYPEVXM(I,J)==-1 .AND. TYPEVXM(I+1,J)==-10)THEN
                    IF     (IB_ITSCT_VX(I,J)<IB_ITSCT_VXN(I,J)-CRITERIA)THEN!-1-1-1-1-1
                        CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),IB_IPSVL_VXN(I,J),XPV(I),IB_ITSCT_VX(I,J),IB_ITSCT_VXN(I,J))
                    ELSE IF(IB_ITSCT_VX(I,J)>IB_ITSCT_VXN(I,J)+CRITERIA)THEN!-3-3-3-3-3
                        VALUE_X=IB_ITSCT_VX(I,J)
                        VALUE_Y=Y(J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_VXN(I,J)-IB_ITSCT_VX(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_VXN(I,J)=IB_IPSVL_VXN(I,J)!0000000
                    END IF
                ELSE IF(TYPEVX(I+1,J)==-10 .AND. TYPEVXM(I+1,J)==-1 .AND. TYPEVXM(I+2,J)==-10)THEN!-2-2-2-2-2
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I+1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I+1))
                ELSE IF(TYPEVX(I+1,J)==-10 .AND. TYPEVXM(I-1,J)==-1 .AND. TYPEVXM(I  ,J)==-10)THEN!-4-4-4-4-4
                    VALUE_X=IB_ITSCT_VX(I,J)
                    VALUE_Y=Y(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEVXM(I-1,J)==10 .AND. TYPEVXM(I,J)==10 .AND. TYPEVXM(I+1,J)==10)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I+1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I+1))
                    !��10���������һ��I/I+1�������
                ELSE IF(TYPEVXM(I,J)==0 .OR. TYPEVXM(I+1,J)==0)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VXN(I,J),VN(I+1,J),XPV(I),IB_ITSCT_VX(I,J),XPV(I+1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEVXM(I-1,J)==-10 .AND. TYPEVXM(I,J)==-10 .AND. TYPEVXM(I+1,J)==-10)THEN
                    VALUE_X=IB_ITSCT_VX(I,J)
                    VALUE_Y=Y(J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VXN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VXN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  VX-1",I,J,TYPEVX(I-1,J),TYPEVX(I,J),TYPEVX(I+1,J),TYPEVXM(I-2,J),TYPEVXM(I-1,J),TYPEVXM(I,J),TYPEVXM(I+1,J),TYPEVXM(I+2,J)
                END IF

            END IF

            !------Y------!
            IF     (TYPEVY(I,J)==0)THEN
                IBN1_IPSVL_VYN(I,J)=VN(I,J)

            ELSE IF(TYPEVY(I,J)==1)THEN
                IF     (TYPEVY(I,J-1)==-10 .AND. TYPEVYM(I,J)==1 .AND. TYPEVYM(I,J-1)==-10)THEN
                    IF     (IB_ITSCT_VYN(I,J)<IB_ITSCT_VY(I,J)-CRITERIA)THEN!+1+1+1+1+1
                        CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),IB_IPSVL_VYN(I,J),Y(J),IB_ITSCT_VY(I,J),IB_ITSCT_VYN(I,J))
                    ELSE IF(IB_ITSCT_VYN(I,J)>IB_ITSCT_VY(I,J)+CRITERIA)THEN!+3+3+3+3+3
                        VALUE_X=XPV(I)
                        VALUE_Y=IB_ITSCT_VY(I,J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_VYN(I,J)-IB_ITSCT_VY(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_VYN(I,J)=IB_IPSVL_VYN(I,J)!0000000
                    END IF
                ELSE IF(TYPEVY(I,J-1)==-10 .AND. TYPEVYM(I,J-1)==1 .AND. TYPEVYM(I,J-2)==-10)THEN!+2+2+2+2+2
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J-1),Y(J),IB_ITSCT_VY(I,J),Y(J-1))
                ELSE IF(TYPEVY(I,J-1)==-10 .AND. TYPEVYM(I,J+1)==1 .AND. TYPEVYM(I,J  )==-10)THEN!+4+4+4+4+4
                    VALUE_X=XPV(I)
                    VALUE_Y=IB_ITSCT_VY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEVYM(I,J-1)==10 .AND. TYPEVYM(I,J)==10 .AND. TYPEVYM(I,J+1)==10)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J-1),Y(J),IB_ITSCT_VY(I,J),Y(J-1))
                    !��10���������һ��I/I-1�������
                ELSE IF(TYPEVYM(I,J)==0 .OR. TYPEVYM(I,J-1)==0)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J-1),Y(J),IB_ITSCT_VY(I,J),Y(J-1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEVYM(I,J-1)==-10 .AND. TYPEVYM(I,J)==-10 .AND. TYPEVYM(I,J+1)==-10)THEN
                    VALUE_X=XPV(I)
                    VALUE_Y=IB_ITSCT_VY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  VY+1",I,J,TYPEVY(I,J-1),TYPEVY(I,J),TYPEVY(I,J+1),TYPEVYM(I,J-2),TYPEVYM(I,J-1),TYPEVYM(I,J),TYPEVYM(I,J+1),TYPEVYM(I,J+2)
                END IF

            ELSE IF(TYPEVY(I,J)==-1)THEN!(������һ�£�����ֻ���������)
                IF     (TYPEVY(I,J+1)==-10 .AND. TYPEVYM(I,J)==-1 .AND. TYPEVYM(I,J+1)==-10)THEN
                    IF     (IB_ITSCT_VY(I,J)<IB_ITSCT_VYN(I,J)-CRITERIA)THEN!-1-1-1-1-1
                        CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),IB_IPSVL_VYN(I,J),Y(J),IB_ITSCT_VY(I,J),IB_ITSCT_VYN(I,J))
                    ELSE IF(IB_ITSCT_VY(I,J)>IB_ITSCT_VYN(I,J)+CRITERIA)THEN!-3-3-3-3-3
                        VALUE_X=XPV(I)
                        VALUE_Y=IB_ITSCT_VY(I,J)
                        CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                        IF     ( BOUNDARY_ID==1 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                        ELSE IF( BOUNDARY_ID==2 )THEN
                            CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                        ELSE IF( BOUNDARY_ID==0 )THEN
                            WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                        END IF
                    ELSE IF(DABS(IB_ITSCT_VYN(I,J)-IB_ITSCT_VY(I,J))<=CRITERIA)THEN
                        IBN1_IPSVL_VYN(I,J)=IB_IPSVL_VYN(I,J)!0000000
                    END IF
                ELSE IF(TYPEVY(I,J+1)==-10 .AND. TYPEVYM(I,J+1)==-1 .AND. TYPEVYM(I,J+2)==-10)THEN!-2-2-2-2-2
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J+1),Y(J),IB_ITSCT_VY(I,J),Y(J+1))
                ELSE IF(TYPEVY(I,J+1)==-10 .AND. TYPEVYM(I,J-1)==-1 .AND. TYPEVYM(I,J  )==-10)THEN!-4-4-4-4-4
                    VALUE_X=XPV(I)
                    VALUE_Y=IB_ITSCT_VY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                    END IF
                    !��9���������һ��I-1,I,I+1���������ڵ�
                ELSE IF(TYPEVYM(I,J-1)==10 .AND. TYPEVYM(I,J)==10 .AND. TYPEVYM(I,J+1)==10)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J+1),Y(J),IB_ITSCT_VY(I,J),Y(J+1))
                    !��10���������һ��I/I+1�������
                ELSE IF(TYPEVYM(I,J)==0 .OR. TYPEVYM(I,J+1)==0)THEN
                    CALL LINEAR_INTERPOLATION(VN(I,J),IBN1_IPSVL_VYN(I,J),VN(I,J+1),Y(J),IB_ITSCT_VY(I,J),Y(J+1))
                    !��11���������һ��I-1,I,I+1���ǹ������ڵ�
                ELSE IF(TYPEVYM(I,J-1)==-10 .AND. TYPEVYM(I,J)==-10 .AND. TYPEVYM(I,J+1)==-10)THEN
                    VALUE_X=XPV(I)
                    VALUE_Y=IB_ITSCT_VY(I,J)
                    CALL DETERMINE_BOUNDARY_ID(VALUE_X,VALUE_Y,BOUNDARY_ID)
                    IF     ( BOUNDARY_ID==1 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N1,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==2 )THEN
                        CALL VELOCITY_LB(QUADRIC_KINETIC_N2,VALUE_X,VALUE_Y,UBTEMP,IBN1_IPSVL_VYN(I,J))
                    ELSE IF( BOUNDARY_ID==0 )THEN
                        WRITE(992,*)"ERROR:IBN1_IPSVL_VYN",I,J
                    END IF
                ELSE
                    WRITE(992,*)NSTEP,"  VY-1",I,J,TYPEVY(I,J-1),TYPEVY(I,J),TYPEVY(I,J+1),TYPEVYM(I,J-2),TYPEVYM(I,J-1),TYPEVYM(I,J),TYPEVYM(I,J+1),TYPEVYM(I,J+2)
                END IF

            END IF

        END DO
    END DO

    RETURN
    END SUBROUTINE