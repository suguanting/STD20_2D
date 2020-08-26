    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************��RK-CN��ɢ��N-S�������̽���ʱ���ƽ�********************************************!
    SUBROUTINE TIMEADVANCE1_ME_RK_CN
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::RU(:,:),RV(:,:)
    REAL(KIND=8),ALLOCATABLE::A1(:,:),A2(:,:),A3(:,:),B1(:,:),B2(:,:),B3(:,:)!C1,C2,C3
    REAL(KIND=8),ALLOCATABLE::DUH(:,:),DVH(:,:)
    REAL(KIND=8),ALLOCATABLE::DUP(:,:),DVP(:,:)
    REAL(KIND=8)::UR1,UL1,UU1,UD1,VR1,VL1,VU1,VD1!������ϳ��ٶȣ�k-1ʱ���
    REAL(KIND=8)::UR2,UL2,UU2,UD2,VR2,VL2,VU2,VD2!������ϳ��ٶȣ�k-2ʱ���
    REAL(KIND=8)::RPK1,RCXK1,RCYK1,RCXK2,RCYK2,RVXK1,RVYK1!�Ҷ��������
    REAL(KIND=8),ALLOCATABLE::TRIA(:),TRIB(:),TRIC(:),TRID(:)
    REAL(KIND=8),ALLOCATABLE::TRIC1(:),TRID1(:)

    REAL(KIND=8)::D1,D2
    !------���ڱ������------!
    REAL(KIND=8),ALLOCATABLE::BA_L(:),BA_R(:),BA_T(:),BA_B(:)
    REAL(KIND=8)::VELO_CORRECTION
    INTEGER::OUTLET_R,OUTLET_L,OUTLET_T,OUTLET_B


    ALLOCATE( RU(IM,0:JM),RV(0:IM,JM),DUP(IM,0:JM),DVP(0:IM,JM),DUH(IM,0:JM),DVH(0:IM,JM) )
    ALLOCATE( BA_L(JM-1),BA_R(JM-1),BA_T(IM-1),BA_B(IM-1) )

    RU=0.0D0
    RV=0.0D0
    DUH=0.0D0
    DUP=0.0D0
    DVH=0.0D0
    DVP=0.0D0

    !-----------------------------------------------V----------------------------------------------!

    !------�����ɢϵ������ֱ�������������㣬������ճ����ʱ���-2.0D0------!
    IF( ALLOCATED(A1) )THEN
        DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    END IF
    ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    A1=0.0D0
    A2=0.0D0
    A3=0.0D0
    B1=0.0D0
    B2=0.0D0
    B3=0.0D0
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            A1(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
            A2(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
            A3(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
            B1(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
            B2(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
            B3(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
        END DO
    END DO

    !------RHS------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            IF( TYPEVX(I,J)/=-10 .AND. TYPEVY(I,J)/=-10 )THEN
                VU1=0.5D0*( VK1(I,J+1)+VK1(I,J) )
                VD1=0.5D0*( VK1(I,J-1)+VK1(I,J) )
                CALL LINEAR_INTERPOLATION( VK1(I+1,J),VR1,VK1(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
                CALL LINEAR_INTERPOLATION( VK1(I  ,J),VL1,VK1(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( UK1(I+1,J),UR1,UK1(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( UK1(I  ,J),UL1,UK1(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                VU2=0.5D0*( VK2(I,J+1)+VK2(I,J) )
                VD2=0.5D0*( VK2(I,J-1)+VK2(I,J) )
                CALL LINEAR_INTERPOLATION( VK2(I+1,J),VR2,VK2(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
                CALL LINEAR_INTERPOLATION( VK2(I  ,J),VL2,VK2(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( UK2(I+1,J),UR2,UK2(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( UK2(I  ,J),UL2,UK2(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                D1=YPU(J-1)-Y(J)
                D2=YPU(J)-Y(J)
                !------P/Y------!
                RPK1=( P(I,J)-P(I,J-1) )/( YPU(J)-YPU(J-1) )
                !------UV/X------!
                IF( .NOT. ISNAN(CONVECT_VXK1(I,J)) )THEN
                    RCXK1=CONVECT_VXK1(I,J)
                ELSE
                    RCXK1=( UR1*VR1-UL1*VL1 )/( X(I+1)-X(I) )
                END IF
                IF( .NOT. ISNAN(CONVECT_VXK2(I,J)) )THEN
                    RCXK2=CONVECT_VXK2(I,J)
                ELSE
                    RCXK2=( UR2*VR2-UL2*VL2 )/( X(I+1)-X(I) )
                END IF
                !------VV/Y------!
                IF( .NOT. ISNAN(CONVECT_VYK1(I,J)) )THEN
                    RCYK1=CONVECT_VYK1(I,J)
                ELSE
                    RCYK1=( D2*D2*VD1*VD1+(D1*D1-D2*D2)*VK1(I,J)*VK1(I,J)-D1*D1*VU1*VU1 )/( D1*D2*(D2-D1) )
                END IF
                IF( .NOT. ISNAN(CONVECT_VYK2(I,J)) )THEN
                    RCYK2=CONVECT_VYK2(I,J)
                ELSE
                    RCYK2=( D2*D2*VD2*VD2+(D1*D1-D2*D2)*VK2(I,J)*VK2(I,J)-D1*D1*VU2*VU2 )/( D1*D2*(D2-D1) )
                END IF
                !------V/XX------!
                IF( ISNAN(LAPLACE_VXK1(I,J)) )THEN
                    RVXK1=-2.0D0*A1(I,J)*VK1(I-1,J) + 2.0D0*A2(I,J)*VK1(I,J) - 2.0D0*A3(I,J)*VK1(I+1,J)
                ELSE
                    RVXK1=ALPHA(NSUBSTEP)*DT/Re*LAPLACE_VXK1(I,J)
                END IF
                !------V/YY------!
                IF( ISNAN(LAPLACE_VYK1(I,J)) )THEN
                    RVYK1=-2.0D0*B1(I,J)*VK1(I,J-1) + 2.0D0*B2(I,J)*VK1(I,J) - 2.0D0*B3(I,J)*VK1(I,J+1)
                ELSE
                    RVYK1=ALPHA(NSUBSTEP)*DT/Re*LAPLACE_VYK1(I,J)
                END IF

                RV(I,J)=DT*(-ALPHA(NSUBSTEP)*RPK1-GAMA(NSUBSTEP)*(RCXK1+RCYK1)-RHO(NSUBSTEP)*(RCXK2+RCYK2))+(RVXK1+RVYK1)
            END IF
        END DO
    END DO

    !------ADI------!
    !------XSWEEP------!
    IF( ALLOCATED(TRIA) )THEN
        DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    END IF
    ALLOCATE( TRIA(IM-1),TRIB(IM-1),TRIC(IM-1),TRID(IM-1) )
    ALLOCATE( TRIC1(IM-1),TRID1(IM-1) )
    TRIA=0.0D0
    TRIB=0.0D0
    TRIC=0.0D0
    TRID=0.0D0
    TRIC1=0.0D0
    TRID1=0.0D0

    DO J=2,JM-1,1

        IF(J==2)THEN!�±߽�
            DO I=2,IM-2,1
                TRIA(I)=A1(I,J)
                TRIB(I)=1.0D0-A2(I,J)
                TRIC(I)=A3(I,J)
                TRID(I)=RV(I,J)
            END DO

            I=1
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RV(I,J)-A1(I,J)*BCV_CL*(VK1(1 ,J)-VK1(0   ,J))/(XPV(1 )-XPV(0   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RV(I,J)-A3(I,J)*BCV_CR*(VK1(IM,J)-VK1(IM-1,J))/(XPV(IM)-XPV(IM-1))

        ELSE IF(J==JM-1)THEN!�ϱ߽�
            DO I=2,IM-2,1
                TRIA(I)=A1(I,J)
                TRIB(I)=1.0D0-A2(I,J)
                TRIC(I)=A3(I,J)
                TRID(I)=RV(I,J)
            END DO

            I=1
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RV(I,J)-A1(I,J)*BCV_CL*(VK1(1 ,J)-VK1(0   ,J))/(XPV(1 )-XPV(0   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RV(I,J)-A3(I,J)*BCV_CR*(VK1(IM,J)-VK1(IM-1,J))/(XPV(IM)-XPV(IM-1))

        ELSE!�ڲ�
            DO I=2,IM-2,1
                IF(TYPEVX(I,J)==10)THEN
                    TRIA(I)=A1(I,J)
                    TRIB(I)=1.0D0-A2(I,J)
                    TRIC(I)=A3(I,J)
                    TRID(I)=RV(I,J)
                ELSE IF(TYPEVX(I,J)==0)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0
                    TRIC(I)=0.0D0
                    TRID(I)=IB_RVX(I,J)
                ELSE IF(TYPEVX(I,J)==1)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0+IB_CVX(I,J)
                    TRIC(I)=IB_AVX(I,J)
                    TRID(I)=RV(I,J)+IB_RVX(I,J)
                ELSE IF(TYPEVX(I,J)==-1)THEN
                    TRIA(I)=IB_AVX(I,J)
                    TRIB(I)=1.0D0+IB_CVX(I,J)
                    TRIC(I)=0.0D0
                    TRID(I)=RV(I,J)+IB_RVX(I,J)
                ELSE IF(TYPEVX(I,J)==-10)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0
                    TRIC(I)=0.0D0
                    TRID(I)=IB_RVX(I,J)
                END IF
            END DO

            I=1
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RV(I,J)-A1(I,J)*BCV_CL*(VK1(1 ,J)-VK1(0   ,J))/(XPV(1 )-XPV(0   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RV(I,J)-A3(I,J)*BCV_CR*(VK1(IM,J)-VK1(IM-1,J))/(XPV(IM)-XPV(IM-1))

        END IF

        !------FORWARD SWEEP------!
        TRIC1(1)=TRIC(1)/TRIB(1)
        DO I=2,IM-2,1
            TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
        END DO
        TRID1(1)=TRID(1)/TRIB(1)
        DO I=2,IM-1,1
            TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
        END DO
        !------BACK SUBSTITUTION------!
        DVP(IM-1,J)=TRID1(IM-1)
        DO I=IM-2,1,-1
            DVP(I,J)=TRID1(I)-TRIC1(I)*DVP(I+1,J)
        END DO

    END DO
    !------SWEEP����------!
    !------ʩ�ӱ߽�����------!
    DVP(0       ,2:JM-1:1)=BCV_AL*DVP(1       ,2:JM-1:1)+BCV_CL*(VK1(2 ,2:JM-1:1)-VK1(1   ,2:JM-1:1))/(XPV(1 )-XPV(0   ))
    DVP(IM      ,2:JM-1:1)=BCV_AR*DVP(IM-1    ,2:JM-1:1)+BCV_CR*(VK1(IM,2:JM-1:1)-VK1(IM-1,2:JM-1:1))/(XPV(IM)-XPV(IM-1))
    DVP(:       ,1       )=BCV_AB*DVP(:       ,2       )+BCV_CB*(VK1(:,2 )-VK1(:,1   ))/(Y(2 )-Y(1   ))
    DVP(:       ,JM      )=BCV_AT*DVP(:       ,JM-1    )+BCV_CT*(VK1(:,JM)-VK1(:,JM-1))/(Y(JM)-Y(JM-1))
    !------YSWEEP------!
    IF( ALLOCATED(TRIA) )THEN
        DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    END IF
    ALLOCATE( TRIA(2:JM-1),TRIB(2:JM-1),TRIC(2:JM-1),TRID(2:JM-1) )
    ALLOCATE( TRIC1(2:JM-1),TRID1(2:JM-1) )
    TRIA=0.0D0
    TRIB=0.0D0
    TRIC=0.0D0
    TRID=0.0D0
    TRIC1=0.0D0
    TRID1=0.0D0

    DO I=1,IM-1,1
        IF(I==1)THEN!��߽�
            DO J=3,JM-2,1
                TRIA(J)=B1(I,J)
                TRIB(J)=1.0D0-B2(I,J)
                TRIC(J)=B3(I,J)
                TRID(J)=DVP(I,J)
            END DO

            J=2
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DVP(I,J)-B1(I,J)*BCV_CB*(VK1(I,2 )-VK1(I,1   ))/(Y(2 )-Y(1   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DVP(I,J)-B3(I,J)*BCV_CT*(VK1(I,JM)-VK1(I,JM-1))/(Y(JM)-Y(JM-1))

        ELSE IF(I==IM-1)THEN!�ұ߽�
            DO J=3,JM-2,1
                TRIA(J)=B1(I,J)
                TRIB(J)=1.0D0-B2(I,J)
                TRIC(J)=B3(I,J)
                TRID(J)=DVP(I,J)
            END DO

            J=2
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DVP(I,J)-B1(I,J)*BCV_CB*(VK1(I,2 )-VK1(I,1   ))/(Y(2 )-Y(1   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DVP(I,J)-B3(I,J)*BCV_CT*(VK1(I,JM)-VK1(I,JM-1))/(Y(JM)-Y(JM-1))

        ELSE!�ڲ�
            DO J=3,JM-2,1
                IF(TYPEVY(I,J)==10)THEN
                    TRIA(J)=B1(I,J)
                    TRIB(J)=1.0D0-B2(I,J)
                    TRIC(J)=B3(I,J)
                    TRID(J)=DVP(I,J)
                ELSE IF(TYPEVY(I,J)==0)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0
                    TRIC(J)=0.0D0
                    TRID(J)=IB_RVY(I,J)
                ELSE IF(TYPEVY(I,J)==1)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0+IB_CVY(I,J)
                    TRIC(J)=IB_BVY(I,J)
                    TRID(J)=DVP(I,J)+IB_RVY(I,J)
                ELSE IF(TYPEVY(I,J)==-1)THEN
                    TRIA(J)=IB_BVY(I,J)
                    TRIB(J)=1.0D0+IB_CVY(I,J)
                    TRIC(J)=0.0D0
                    TRID(J)=DVP(I,J)+IB_RVY(I,J)
                ELSE IF(TYPEVY(I,J)==-10)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0
                    TRIC(J)=0.0D0
                    TRID(J)=IB_RVY(I,J)
                END IF
            END DO

            J=2
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DVP(I,J)-B1(I,J)*BCV_CB*(VK1(I,2 )-VK1(I,1   ))/(Y(2 )-Y(1   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DVP(I,J)-B3(I,J)*BCV_CT*(VK1(I,JM)-VK1(I,JM-1))/(Y(JM)-Y(JM-1))

        END IF
        !------FORWARD SWEEP------!
        TRIC1(2)=TRIC(2)/TRIB(2)
        DO J=3,JM-2,1
            TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
        END DO
        TRID1(2)=TRID(2)/TRIB(2)
        DO J=3,JM-1,1
            TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
        END DO
        !------BACK SUBSTITUTION------!
        DVH(I,JM-1)=TRID1(JM-1)
        DO J=JM-2,2,-1
            DVH(I,J)=TRID1(J)-TRIC1(J)*DVH(I,J+1)
        END DO

    END DO
    !------SWEEP����------!
    !------ʩ�ӱ߽�����------!
    DVH(0       ,2:JM-1:1)=BCV_AL*DVH(1       ,2:JM-1:1)+BCV_CL*(VK1(2 ,2:JM-1:1)-VK1(1   ,2:JM-1:1))/(XPV(1 )-XPV(0   ))
    DVH(IM      ,2:JM-1:1)=BCV_AR*DVH(IM-1    ,2:JM-1:1)+BCV_CR*(VK1(IM,2:JM-1:1)-VK1(IM-1,2:JM-1:1))/(XPV(IM)-XPV(IM-1))
    DVH(:       ,1       )=BCV_AB*DVH(:       ,2       )+BCV_CB*(VK1(:,2 )-VK1(:,1   ))/(Y(2 )-Y(1   ))
    DVH(:       ,JM      )=BCV_AT*DVH(:       ,JM-1    )+BCV_CT*(VK1(:,JM)-VK1(:,JM-1))/(Y(JM)-Y(JM-1))

    !-----------------------------------------------U----------------------------------------------!

    !------�����ɢϵ��------!
    IF( ALLOCATED(A1) )THEN
        DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    END IF
    ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    A1=0.0D0
    A2=0.0D0
    A3=0.0D0
    B1=0.0D0
    B2=0.0D0
    B3=0.0D0
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            A1(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
            A2(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
            A3(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
            B1(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
            B2(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
            B3(I,J)=-ALPHA(NSUBSTEP)*DT/(2.0D0*Re)*2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
        END DO
    END DO

    !------RHS------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            IF( TYPEUX(I,J)/=-10 .AND. TYPEUY(I,J)/=-10 )THEN
                UR1=0.5D0*( UK1(I+1,J)+UK1(I,J) )
                UL1=0.5D0*( UK1(I-1,J)+UK1(I,J) )
                CALL LINEAR_INTERPOLATION( UK1(I,J+1),UU1,UK1(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
                CALL LINEAR_INTERPOLATION( UK1(I,J  ),UD1,UK1(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( VK1(I,J+1),VU1,VK1(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( VK1(I,J  ),VD1,VK1(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
                UR2=0.5D0*( UK2(I+1,J)+UK2(I,J) )
                UL2=0.5D0*( UK2(I-1,J)+UK2(I,J) )
                CALL LINEAR_INTERPOLATION( UK2(I,J+1),UU2,UK2(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
                CALL LINEAR_INTERPOLATION( UK2(I,J  ),UD2,UK2(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( VK2(I,J+1),VU2,VK2(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( VK2(I,J  ),VD2,VK2(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
                D1=XPV(I-1)-X(I)
                D2=XPV(I)-X(I)
                !------P/X------!
                RPK1=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
                !------UU/X------!
                IF( .NOT. ISNAN(CONVECT_UXK1(I,J)) )THEN
                    RCXK1=CONVECT_UXK1(I,J)
                ELSE
                    RCXK1=( D2*D2*UL1*UL1+(D1*D1-D2*D2)*UK1(I,J)*UK1(I,J)-D1*D1*UR1*UR1 )/( D1*D2*(D2-D1) )
                END IF
                IF( .NOT. ISNAN(CONVECT_UXK2(I,J)) )THEN
                    RCXK2=CONVECT_UXK2(I,J)
                ELSE
                    RCXK2=( D2*D2*UL2*UL2+(D1*D1-D2*D2)*UK2(I,J)*UK2(I,J)-D1*D1*UR2*UR2 )/( D1*D2*(D2-D1) )
                END IF
                !------UV/Y------!
                IF( .NOT. ISNAN(CONVECT_UYK1(I,J)) )THEN
                    RCYK1=CONVECT_UYK1(I,J)
                ELSE
                    RCYK1=( UU1*VU1-UD1*VD1 )/( Y(J+1)-Y(J) )
                END IF
                IF( .NOT. ISNAN(CONVECT_UYK2(I,J)) )THEN
                    RCYK2=CONVECT_UYK2(I,J)
                ELSE
                    RCYK2=( UU2*VU2-UD2*VD2 )/( Y(J+1)-Y(J) )
                END IF
                !------U/XX------!
                IF( ISNAN(LAPLACE_UXK1(I,J)) )THEN
                    RVXK1=-2.0D0*A1(I,J)*UK1(I-1,J) + 2.0D0*A2(I,J)*UK1(I,J) - 2.0D0*A3(I,J)*UK1(I+1,J)
                ELSE
                    RVXK1=ALPHA(NSUBSTEP)*DT/Re*LAPLACE_UXK1(I,J)
                END IF

                !------U/YY------!
                IF( ISNAN(LAPLACE_UYK1(I,J)) )THEN
                    RVYK1=-2.0D0*B1(I,J)*UK1(I,J-1) + 2.0D0*B2(I,J)*UK1(I,J) - 2.0D0*B3(I,J)*UK1(I,J+1)
                ELSE
                    RVYK1=ALPHA(NSUBSTEP)*DT/Re*LAPLACE_UYK1(I,J)
                END IF

                RU(I,J)=DT*(-ALPHA(NSUBSTEP)*RPK1-GAMA(NSUBSTEP)*(RCXK1+RCYK1)-RHO(NSUBSTEP)*(RCXK2+RCYK2))+(RVXK1+RVYK1)
            END IF
        END DO
    END DO

    !------ADI------!
    !------XSWEEP------!
    IF( ALLOCATED(TRIA) )THEN
        DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    END IF
    ALLOCATE( TRIA(2:IM-1),TRIB(2:IM-1),TRIC(2:IM-1),TRID(2:IM-1) )
    ALLOCATE( TRIC1(2:IM-1),TRID1(2:IM-1) )
    TRIA=0.0D0
    TRIB=0.0D0
    TRIC=0.0D0
    TRID=0.0D0
    TRIC1=0.0D0
    TRID1=0.0D0

    DO J=1,JM-1,1
        IF(J==1)THEN!�±߽�
            DO I=3,IM-2,1
                TRIA(I)=A1(I,J)
                TRIB(I)=1.0D0-A2(I,J)
                TRIC(I)=A3(I,J)
                TRID(I)=RU(I,J)
            END DO

            I=2
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RU(I,J)-A1(I,J)*BCU_CL*(UK1(2 ,J)-UK1(1   ,J))/(X(2 )-X(1   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RU(I,J)-A3(I,J)*BCU_CR*(UK1(IM,J)-UK1(IM-1,J))/(X(IM)-X(IM-1))

        ELSE IF(J==JM-1)THEN!�ϱ߽�
            DO I=3,IM-2,1
                TRIA(I)=A1(I,J)
                TRIB(I)=1.0D0-A2(I,J)
                TRIC(I)=A3(I,J)
                TRID(I)=RU(I,J)
            END DO

            I=2
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RU(I,J)-A1(I,J)*BCU_CL*(UK1(2 ,J)-UK1(1   ,J))/(X(2 )-X(1   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RU(I,J)-A3(I,J)*BCU_CR*(UK1(IM,J)-UK1(IM-1,J))/(X(IM)-X(IM-1))

        ELSE!�ڲ�
            DO I=3,IM-2,1
                IF(TYPEUX(I,J)==10)THEN
                    TRIA(I)=A1(I,J)
                    TRIB(I)=1.0D0-A2(I,J)
                    TRIC(I)=A3(I,J)
                    TRID(I)=RU(I,J)
                ELSE IF(TYPEUX(I,J)==0)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0
                    TRIC(I)=0.0D0
                    TRID(I)=IB_RUX(I,J)
                ELSE IF(TYPEUX(I,J)==1)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0+IB_CUX(I,J)
                    TRIC(I)=IB_AUX(I,J)
                    TRID(I)=RU(I,J)+IB_RUX(I,J)
                ELSE IF(TYPEUX(I,J)==-1)THEN
                    TRIA(I)=IB_AUX(I,J)
                    TRIB(I)=1.0D0+IB_CUX(I,J)
                    TRIC(I)=0.0D0
                    TRID(I)=RU(I,J)+IB_RUX(I,J)
                ELSE IF(TYPEUX(I,J)==-10)THEN
                    TRIA(I)=0.0D0
                    TRIB(I)=1.0D0
                    TRIC(I)=0.0D0
                    TRID(I)=IB_RUX(I,J)
                END IF
            END DO

            I=2
            TRIA(I)=0.0D0
            TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
            TRIC(I)=A3(I,J)
            TRID(I)=RU(I,J)-A1(I,J)*BCU_CL*(UK1(2 ,J)-UK1(1   ,J))/(X(2 )-X(1   ))

            I=IM-1
            TRIA(I)=A1(I,J)
            TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
            TRIC(I)=0.0D0
            TRID(I)=RU(I,J)-A3(I,J)*BCU_CR*(UK1(IM,J)-UK1(IM-1,J))/(X(IM)-X(IM-1))

        END IF

        !------FORWARD SWEEP------!
        TRIC1(2)=TRIC(2)/TRIB(2)
        DO I=3,IM-2,1
            TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
        END DO
        TRID1(2)=TRID(2)/TRIB(2)
        DO I=3,IM-1,1
            TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
        END DO
        !------BACK SUBSTITUTION------!
        DUP(IM-1,J)=TRID1(IM-1)
        DO I=IM-2,2,-1
            DUP(I,J)=TRID1(I)-TRIC1(I)*DUP(I+1,J)
        END DO

    END DO
    !------SWEEP����------!
    !------ʩ�ӱ߽�����------!���������ƣ����¸����̶�ֵ��
    DUP(1       ,1:JM-1:1)=BCU_AL*DUP(2       ,1:JM-1:1)+BCU_CL*(UK1(2 ,1:JM-1:1)-UK1(1   ,1:JM-1:1))/(X(2 )-X(1   ))
    DUP(IM      ,1:JM-1:1)=BCU_AR*DUP(IM-1    ,1:JM-1:1)+BCU_CR*(UK1(IM,1:JM-1:1)-UK1(IM-1,1:JM-1:1))/(X(IM)-X(IM-1))
    DUP(:       ,0       )=BCU_AB*DUP(:       ,1       )+BCU_CB*(UK1(:,1 )-UK1(:,0   ))/(YPU(1 )-YPU(0   ))
    DUP(:       ,JM      )=BCU_AT*DUP(:       ,JM-1    )+BCU_CT*(UK1(:,JM)-UK1(:,JM-1))/(YPU(JM)-YPU(JM-1))

    !------YSWEEP------!
    IF( ALLOCATED(TRIA) )THEN
        DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    END IF
    ALLOCATE( TRIA(JM-1),TRIB(JM-1),TRIC(JM-1),TRID(JM-1) )
    ALLOCATE( TRIC1(JM-1),TRID1(JM-1) )
    TRIA=0.0D0
    TRIB=0.0D0
    TRIC=0.0D0
    TRID=0.0D0
    TRIC1=0.0D0
    TRID1=0.0D0

    DO I=2,IM-1,1
        IF(I==2)THEN!��߽�
            DO J=2,JM-2,1
                TRIA(J)=B1(I,J)
                TRIB(J)=1.0D0-B2(I,J)
                TRIC(J)=B3(I,J)
                TRID(J)=DUP(I,J)
            END DO

            J=1
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DUP(I,J)-B1(I,J)*BCU_CB*(UK1(I,1 )-UK1(I,0   ))/(YPU(1 )-YPU(0   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DUP(I,J)-B3(I,J)*BCU_CT*(UK1(I,JM)-UK1(I,JM-1))/(YPU(JM)-YPU(JM-1))

        ELSE IF(I==IM-1)THEN!�ұ߽�
            DO J=2,JM-2,1
                TRIA(J)=B1(I,J)
                TRIB(J)=1.0D0-B2(I,J)
                TRIC(J)=B3(I,J)
                TRID(J)=DUP(I,J)
            END DO

            J=1
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DUP(I,J)-B1(I,J)*BCU_CB*(UK1(I,1 )-UK1(I,0   ))/(YPU(1 )-YPU(0   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DUP(I,J)-B3(I,J)*BCU_CT*(UK1(I,JM)-UK1(I,JM-1))/(YPU(JM)-YPU(JM-1))

        ELSE!�ڲ�
            DO J=2,JM-2,1
                IF(TYPEUY(I,J)==10)THEN
                    TRIA(J)=B1(I,J)
                    TRIB(J)=1.0D0-B2(I,J)
                    TRIC(J)=B3(I,J)
                    TRID(J)=DUP(I,J)
                ELSE IF(TYPEUY(I,J)==0)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0
                    TRIC(J)=0.0D0
                    TRID(J)=IB_RUY(I,J)
                ELSE IF(TYPEUY(I,J)==1)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0+IB_CUY(I,J)
                    TRIC(J)=IB_BUY(I,J)
                    TRID(J)=DUP(I,J)+IB_RUY(I,J)
                ELSE IF(TYPEUY(I,J)==-1)THEN
                    TRIA(J)=IB_BUY(I,J)
                    TRIB(J)=1.0D0+IB_CUY(I,J)
                    TRIC(J)=0.0D0
                    TRID(J)=DUP(I,J)+IB_RUY(I,J)
                ELSE IF(TYPEUY(I,J)==-10)THEN
                    TRIA(J)=0.0D0
                    TRIB(J)=1.0D0
                    TRIC(J)=0.0D0
                    TRID(J)=IB_RUY(I,J)
                END IF
            END DO

            J=1
            TRIA(J)=0.0D0
            TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
            TRIC(J)=B3(I,J)
            TRID(J)=DUP(I,J)-B1(I,J)*BCU_CB*(UK1(I,1 )-UK1(I,0   ))/(YPU(1 )-YPU(0   ))

            J=JM-1
            TRIA(J)=B1(I,J)
            TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
            TRIC(J)=0.0D0
            TRID(J)=DUP(I,J)-B3(I,J)*BCU_CT*(UK1(I,JM)-UK1(I,JM-1))/(YPU(JM)-YPU(JM-1))

        END IF
        !------FORWARD SWEEP------!
        TRIC1(1)=TRIC(1)/TRIB(1)
        DO J=2,JM-2,1
            TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
        END DO
        TRID1(1)=TRID(1)/TRIB(1)
        DO J=2,JM-1,1
            TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
        END DO
        !------BACK SUBSTITUTION------!
        DUH(I,JM-1)=TRID1(JM-1)
        DO J=JM-2,1,-1
            DUH(I,J)=TRID1(J)-TRIC1(J)*DUH(I,J+1)
        END DO

    END DO
    !------SWEEP����------!
    !------ʩ�ӱ߽�����------!
    DUH(1       ,1:JM-1:1)=BCU_AL*DUH(2       ,1:JM-1:1)+BCU_CL*(UK1(2 ,1:JM-1:1)-UK1(1   ,1:JM-1:1))/(X(2 )-X(1   ))
    DUH(IM      ,1:JM-1:1)=BCU_AR*DUH(IM-1    ,1:JM-1:1)+BCU_CR*(UK1(IM,1:JM-1:1)-UK1(IM-1,1:JM-1:1))/(X(IM)-X(IM-1))
    DUH(:       ,0       )=BCU_AB*DUH(:       ,1       )+BCU_CB*(UK1(:,1 )-UK1(:,0   ))/(YPU(1 )-YPU(0   ))
    DUH(:       ,JM      )=BCU_AT*DUH(:       ,JM-1    )+BCU_CT*(UK1(:,JM)-UK1(:,JM-1))/(YPU(JM)-YPU(JM-1))

    UHAT=UK1+DUH
    VHAT=VK1+DVH

    !!------ʩ�ӳ��ڱ߽�����(����߽�cell�������غ�)------!
    !IF(BCTYPE_L==2)THEN
    !    DO J=1,JM-1,1
    !        UHAT(1 ,J)=UHAT(2   ,J)+(X(2 )-X(1   ))/(Y(J+1)-Y(J))*(VHAT(1   ,J+1)-VHAT(1   ,J))
    !    END DO
    !END IF
    !IF(BCTYPE_R==2)THEN
    !    DO J=1,JM-1,1
    !        UHAT(IM,J)=UHAT(IM-1,J)-(X(IM)-X(IM-1))/(Y(J+1)-Y(J))*(VHAT(IM-1,J+1)-VHAT(IM-1,J))
    !    END DO
    !END IF
    !IF(BCTYPE_B==2)THEN
    !    DO I=1,IM-1,1
    !        VHAT(I,1 )=VHAT(I,2   )+(Y(2 )-Y(1   ))/(X(I+1)-X(I))*(UHAT(I+1,1   )-UHAT(I,1   ))
    !    END DO
    !END IF
    !IF(BCTYPE_T==2)THEN
    !    DO I=1,IM-1,1
    !        VHAT(I,JM)=VHAT(I,JM-1)-(Y(JM)-Y(JM-1))/(X(I+1)-X(I))*(UHAT(I+1,JM-1)-UHAT(I,JM-1))
    !    END DO
    !END IF

    !------ʩ�ӳ��ڱ߽�����(���������ٶȣ�����ȫ�������غ�)------!
    !------��߽��С------!
    DO J=1,JM-1,1
        BA_R(J)=Y(J+1)-Y(J)
        BA_L(J)=Y(J+1)-Y(J)
    END DO
    DO I=1,IM-1,1
        BA_T(I)=X(I+1)-X(I)
        BA_B(I)=X(I+1)-X(I)
    END DO
    !------�жϳ��ڱ߽�------!
    OUTLET_R=0
    OUTLET_L=0
    OUTLET_T=0
    OUTLET_B=0

    IF(BCTYPE_R==2 .OR. BCTYPE_R==5)THEN
        OUTLET_R=1
    END IF
    IF(BCTYPE_L==2 .OR. BCTYPE_L==5)THEN
        OUTLET_L=1
    END IF
    IF(BCTYPE_B==2 .OR. BCTYPE_B==5)THEN
        OUTLET_B=1
    END IF
    IF(BCTYPE_T==2 .OR. BCTYPE_T==5)THEN
        OUTLET_T=1
    END IF
    !------������ֵ------!
    IF(OUTLET_R+OUTLET_L+OUTLET_T+OUTLET_B>0)THEN
        VELO_CORRECTION=&
            & ( SUM(BA_R*UHAT(IM,1:JM-1)) &
            & - SUM(BA_L*UHAT(1 ,1:JM-1)) &
            & + SUM(BA_T*VHAT(1:IM-1,JM)) &
            & - SUM(BA_B*VHAT(1:IM-1,1 )) ) / &
            & ( SUM(BA_R)*DBLE(OUTLET_R) &
            & + SUM(BA_L)*DBLE(OUTLET_L) &
            & + SUM(BA_T)*DBLE(OUTLET_T) &
            & + SUM(BA_B)*DBLE(OUTLET_B) )
    ELSE
        VELO_CORRECTION=0.0D0
    END IF
    !------�����ٶ�------!ע��������������ķ���ͬ����������������
    UHAT(1 ,1:JM-1)=UHAT(1 ,1:JM-1)+VELO_CORRECTION*DBLE(OUTLET_L)
    UHAT(IM,1:JM-1)=UHAT(IM,1:JM-1)-VELO_CORRECTION*DBLE(OUTLET_R)
    VHAT(1:IM-1,1 )=VHAT(1:IM-1,1 )+VELO_CORRECTION*DBLE(OUTLET_B)
    VHAT(1:IM-1,JM)=VHAT(1:IM-1,JM)-VELO_CORRECTION*DBLE(OUTLET_T)

    RETURN
    END SUBROUTINE

    !!***************************************��AB-CN��ɢ��N-S�������̽���ʱ���ƽ�********************************************!
    !SUBROUTINE TIMEADVANCE1_ME_AB_CN
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !REAL(KIND=8),ALLOCATABLE::RU(:,:),RV(:,:)
    !REAL(KIND=8),ALLOCATABLE::A1(:,:),A2(:,:),A3(:,:),B1(:,:),B2(:,:),B3(:,:)!C1,C2,C3
    !REAL(KIND=8),ALLOCATABLE::DUH(:,:),DVH(:,:)
    !REAL(KIND=8),ALLOCATABLE::DUP(:,:),DVP(:,:)
    !REAL(KIND=8)::UR,UL,UU,UD,VR,VL,VU,VD!������ϳ��ٶȣ�nʱ���
    !REAL(KIND=8)::UR1,UL1,UU1,UD1,VR1,VL1,VU1,VD1!������ϳ��ٶȣ�n-1ʱ���
    !REAL(KIND=8)::RPN,RCN,RCN1,RVXN,RVYN!�Ҷ��������
    !REAL(KIND=8),ALLOCATABLE::TRIA(:),TRIB(:),TRIC(:),TRID(:)
    !REAL(KIND=8),ALLOCATABLE::TRIC1(:),TRID1(:)
    !
    !REAL(KIND=8)::D1,D2
    !
    !ALLOCATE( RU(IM,0:JM),RV(0:IM,JM),DUP(IM,0:JM),DVP(0:IM,JM),DUH(IM,0:JM),DVH(0:IM,JM) )
    !
    !RU=0.0D0
    !RV=0.0D0
    !DUH=0.0D0
    !DUP=0.0D0
    !DVH=0.0D0
    !DVP=0.0D0
    !
    !!-----------------------------------------------V----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        IF( TYPEVX(I,J)/=-10 .AND. TYPEVY(I,J)/=-10 )THEN
    !            VU=0.5D0*( VN(I,J+1)+VN(I,J) )
    !            VD=0.5D0*( VN(I,J-1)+VN(I,J) )
    !            CALL LINEAR_INTERPOLATION( VN(I+1,J),VR,VN(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
    !            CALL LINEAR_INTERPOLATION( VN(I  ,J),VL,VN(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( UN(I+1,J),UR,UN(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( UN(I  ,J),UL,UN(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            VU1=0.5D0*( VN1(I,J+1)+VN1(I,J) )
    !            VD1=0.5D0*( VN1(I,J-1)+VN1(I,J) )
    !            CALL LINEAR_INTERPOLATION( VN1(I+1,J),VR1,VN1(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
    !            CALL LINEAR_INTERPOLATION( VN1(I  ,J),VL1,VN1(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( UN1(I+1,J),UR1,UN1(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( UN1(I  ,J),UL1,UN1(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            D1=YPU(J-1)-Y(J)
    !            D2=YPU(J)-Y(J)
    !            RPN=( P(I,J)-P(I,J-1) )/( YPU(J)-YPU(J-1) )
    !            RCN =( UR *VR -UL *VL  )/( X(I+1)-X(I) )+( D2*D2*VD *VD +(D1*D1-D2*D2)*VN (I,J)*VN (I,J)-D1*D1*VU *VU  )/( D1*D2*(D2-D1) )
    !            RCN1=( UR1*VR1-UL1*VL1 )/( X(I+1)-X(I) )+( D2*D2*VD1*VD1+(D1*D1-D2*D2)*VN1(I,J)*VN1(I,J)-D1*D1*VU1*VU1 )/( D1*D2*(D2-D1) )
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*VN(I-1,J) + 2.0D0*A2(I,J)*VN(I,J) - 2.0D0*A3(I,J)*VN(I+1,J)
    !            ELSE
    !                RVXN=DT/Re*VISCOUS_VXN(I,J)
    !            END IF
    !
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*VN(I,J-1) + 2.0D0*B2(I,J)*VN(I,J) - 2.0D0*B3(I,J)*VN(I,J+1)
    !            ELSE
    !                RVYN=DT/Re*VISCOUS_VYN(I,J)
    !            END IF
    !
    !            RV(I,J)=DT*(-RPN-1.5D0*RCN+0.5D0*RCN1)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(IM-1),TRIB(IM-1),TRIC(IM-1),TRID(IM-1) )
    !ALLOCATE( TRIC1(IM-1),TRID1(IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=2,JM-1,1
    !    IF(J==2)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=2,IM-2,1
    !            IF(TYPEVX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RV(I,J)
    !            ELSE IF(TYPEVX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=IB_AVX(I,J)
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-1)THEN
    !                TRIA(I)=IB_AVX(I,J)
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            END IF
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO I=2,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO I=2,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,1,-1
    !        DVP(I,J)=TRID1(I)-TRIC1(I)*DVP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVP(0       ,2:JM-1:1)=BCV_AL*DVP(1       ,2:JM-1:1)
    !DVP(IM      ,2:JM-1:1)=BCV_AR*DVP(IM-1    ,2:JM-1:1)
    !DVP(:       ,1       )=BCV_AB*DVP(:       ,2       )
    !DVP(:       ,JM      )=BCV_AT*DVP(:       ,JM-1    )
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:JM-1),TRIB(2:JM-1),TRIC(2:JM-1),TRID(2:JM-1) )
    !ALLOCATE( TRIC1(2:JM-1),TRID1(2:JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=1,IM-1,1
    !    IF(I==1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=3,JM-2,1
    !            IF(TYPEVY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DVP(I,J)
    !            ELSE IF(TYPEVY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=IB_BVY(I,J)
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-1)THEN
    !                TRIA(J)=IB_BVY(I,J)
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            END IF
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO J=3,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO J=3,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,2,-1
    !        DVH(I,J)=TRID1(J)-TRIC1(J)*DVH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVH(0       ,2:JM-1:1)=BCV_AL*DVH(1       ,2:JM-1:1)
    !DVH(IM      ,2:JM-1:1)=BCV_AR*DVH(IM-1    ,2:JM-1:1)
    !DVH(:       ,1       )=BCV_AB*DVH(:       ,2       )
    !DVH(:       ,JM      )=BCV_AT*DVH(:       ,JM-1    )
    !
    !!-----------------------------------------------U----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        IF( TYPEUX(I,J)/=-10 .AND. TYPEUY(I,J)/=-10 )THEN
    !            UR=0.5D0*( UN(I+1,J)+UN(I,J) )
    !            UL=0.5D0*( UN(I-1,J)+UN(I,J) )
    !            CALL LINEAR_INTERPOLATION( UN(I,J+1),UU,UN(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !            CALL LINEAR_INTERPOLATION( UN(I,J  ),UD,UN(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J+1),VU,VN(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J  ),VD,VN(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !            UR1=0.5D0*( UN1(I+1,J)+UN1(I,J) )
    !            UL1=0.5D0*( UN1(I-1,J)+UN1(I,J) )
    !            CALL LINEAR_INTERPOLATION( UN1(I,J+1),UU1,UN1(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !            CALL LINEAR_INTERPOLATION( UN1(I,J  ),UD1,UN1(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( VN1(I,J+1),VU1,VN1(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( VN1(I,J  ),VD1,VN1(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !            D1=XPV(I-1)-X(I)
    !            D2=XPV(I)-X(I)
    !            RPN=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
    !            RCN =( D2*D2*UL *UL +(D1*D1-D2*D2)*UN (I,J)*UN (I,J)-D1*D1*UR *UR  )/( D1*D2*(D2-D1) )+( UU *VU -UD *VD  )/( Y(J+1)-Y(J) )
    !            RCN1=( D2*D2*UL1*UL1+(D1*D1-D2*D2)*UN1(I,J)*UN1(I,J)-D1*D1*UR1*UR1 )/( D1*D2*(D2-D1) )+( UU1*VU1-UD1*VD1 )/( Y(J+1)-Y(J) )
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*UN(I-1,J) + 2.0D0*A2(I,J)*UN(I,J) - 2.0D0*A3(I,J)*UN(I+1,J)
    !            ELSE
    !                RVXN=DT/Re*VISCOUS_UXN(I,J)
    !            END IF
    !
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*UN(I,J-1) + 2.0D0*B2(I,J)*UN(I,J) - 2.0D0*B3(I,J)*UN(I,J+1)
    !            ELSE
    !                RVYN=DT/Re*VISCOUS_UYN(I,J)
    !            END IF
    !
    !            RU(I,J)=DT*(-RPN-1.5D0*RCN+0.5D0*RCN1)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:IM-1),TRIB(2:IM-1),TRIC(2:IM-1),TRID(2:IM-1) )
    !ALLOCATE( TRIC1(2:IM-1),TRID1(2:IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=1,JM-1,1
    !    IF(J==1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=3,IM-2,1
    !            IF(TYPEUX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RU(I,J)
    !            ELSE IF(TYPEUX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=IB_AUX(I,J)
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-1)THEN
    !                TRIA(I)=IB_AUX(I,J)
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            END IF
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO I=3,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO I=3,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,2,-1
    !        DUP(I,J)=TRID1(I)-TRIC1(I)*DUP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!���������ƣ����¸����̶�ֵ��
    !DUP(1       ,1:JM-1:1)=BCU_AL*DUP(2       ,1:JM-1:1)
    !DUP(IM      ,1:JM-1:1)=BCU_AR*DUP(IM-1    ,1:JM-1:1)
    !DUP(:       ,0       )=BCU_AB*DUP(:       ,1       )
    !DUP(:       ,JM      )=BCU_AT*DUP(:       ,JM-1    )
    !
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(JM-1),TRIB(JM-1),TRIC(JM-1),TRID(JM-1) )
    !ALLOCATE( TRIC1(JM-1),TRID1(JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=2,IM-1,1
    !    IF(I==2)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=2,JM-2,1
    !            IF(TYPEUY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DUP(I,J)
    !            ELSE IF(TYPEUY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=IB_BUY(I,J)
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-1)THEN
    !                TRIA(J)=IB_BUY(I,J)
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            END IF
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO J=2,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO J=2,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,1,-1
    !        DUH(I,J)=TRID1(J)-TRIC1(J)*DUH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DUH(1       ,1:JM-1:1)=BCU_AL*DUH(2       ,1:JM-1:1)
    !DUH(IM      ,1:JM-1:1)=BCU_AR*DUH(IM-1    ,1:JM-1:1)
    !DUH(:       ,0       )=BCU_AB*DUH(:       ,1       )
    !DUH(:       ,JM      )=BCU_AT*DUH(:       ,JM-1    )
    !
    !UHAT=UN+DUH
    !VHAT=VN+DVH
    !
    !RETURN
    !END SUBROUTINE
    !
    !!***************************************��AB-CN��ɢ��N-S�������̽���ʱ���ƽ�(�����������⴦��)********************************************!
    !SUBROUTINE TIMEADVANCE1_ME_AB_CN_CONVECT
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !REAL(KIND=8),ALLOCATABLE::RU(:,:),RV(:,:)
    !REAL(KIND=8),ALLOCATABLE::A1(:,:),A2(:,:),A3(:,:),B1(:,:),B2(:,:),B3(:,:)!C1,C2,C3
    !REAL(KIND=8),ALLOCATABLE::DUH(:,:),DVH(:,:)
    !REAL(KIND=8),ALLOCATABLE::DUP(:,:),DVP(:,:)
    !REAL(KIND=8)::UR,UL,UU,UD,VR,VL,VU,VD!������ϳ��ٶȣ�nʱ���
    !REAL(KIND=8)::UR1,UL1,UU1,UD1,VR1,VL1,VU1,VD1!������ϳ��ٶȣ�n-1ʱ���
    !REAL(KIND=8)::RPN,RCXN,RCYN,RCN,RCN1,RVXN,RVYN!�Ҷ��������
    !REAL(KIND=8),ALLOCATABLE::TRIA(:),TRIB(:),TRIC(:),TRID(:)
    !REAL(KIND=8),ALLOCATABLE::TRIC1(:),TRID1(:)
    !
    !REAL(KIND=8)::D1,D2
    !
    !ALLOCATE( RU(IM,0:JM),RV(0:IM,JM),DUP(IM,0:JM),DVP(0:IM,JM),DUH(IM,0:JM),DVH(0:IM,JM) )
    !
    !RU=0.0D0
    !RV=0.0D0
    !DUH=0.0D0
    !DUP=0.0D0
    !DVH=0.0D0
    !DVP=0.0D0
    !
    !!-----------------------------------------------V----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        IF( TYPEVX(I,J)/=-10 .AND. TYPEVY(I,J)/=-10 )THEN
    !            VU=0.5D0*( VN(I,J+1)+VN(I,J) )
    !            VD=0.5D0*( VN(I,J-1)+VN(I,J) )
    !            CALL LINEAR_INTERPOLATION( VN(I+1,J),VR,VN(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
    !            CALL LINEAR_INTERPOLATION( VN(I  ,J),VL,VN(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( UN(I+1,J),UR,UN(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( UN(I  ,J),UL,UN(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            VU1=0.5D0*( VN1(I,J+1)+VN1(I,J) )
    !            VD1=0.5D0*( VN1(I,J-1)+VN1(I,J) )
    !            CALL LINEAR_INTERPOLATION( VN1(I+1,J),VR1,VN1(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
    !            CALL LINEAR_INTERPOLATION( VN1(I  ,J),VL1,VN1(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( UN1(I+1,J),UR1,UN1(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( UN1(I  ,J),UL1,UN1(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            D1=YPU(J-1)-Y(J)
    !            D2=YPU(J)-Y(J)
    !            RPN=( P(I,J)-P(I,J-1) )/( YPU(J)-YPU(J-1) )
    !
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCXN=( UR *VR -UL *VL  )/( X(I+1)-X(I) )
    !            ELSE
    !                RCXN=CONVECT_VXN(I,J)
    !            END IF
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCYN=( D2*D2*VD *VD +(D1*D1-D2*D2)*VN (I,J)*VN (I,J)-D1*D1*VU *VU  )/( D1*D2*(D2-D1) )
    !            ELSE
    !                RCYN=CONVECT_VYN(I,J)
    !            END IF
    !
    !            RCN =RCXN+RCYN
    !
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCXN=( UR1*VR1-UL1*VL1 )/( X(I+1)-X(I) )
    !            ELSE
    !                RCXN=CONVECT_VXN(I,J)
    !            END IF
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCYN=( D2*D2*VD1*VD1+(D1*D1-D2*D2)*VN1(I,J)*VN1(I,J)-D1*D1*VU1*VU1 )/( D1*D2*(D2-D1) )
    !            ELSE
    !                RCYN=CONVECT_VYN(I,J)
    !            END IF
    !
    !            RCN1=RCXN+RCYN
    !
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*VN(I-1,J) + 2.0D0*A2(I,J)*VN(I,J) - 2.0D0*A3(I,J)*VN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_VXN(I,J)
    !            END IF
    !
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*VN(I,J-1) + 2.0D0*B2(I,J)*VN(I,J) - 2.0D0*B3(I,J)*VN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_VYN(I,J)
    !            END IF
    !
    !            RV(I,J)=DT*(-RPN-1.5D0*RCN+0.5D0*RCN1)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(IM-1),TRIB(IM-1),TRIC(IM-1),TRID(IM-1) )
    !ALLOCATE( TRIC1(IM-1),TRID1(IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=2,JM-1,1
    !    IF(J==2)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=2,IM-2,1
    !            IF(TYPEVX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RV(I,J)
    !            ELSE IF(TYPEVX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=IB_AVX(I,J)
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-1)THEN
    !                TRIA(I)=IB_AVX(I,J)
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            END IF
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO I=2,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO I=2,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,1,-1
    !        DVP(I,J)=TRID1(I)-TRIC1(I)*DVP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVP(0       ,2:JM-1:1)=BCV_AL*DVP(1       ,2:JM-1:1)
    !DVP(IM      ,2:JM-1:1)=BCV_AR*DVP(IM-1    ,2:JM-1:1)
    !DVP(:       ,1       )=BCV_AB*DVP(:       ,2       )
    !DVP(:       ,JM      )=BCV_AT*DVP(:       ,JM-1    )
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:JM-1),TRIB(2:JM-1),TRIC(2:JM-1),TRID(2:JM-1) )
    !ALLOCATE( TRIC1(2:JM-1),TRID1(2:JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=1,IM-1,1
    !    IF(I==1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=3,JM-2,1
    !            IF(TYPEVY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DVP(I,J)
    !            ELSE IF(TYPEVY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=IB_BVY(I,J)
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-1)THEN
    !                TRIA(J)=IB_BVY(I,J)
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            END IF
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO J=3,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO J=3,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,2,-1
    !        DVH(I,J)=TRID1(J)-TRIC1(J)*DVH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVH(0       ,2:JM-1:1)=BCV_AL*DVH(1       ,2:JM-1:1)
    !DVH(IM      ,2:JM-1:1)=BCV_AR*DVH(IM-1    ,2:JM-1:1)
    !DVH(:       ,1       )=BCV_AB*DVH(:       ,2       )
    !DVH(:       ,JM      )=BCV_AT*DVH(:       ,JM-1    )
    !
    !!-----------------------------------------------U----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        IF( TYPEUX(I,J)/=-10 .AND. TYPEUY(I,J)/=-10 )THEN
    !            UR=0.5D0*( UN(I+1,J)+UN(I,J) )
    !            UL=0.5D0*( UN(I-1,J)+UN(I,J) )
    !            CALL LINEAR_INTERPOLATION( UN(I,J+1),UU,UN(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !            CALL LINEAR_INTERPOLATION( UN(I,J  ),UD,UN(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J+1),VU,VN(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J  ),VD,VN(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !            UR1=0.5D0*( UN1(I+1,J)+UN1(I,J) )
    !            UL1=0.5D0*( UN1(I-1,J)+UN1(I,J) )
    !            CALL LINEAR_INTERPOLATION( UN1(I,J+1),UU1,UN1(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !            CALL LINEAR_INTERPOLATION( UN1(I,J  ),UD1,UN1(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( VN1(I,J+1),VU1,VN1(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( VN1(I,J  ),VD1,VN1(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !            D1=XPV(I-1)-X(I)
    !            D2=XPV(I)-X(I)
    !            RPN=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCXN=( D2*D2*UL *UL +(D1*D1-D2*D2)*UN (I,J)*UN (I,J)-D1*D1*UR *UR  )/( D1*D2*(D2-D1) )
    !            ELSE
    !                RCXN=CONVECT_UXN(I,J)
    !            END IF
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCYN=( UU *VU -UD *VD  )/( Y(J+1)-Y(J) )
    !            ELSE
    !                RCYN=CONVECT_UYN(I,J)
    !            END IF
    !
    !            RCN =RCXN+RCYN
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCXN=( D2*D2*UL1*UL1+(D1*D1-D2*D2)*UN1(I,J)*UN1(I,J)-D1*D1*UR1*UR1 )/( D1*D2*(D2-D1) )
    !            ELSE
    !                RCXN=CONVECT_UXN(I,J)
    !            END IF
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RCYN=( UU1*VU1-UD1*VD1 )/( Y(J+1)-Y(J) )
    !            ELSE
    !                RCYN=CONVECT_UYN(I,J)
    !            END IF
    !
    !            RCN =RCXN+RCYN
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*UN(I-1,J) + 2.0D0*A2(I,J)*UN(I,J) - 2.0D0*A3(I,J)*UN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_UXN(I,J)
    !            END IF
    !
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*UN(I,J-1) + 2.0D0*B2(I,J)*UN(I,J) - 2.0D0*B3(I,J)*UN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_UYN(I,J)
    !            END IF
    !
    !            RU(I,J)=DT*(-RPN-1.5D0*RCN+0.5D0*RCN1)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:IM-1),TRIB(2:IM-1),TRIC(2:IM-1),TRID(2:IM-1) )
    !ALLOCATE( TRIC1(2:IM-1),TRID1(2:IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=1,JM-1,1
    !    IF(J==1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=3,IM-2,1
    !            IF(TYPEUX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RU(I,J)
    !            ELSE IF(TYPEUX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=IB_AUX(I,J)
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-1)THEN
    !                TRIA(I)=IB_AUX(I,J)
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            END IF
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO I=3,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO I=3,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,2,-1
    !        DUP(I,J)=TRID1(I)-TRIC1(I)*DUP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!���������ƣ����¸����̶�ֵ��
    !DUP(1       ,1:JM-1:1)=BCU_AL*DUP(2       ,1:JM-1:1)
    !DUP(IM      ,1:JM-1:1)=BCU_AR*DUP(IM-1    ,1:JM-1:1)
    !DUP(:       ,0       )=BCU_AB*DUP(:       ,1       )
    !DUP(:       ,JM      )=BCU_AT*DUP(:       ,JM-1    )
    !
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(JM-1),TRIB(JM-1),TRIC(JM-1),TRID(JM-1) )
    !ALLOCATE( TRIC1(JM-1),TRID1(JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=2,IM-1,1
    !    IF(I==2)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=2,JM-2,1
    !            IF(TYPEUY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DUP(I,J)
    !            ELSE IF(TYPEUY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=IB_BUY(I,J)
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-1)THEN
    !                TRIA(J)=IB_BUY(I,J)
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            END IF
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO J=2,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO J=2,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,1,-1
    !        DUH(I,J)=TRID1(J)-TRIC1(J)*DUH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DUH(1       ,1:JM-1:1)=BCU_AL*DUH(2       ,1:JM-1:1)
    !DUH(IM      ,1:JM-1:1)=BCU_AR*DUH(IM-1    ,1:JM-1:1)
    !DUH(:       ,0       )=BCU_AB*DUH(:       ,1       )
    !DUH(:       ,JM      )=BCU_AT*DUH(:       ,JM-1    )
    !
    !UHAT=UN+DUH
    !VHAT=VN+DVH
    !
    !RETURN
    !END SUBROUTINE
    !
    !
    !!***************************************��CN��ɢ��N-S�������̣��޶��������ʱ���ƽ�********************************************!
    !SUBROUTINE TIMEADVANCE1_ME_NOCONVECT_CN
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !REAL(KIND=8),ALLOCATABLE::RU(:,:),RV(:,:)
    !REAL(KIND=8),ALLOCATABLE::A1(:,:),A2(:,:),A3(:,:),B1(:,:),B2(:,:),B3(:,:)!C1,C2,C3
    !REAL(KIND=8),ALLOCATABLE::DUH(:,:),DVH(:,:)
    !REAL(KIND=8),ALLOCATABLE::DUP(:,:),DVP(:,:)
    !REAL(KIND=8)::RPN,RVXN,RVYN!�Ҷ��������
    !REAL(KIND=8),ALLOCATABLE::TRIA(:),TRIB(:),TRIC(:),TRID(:)
    !REAL(KIND=8),ALLOCATABLE::TRIC1(:),TRID1(:)
    !
    !REAL(KIND=8)::D1,D2
    !
    !ALLOCATE( RU(IM,0:JM),RV(0:IM,JM),DUP(IM,0:JM),DVP(0:IM,JM),DUH(IM,0:JM),DVH(0:IM,JM) )
    !
    !RU=0.0D0
    !RV=0.0D0
    !DUH=0.0D0
    !DUP=0.0D0
    !DVH=0.0D0
    !DVP=0.0D0
    !
    !!-----------------------------------------------V----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        IF( TYPEVX(I,J)/=-10 .AND. TYPEVY(I,J)/=-10 )THEN
    !            RPN=( P(I,J)-P(I,J-1) )/( YPU(J)-YPU(J-1) )
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*VN(I-1,J) + 2.0D0*A2(I,J)*VN(I,J) - 2.0D0*A3(I,J)*VN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_VXN(I,J)
    !            END IF
    !
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*VN(I,J-1) + 2.0D0*B2(I,J)*VN(I,J) - 2.0D0*B3(I,J)*VN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_VYN(I,J)
    !            END IF
    !
    !            RV(I,J)=DT*(-RPN)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(IM-1),TRIB(IM-1),TRIC(IM-1),TRID(IM-1) )
    !ALLOCATE( TRIC1(IM-1),TRID1(IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=2,JM-1,1
    !    IF(J==2)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=2,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RV(I,J)
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=2,IM-2,1
    !            IF(TYPEVX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RV(I,J)
    !            ELSE IF(TYPEVX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=IB_AVX(I,J)
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-1)THEN
    !                TRIA(I)=IB_AVX(I,J)
    !                TRIB(I)=1.0D0+IB_CVX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RV(I,J)+IB_RVX(I,J)
    !            ELSE IF(TYPEVX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RVX(I,J)
    !            END IF
    !        END DO
    !
    !        I=1
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RV(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCV_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RV(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO I=2,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO I=2,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,1,-1
    !        DVP(I,J)=TRID1(I)-TRIC1(I)*DVP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVP(0       ,2:JM-1:1)=BCV_AL*DVP(1       ,2:JM-1:1)
    !DVP(IM      ,2:JM-1:1)=BCV_AR*DVP(IM-1    ,2:JM-1:1)
    !DVP(:       ,1       )=BCV_AB*DVP(:       ,2       )
    !DVP(:       ,JM      )=BCV_AT*DVP(:       ,JM-1    )
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:JM-1),TRIB(2:JM-1),TRIC(2:JM-1),TRID(2:JM-1) )
    !ALLOCATE( TRIC1(2:JM-1),TRID1(2:JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=1,IM-1,1
    !    IF(I==1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=3,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DVP(I,J)
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=3,JM-2,1
    !            IF(TYPEVY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DVP(I,J)
    !            ELSE IF(TYPEVY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=IB_BVY(I,J)
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-1)THEN
    !                TRIA(J)=IB_BVY(I,J)
    !                TRIB(J)=1.0D0+IB_CVY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DVP(I,J)+IB_RVY(I,J)
    !            ELSE IF(TYPEVY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RVY(I,J)
    !            END IF
    !        END DO
    !
    !        J=2
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DVP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCV_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DVP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO J=3,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO J=3,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DVH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,2,-1
    !        DVH(I,J)=TRID1(J)-TRIC1(J)*DVH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DVH(0       ,2:JM-1:1)=BCV_AL*DVH(1       ,2:JM-1:1)
    !DVH(IM      ,2:JM-1:1)=BCV_AR*DVH(IM-1    ,2:JM-1:1)
    !DVH(:       ,1       )=BCV_AB*DVH(:       ,2       )
    !DVH(:       ,JM      )=BCV_AT*DVH(:       ,JM-1    )
    !
    !!-----------------------------------------------U----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        IF( TYPEUX(I,J)/=-10 .AND. TYPEUY(I,J)/=-10 )THEN
    !            RPN=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*UN(I-1,J) + 2.0D0*A2(I,J)*UN(I,J) - 2.0D0*A3(I,J)*UN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_UXN(I,J)
    !            END IF
    !
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*UN(I,J-1) + 2.0D0*B2(I,J)*UN(I,J) - 2.0D0*B3(I,J)*UN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_UYN(I,J)
    !            END IF
    !
    !            RU(I,J)=DT*(-RPN)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !!------ADI------!
    !!------XSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(2:IM-1),TRIB(2:IM-1),TRIC(2:IM-1),TRID(2:IM-1) )
    !ALLOCATE( TRIC1(2:IM-1),TRID1(2:IM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO J=1,JM-1,1
    !    IF(J==1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ϱ߽�
    !    ELSE IF(J==JM-1)THEN
    !        DO I=3,IM-2,1
    !            TRIA(I)=A1(I,J)
    !            TRIB(I)=1.0D0-A2(I,J)
    !            TRIC(I)=A3(I,J)
    !            TRID(I)=RU(I,J)
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO I=3,IM-2,1
    !            IF(TYPEUX(I,J)==10)THEN
    !                TRIA(I)=A1(I,J)
    !                TRIB(I)=1.0D0-A2(I,J)
    !                TRIC(I)=A3(I,J)
    !                TRID(I)=RU(I,J)
    !            ELSE IF(TYPEUX(I,J)==0)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==1)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=IB_AUX(I,J)
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-1)THEN
    !                TRIA(I)=IB_AUX(I,J)
    !                TRIB(I)=1.0D0+IB_CUX(I,J)
    !                TRIC(I)=0.0D0
    !                TRID(I)=RU(I,J)+IB_RUX(I,J)
    !            ELSE IF(TYPEUX(I,J)==-10)THEN
    !                TRIA(I)=0.0D0
    !                TRIB(I)=1.0D0
    !                TRIC(I)=0.0D0
    !                TRID(I)=IB_RUX(I,J)
    !            END IF
    !        END DO
    !
    !        I=2
    !        TRIA(I)=0.0D0
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AL*A1(I,J)
    !        TRIC(I)=A3(I,J)
    !        TRID(I)=RU(I,J)
    !
    !        I=IM-1
    !        TRIA(I)=A1(I,J)
    !        TRIB(I)=1.0D0-A2(I,J)+BCU_AR*A3(I,J)
    !        TRIC(I)=0.0D0
    !        TRID(I)=RU(I,J)
    !
    !    END IF
    !
    !    !------FORWARD SWEEP------!
    !    TRIC1(2)=TRIC(2)/TRIB(2)
    !    DO I=3,IM-2,1
    !        TRIC1(I)=TRIC(I)/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    TRID1(2)=TRID(2)/TRIB(2)
    !    DO I=3,IM-1,1
    !        TRID1(I)=( TRID(I)-TRIA(I)*TRID1(I-1) )/( TRIB(I)-TRIA(I)*TRIC1(I-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUP(IM-1,J)=TRID1(IM-1)
    !    DO I=IM-2,2,-1
    !        DUP(I,J)=TRID1(I)-TRIC1(I)*DUP(I+1,J)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!���������ƣ����¸����̶�ֵ��
    !DUP(1       ,1:JM-1:1)=BCU_AL*DUP(2       ,1:JM-1:1)
    !DUP(IM      ,1:JM-1:1)=BCU_AR*DUP(IM-1    ,1:JM-1:1)
    !DUP(:       ,0       )=BCU_AB*DUP(:       ,1       )
    !DUP(:       ,JM      )=BCU_AT*DUP(:       ,JM-1    )
    !
    !!------YSWEEP------!
    !IF( ALLOCATED(TRIA) )THEN
    !    DEALLOCATE(TRIA,TRIB,TRIC,TRID,TRIC1,TRID1)
    !END IF
    !ALLOCATE( TRIA(JM-1),TRIB(JM-1),TRIC(JM-1),TRID(JM-1) )
    !ALLOCATE( TRIC1(JM-1),TRID1(JM-1) )
    !TRIA=0.0D0
    !TRIB=0.0D0
    !TRIC=0.0D0
    !TRID=0.0D0
    !TRIC1=0.0D0
    !TRID1=0.0D0
    !
    !DO I=2,IM-1,1
    !    IF(I==2)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ұ߽�
    !    ELSE IF(I==IM-1)THEN
    !        DO J=2,JM-2,1
    !            TRIA(J)=B1(I,J)
    !            TRIB(J)=1.0D0-B2(I,J)
    !            TRIC(J)=B3(I,J)
    !            TRID(J)=DUP(I,J)
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !        !�ڲ�
    !    ELSE
    !        DO J=2,JM-2,1
    !            IF(TYPEUY(I,J)==10)THEN
    !                TRIA(J)=B1(I,J)
    !                TRIB(J)=1.0D0-B2(I,J)
    !                TRIC(J)=B3(I,J)
    !                TRID(J)=DUP(I,J)
    !            ELSE IF(TYPEUY(I,J)==0)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==1)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=IB_BUY(I,J)
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-1)THEN
    !                TRIA(J)=IB_BUY(I,J)
    !                TRIB(J)=1.0D0+IB_CUY(I,J)
    !                TRIC(J)=0.0D0
    !                TRID(J)=DUP(I,J)+IB_RUY(I,J)
    !            ELSE IF(TYPEUY(I,J)==-10)THEN
    !                TRIA(J)=0.0D0
    !                TRIB(J)=1.0D0
    !                TRIC(J)=0.0D0
    !                TRID(J)=IB_RUY(I,J)
    !            END IF
    !        END DO
    !
    !        J=1
    !        TRIA(J)=0.0D0
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AB*B1(I,J)
    !        TRIC(J)=B3(I,J)
    !        TRID(J)=DUP(I,J)
    !
    !        J=JM-1
    !        TRIA(J)=B1(I,J)
    !        TRIB(J)=1.0D0-B2(I,J)+BCU_AT*B3(I,J)
    !        TRIC(J)=0.0D0
    !        TRID(J)=DUP(I,J)
    !
    !    END IF
    !    !------FORWARD SWEEP------!
    !    TRIC1(1)=TRIC(1)/TRIB(1)
    !    DO J=2,JM-2,1
    !        TRIC1(J)=TRIC(J)/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    TRID1(1)=TRID(1)/TRIB(1)
    !    DO J=2,JM-1,1
    !        TRID1(J)=( TRID(J)-TRIA(J)*TRID1(J-1) )/( TRIB(J)-TRIA(J)*TRIC1(J-1) )
    !    END DO
    !    !------BACK SUBSTITUTION------!
    !    DUH(I,JM-1)=TRID1(JM-1)
    !    DO J=JM-2,1,-1
    !        DUH(I,J)=TRID1(J)-TRIC1(J)*DUH(I,J+1)
    !    END DO
    !
    !END DO
    !!------SWEEP����------!
    !!------ʩ�ӱ߽�����------!
    !DUH(1       ,1:JM-1:1)=BCU_AL*DUH(2       ,1:JM-1:1)
    !DUH(IM      ,1:JM-1:1)=BCU_AR*DUH(IM-1    ,1:JM-1:1)
    !DUH(:       ,0       )=BCU_AB*DUH(:       ,1       )
    !DUH(:       ,JM      )=BCU_AT*DUH(:       ,JM-1    )
    !
    !UHAT=UN+DUH
    !VHAT=VN+DVH
    !
    !RETURN
    !END SUBROUTINE
    !
    !!***************************************����ʽEULER��ɢ��N-S�������̣��޶��������ʱ���ƽ�********************************************!
    !SUBROUTINE TIMEADVANCE1_ME_NOCONVECT_EXPLICIT
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !REAL(KIND=8),ALLOCATABLE::RU(:,:),RV(:,:)
    !REAL(KIND=8),ALLOCATABLE::A1(:,:),A2(:,:),A3(:,:),B1(:,:),B2(:,:),B3(:,:)!C1,C2,C3
    !REAL(KIND=8),ALLOCATABLE::DUH(:,:),DVH(:,:)
    !REAL(KIND=8),ALLOCATABLE::DUP(:,:),DVP(:,:)
    !REAL(KIND=8)::RPN,RVXN,RVYN!�Ҷ��������
    !REAL(KIND=8),ALLOCATABLE::TRIA(:),TRIB(:),TRIC(:),TRID(:)
    !REAL(KIND=8),ALLOCATABLE::TRIC1(:),TRID1(:)
    !
    !REAL(KIND=8)::D1,D2
    !
    !ALLOCATE( RU(IM,0:JM),RV(0:IM,JM),DUP(IM,0:JM),DVP(0:IM,JM),DUH(IM,0:JM),DVH(0:IM,JM) )
    !
    !RU=0.0D0
    !RV=0.0D0
    !DUH=0.0D0
    !DUP=0.0D0
    !DVH=0.0D0
    !DVP=0.0D0
    !
    !!-----------------------------------------------V----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=2,JM-1,1
    !    DO I=1,IM-1,1
    !        IF( TYPEVX(I,J)/=-10 .AND. TYPEVY(I,J)/=-10 )THEN
    !            RPN=( P(I,J)-P(I,J-1) )/( YPU(J)-YPU(J-1) )
    !            IF(TYPEVXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*VN(I-1,J) + 2.0D0*A2(I,J)*VN(I,J) - 2.0D0*A3(I,J)*VN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_VXN(I,J)
    !            END IF
    !
    !            IF(TYPEVYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*VN(I,J-1) + 2.0D0*B2(I,J)*VN(I,J) - 2.0D0*B3(I,J)*VN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_VYN(I,J)
    !            END IF
    !
    !            RV(I,J)=DT*(-RPN)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !DVH=RV
    !!------ʩ�ӱ߽�����------!
    !DVH(0       ,2:JM-1:1)=BCV_AL*DVH(1       ,2:JM-1:1)
    !DVH(IM      ,2:JM-1:1)=BCV_AR*DVH(IM-1    ,2:JM-1:1)
    !DVH(:       ,1       )=BCV_AB*DVH(:       ,2       )
    !DVH(:       ,JM      )=BCV_AT*DVH(:       ,JM-1    )
    !
    !!-----------------------------------------------U----------------------------------------------!
    !
    !!------�����ɢϵ��------!
    !IF( ALLOCATED(A1) )THEN
    !    DEALLOCATE( A1,A2,A3,B1,B2,B3 )
    !END IF
    !ALLOCATE( A1(IM-1,JM-1),A2(IM-1,JM-1),A3(IM-1,JM-1),B1(IM-1,JM-1),B2(IM-1,JM-1),B3(IM-1,JM-1) )
    !A1=0.0D0
    !A2=0.0D0
    !A3=0.0D0
    !B1=0.0D0
    !B2=0.0D0
    !B3=0.0D0
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        A1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
    !        A2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        A3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
    !        B1(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
    !        B2(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !        B3(I,J)=-DT/(2.0D0*Re)*2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
    !    END DO
    !END DO
    !
    !!------RHS------!
    !DO J=1,JM-1,1
    !    DO I=2,IM-1,1
    !        IF( TYPEUX(I,J)/=-10 .AND. TYPEUY(I,J)/=-10 )THEN
    !            RPN=( P(I,J)-P(I-1,J) )/( XPV(I)-XPV(I-1) )
    !
    !            IF(TYPEUXN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVXN=-2.0D0*A1(I,J)*UN(I-1,J) + 2.0D0*A2(I,J)*UN(I,J) - 2.0D0*A3(I,J)*UN(I+1,J)
    !            ELSE
    !                RVXN=VISCOUS_UXN(I,J)
    !            END IF
    !
    !            IF(TYPEUYN(I,J)==10 .OR. VISCOUS_TERM_METHOD==1)THEN
    !                RVYN=-2.0D0*B1(I,J)*UN(I,J-1) + 2.0D0*B2(I,J)*UN(I,J) - 2.0D0*B3(I,J)*UN(I,J+1)
    !            ELSE
    !                RVYN=VISCOUS_UYN(I,J)
    !            END IF
    !
    !            RU(I,J)=DT*(-RPN)+(RVXN+RVYN)
    !        END IF
    !    END DO
    !END DO
    !
    !DUH=RU
    !!------ʩ�ӱ߽�����------!
    !DUH(1       ,1:JM-1:1)=BCU_AL*DUH(2       ,1:JM-1:1)
    !DUH(IM      ,1:JM-1:1)=BCU_AR*DUH(IM-1    ,1:JM-1:1)
    !DUH(:       ,0       )=BCU_AB*DUH(:       ,1       )
    !DUH(:       ,JM      )=BCU_AT*DUH(:       ,JM-1    )
    !
    !UHAT=UN+DUH
    !VHAT=VN+DVH
    !
    !RETURN
    !END SUBROUTINE