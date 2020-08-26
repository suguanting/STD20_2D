    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************ʹ��PHI���µõ��ٶȳ���ѹ���� R-K��ʽ********************************************!
    SUBROUTINE TIMEADVANCE3_UPDATEUP_RK
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::A1,A2,A3,B1,B2,B3,CORRECTION

    !------����VK------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            VK(I,J)=VHAT(I,J)-ALPHA(NSUBSTEP)*DT*( PHI(I,J)-PHI(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    !------ʩ�ӱ߽�����------!
    VK(1:IM-1:1,1       )=BCV_AB*VK(1:IM-1:1,2       )+BCV_BB
    VK(1:IM-1:1,JM      )=BCV_AT*VK(1:IM-1:1,JM-1    )+BCV_BT
    VK(0       ,:       )=BCV_AL*VK(1       ,:       )+BCV_BL
    VK(IM      ,:       )=BCV_AR*VK(IM-1    ,:       )+BCV_BR
    
    
    !------����UK------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            UK(I,J)=UHAT(I,J)-ALPHA(NSUBSTEP)*DT*( PHI(I,J)-PHI(I-1,J) )/( XPV(I)-XPV(I-1) )
        END DO
    END DO
    !------ʩ�ӱ߽�����------!
    UK(1       ,1:JM-1:1)=BCU_AL*UK(2       ,1:JM-1:1)+BCU_BL
    UK(IM      ,1:JM-1:1)=BCU_AR*UK(IM-1    ,1:JM-1:1)+BCU_BR
    UK(:       ,0       )=BCU_AB*UK(:       ,1       )+BCU_BB
    UK(:       ,JM      )=BCU_AT*UK(:       ,JM-1    )+BCU_BT
    
    !!------ʩ�ӳ��ڱ߽�����(����߽�cell�������غ�)------!
    !IF(BCTYPE_L==2)THEN
    !    DO J=1,JM-1,1
    !        UK(1 ,J)=UK(2   ,J)+(X(2 )-X(1   ))/(Y(J+1)-Y(J))*(VK(1   ,J+1)-VK(1   ,J))
    !    END DO
    !END IF
    !IF(BCTYPE_R==2)THEN
    !    DO J=1,JM-1,1
    !        UK(IM,J)=UK(IM-1,J)-(X(IM)-X(IM-1))/(Y(J+1)-Y(J))*(VK(IM-1,J+1)-VK(IM-1,J))
    !    END DO
    !END IF
    !IF(BCTYPE_B==2)THEN
    !    DO I=1,IM-1,1
    !        VK(I,1 )=VK(I,2   )+(Y(2 )-Y(1   ))/(X(I+1)-X(I))*(UK(I+1,1   )-UK(I,1   ))
    !    END DO
    !END IF
    !IF(BCTYPE_T==2)THEN
    !    DO I=1,IM-1,1
    !        VK(I,JM)=VK(I,JM-1)-(Y(JM)-Y(JM-1))/(X(I+1)-X(I))*(UK(I+1,JM-1)-UK(I,JM-1))
    !    END DO
    !END IF
    
    !------ʩ�ӳ��ڱ߽�����(���������ٶȣ�����ȫ�������غ�)------!
    IF(BCTYPE_L==2 .OR. BCTYPE_L==5)THEN
        DO J=1,JM-1,1
            UK(1 ,J)=UHAT(1 ,J)
        END DO
    END IF
    IF(BCTYPE_R==2 .OR. BCTYPE_R==5)THEN
        DO J=1,JM-1,1
            UK(IM,J)=UHAT(IM,J)
        END DO
    END IF
    IF(BCTYPE_B==2 .OR. BCTYPE_B==5)THEN
        DO I=1,IM-1,1
            VK(I,1 )=VHAT(I,1 )
        END DO
    END IF
    IF(BCTYPE_T==2 .OR. BCTYPE_T==5)THEN
        DO I=1,IM-1,1
            VK(I,JM)=VHAT(I,JM)
        END DO
    END IF
    !------ʩ�ӳ��ڱ߽�����(������������)------!
    IF(BCTYPE_L==5)THEN
        DO J=1,JM,1
            VK(0 ,J)=VHAT(0 ,J)
        END DO
    END IF
    IF(BCTYPE_R==5)THEN
        DO J=1,JM,1
            VK(IM,J)=VHAT(IM,J)
        END DO
    END IF
    IF(BCTYPE_B==5)THEN
        DO I=1,IM,1
            UK(I,0 )=UHAT(I,0 )
        END DO
    END IF
    IF(BCTYPE_T==5)THEN
        DO I=1,IM,1
            UK(I,JM)=UHAT(I,JM)
        END DO
    END IF
    


    !------����P------!
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            A1=2.0D0*( XPV(I+1)-XPV(I  ) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            A2=2.0D0*( XPV(I+1)-XPV(I-1) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            A3=2.0D0*( XPV(I  )-XPV(I-1) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            B1=2.0D0*( YPU(J+1)-YPU(J  ) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            B2=2.0D0*( YPU(J+1)-YPU(J-1) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            B3=2.0D0*( YPU(J  )-YPU(J-1) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            P(I,J)=P(I,J)+PHI(I,J)-0.5D0*ALPHA(NSUBSTEP)*DT/Re*( A1*PHI(I-1,J) - A2*PHI(I,J) + A3*PHI(I+1,J) + B1*PHI(I,J-1) - B2*PHI(I,J) + B3*PHI(I,J+1) )
        END DO
    END DO
    !------��ֵê------!
    !ò��ֻ�н�����������Ҫ�ˣ���ΪҪ�ͽ������ѹ�����ϣ�����������ѹ������PHI���õ�������
    IF( DABS(BCPHI_AL*BCPHI_AR*BCPHI_AB*BCPHI_AT-1.0D0)<=CRITERIA )THEN
        I=1
        J=1
        IF(IB_LOCOMOTION==-2)THEN
            !֮ǰ�����������ʽ�����ݽ�����ĸ�ֵ����
            !CORRECTION=P(I,J)-0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*(T-0.5D0*DT)/Re)
            CORRECTION=P(I,J)-0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*T/Re)
            !����һ����ʽ�������κν����������ʽ
            !CORRECTION=SUM(P)/(IM-1)/(JM-1)

            P=P-CORRECTION
            WRITE(1234566,*) "������������",NSTEP,CORRECTION
        END IF
    END IF
    PMAX=MAXVAL(P)
    PLOC=MAXLOC(P)

    RETURN
    END

    !***************************************ʹ��PHI���µõ��ٶȳ���ѹ����������仯�� A-B��ʽ********************************************!
    SUBROUTINE TIMEADVANCE3_UPDATEUP_AB
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::ERRORU(:,:),ERRORV(:,:)!ERRORVELO(:,:)
    REAL(KIND=8)::A1,A2,A3,B1,B2,B3,CORRECTION

    ALLOCATE( ERRORU(IM,0:JM),ERRORV(0:IM,JM) )

    !------����V------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            V(I,J)=VHAT(I,J)-DT*( PHI(I,J)-PHI(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    !------ʩ�ӱ߽�����------!
    V(1:IM-1:1,1       )=BCV_AB*V(1:IM-1:1,2       )+BCV_BB
    V(1:IM-1:1,JM      )=BCV_AT*V(1:IM-1:1,JM-1    )+BCV_BT
    !���³���
    IF(BCU_AB==1.0D0 .AND. BCU_BB==0.0D0 .AND. BCV_AB==1.0D0 .AND. BCV_BB==0.0D0 .AND. BCPHI_AB==0.0D0 .AND. BCPHI_BB==0.0D0)THEN
        DO I=1,IM-1,1
            V(I,1       )=VHAT(I,1       )-DT*( PHI(I,1 )-PHI(I,0   ) )/( YPU(1 )-YPU(0   ) )
        END DO
    END IF
    !���ϳ���
    IF(BCU_AT==1.0D0 .AND. BCU_BT==0.0D0 .AND. BCV_AT==1.0D0 .AND. BCV_BT==0.0D0 .AND. BCPHI_AT==0.0D0 .AND. BCPHI_BT==0.0D0)THEN
        DO I=1,IM-1,1
            V(I,JM      )=VHAT(I,JM      )-DT*( PHI(I,JM)-PHI(I,JM-1) )/( YPU(JM)-YPU(JM-1) )
        END DO
    END IF
    V(0       ,:       )=BCV_AL*V(1       ,:       )+BCV_BL
    V(IM      ,:       )=BCV_AR*V(IM-1    ,:       )+BCV_BR

    !------����仯�����ֵ��ֵn/n-1ʱ���------!
    ERRORV=DABS(VN-V)
    EVMAX=MAXVAL(ERRORV)
    EVLOC=MAXLOC(ERRORV)
    VMAX=MAXVAL(V)
    VLOC=MAXLOC(V)
    VN1=VN
    VN =V

    !------����U------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            U(I,J)=UHAT(I,J)-DT*( PHI(I,J)-PHI(I-1,J) )/( XPV(I)-XPV(I-1) )
        END DO
    END DO
    !------ʩ�ӱ߽�����------!
    U(1       ,1:JM-1:1)=BCU_AL*U(2       ,1:JM-1:1)+BCU_BL
    U(IM      ,1:JM-1:1)=BCU_AR*U(IM-1    ,1:JM-1:1)+BCU_BR
    !�������
    IF(BCU_AL==1.0D0 .AND. BCU_BL==0.0D0 .AND. BCV_AL==1.0D0 .AND. BCV_BL==0.0D0 .AND. BCPHI_AL==0.0D0 .AND. BCPHI_BL==0.0D0)THEN
        DO J=1,JM-1,1
            U(1       ,J)=UHAT(1     ,J)-DT*( PHI(1 ,J)-PHI(0   ,J) )/( XPV(1 )-XPV(0   ) )
        END DO
    END IF
    !���ҳ���
    IF(BCU_AR==1.0D0 .AND. BCU_BR==0.0D0 .AND. BCV_AR==1.0D0 .AND. BCV_BR==0.0D0 .AND. BCPHI_AR==0.0D0 .AND. BCPHI_BR==0.0D0)THEN
        DO J=1,JM-1,1
            U(IM      ,J)=UHAT(IM    ,J)-DT*( PHI(IM,J)-PHI(IM-1,J) )/( XPV(IM)-XPV(IM-1) )
        END DO
    END IF
    U(:       ,0       )=BCU_AB*U(:       ,1       )+BCU_BB
    U(:       ,JM      )=BCU_AT*U(:       ,JM-1    )+BCU_BT

    !------����仯�����ֵ��ֵn/n-1ʱ���------!
    ERRORU=DABS(UN-U)
    EUMAX=MAXVAL(ERRORU)
    EULOC=MAXLOC(ERRORU)
    UMAX=MAXVAL(U)
    ULOC=MAXLOC(U)
    UN1=UN
    UN =U

    !------����P------!
    DO J=1,JM-1,1
        DO I=1,IM-1,1
            A1=2.0D0*( XPV(I+1)-XPV(I  ) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            A2=2.0D0*( XPV(I+1)-XPV(I-1) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            A3=2.0D0*( XPV(I  )-XPV(I-1) )/( ( XPV(I)-XPV(I-1) )*( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I) ) )
            B1=2.0D0*( YPU(J+1)-YPU(J  ) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            B2=2.0D0*( YPU(J+1)-YPU(J-1) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            B3=2.0D0*( YPU(J  )-YPU(J-1) )/( ( YPU(J)-YPU(J-1) )*( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J) ) )
            P(I,J)=P(I,J)+PHI(I,J)-0.5D0*DT/Re*( A1*PHI(I-1,J) - A2*PHI(I,J) + A3*PHI(I+1,J) + B1*PHI(I,J-1) - B2*PHI(I,J) + B3*PHI(I,J+1) )
        END DO
    END DO
    !------��ֵê------!
    IF( DABS(BCPHI_AL*BCPHI_AR*BCPHI_AB*BCPHI_AT-1.0D0)<=CRITERIA )THEN
        I=1
        J=1
        IF(IB_LOCOMOTION==-2)THEN
            !CORRECTION=P(I,J)-0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*(T-0.5D0*DT)/Re)
            CORRECTION=SUM(P)/(IM-1)/(JM-1)
        ELSE
            CORRECTION=P(I,J)
        END IF
        P=P-CORRECTION
        WRITE(1234566,*) NSTEP,CORRECTION
    END IF
    PMAX=MAXVAL(P)
    PLOC=MAXLOC(P)

    !------�����ٶȱ仯�����ֵ(����������)------!
    ERRORVELOMAX=DMAX1( EUMAX,EVMAX )
    VELOMAX=MAX( UMAX,VMAX )

    !------�в��ļ����------!
    WRITE(30,"( I6,(1X,F10.6),2(1X,F10.6),2(1X,F12.8),(1X,F12.8) )") NSTEP,T,UMAX,VMAX,EUMAX,EVMAX,PMAX

    WRITE(40,"( '������',I6,'  ʱ�䣺',(1X,F10.6) )") NSTEP,T
    WRITE(40,"( 'U��',3I6 )") ULOC
    WRITE(40,"( 'V��',3I6 )") VLOC
    WRITE(40,"( 'P��',3I6 )") PLOC
    WRITE(40,"( 'ERRORU��',3I6 )") EULOC
    WRITE(40,"( 'ERRORV��',3I6 )") EVLOC

    RETURN
    END