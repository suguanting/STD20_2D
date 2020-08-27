    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !#                  ʹ��FRACTIONAL STEP METHODʱ���ƽ�                #!
    !#             RUNGE-KUTTA-3-CRANK-NICOLSON��ɢ�Ĳ���ѹN-S����        #!
    !#                                                                    #!
    !######################################################################!
    SUBROUTINE TIMEADVANCE_INSE_FSM_RK3CN
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::ERRORU(:,:),ERRORV(:,:)!ERRORVELO(:,:)
    ALLOCATE( ERRORU(IM,0:JM),ERRORV(0:IM,JM) )


    GAMA(1)=8.0D0/15.0D0
    GAMA(2)=5.0D0/12.0D0
    GAMA(3)=3.0D0/4.0D0
    RHO(1)=0.0
    RHO(2)=-17.0D0/60.0D0
    RHO(3)=-5.0D0/12.0D0
    ALPHA(1)=8.0D0/15.0D0
    ALPHA(2)=2.0D0/15.0D0
    ALPHA(3)=1.0D0/3.0D0

    UM=0.0D0
    VM=0.0D0
    UM1=UN
    VM1=VN
    UM2=0.0D0
    VM2=0.0D0


    DO NSUBSTEP=1,3,1
        !------����ʱ�䲽�����ϵ�����Ƶ��Ҷ���Ĵ�С------!
        CALL IBM_PRIMITIVE2DERIVATIVE_TRUECARTESIAN

        !------��RK-CN��ɢ�Ķ������̽���ʱ���ƽ��õ�UHAT------!
        CALL TIMEADVANCE1_ME_RK_CN
        !------���ѹ�����ɷ��̵õ�PHI------!
        CALL TIMEADVANCE2_PPE
        !------ʹ��PHI���µõ��ٶȳ���ѹ����------!
        CALL TIMEADVANCE3_UPDATEUP_RK

        !------�����ʱ���IB�Ҷ���------!
        CALL CAL_NONFLUIDIC_CONVECT_K
        CALL CAL_NONFLUIDIC_LAPLACE_K

        !------��һѭ��------!
        UM2=UM1
        VM2=VM1
        UM1=UM
        VM1=VM

    END DO

    !------ѭ�������������ٶ�------!
    U=UM
    V=VM

    !------�����������------!
    IF( MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        IF(IB_LOCOMOTION==-1)THEN
            CALL CAL_CLCT_CONTROL_VOLUME_SIMPLIFIED
            CALL CAL_CLCT_2DCURVE
        END IF
        CALL CAL_CLCT_CONTROL_VOLUME
    END IF

    !------����仯�����ֵ��ֵnʱ���------!
    ERRORV=DABS(VN-V)
    EVMAX=MAXVAL(ERRORV)
    EVLOC=MAXLOC(ERRORV)
    VMAX=MAXVAL(V)
    VLOC=MAXLOC(V)
    VN =V
    !------����仯�����ֵ��ֵnʱ���------!
    ERRORU=DABS(UN-U)
    EUMAX=MAXVAL(ERRORU)
    EULOC=MAXLOC(ERRORU)
    UMAX=MAXVAL(U)
    ULOC=MAXLOC(U)
    UN =U


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
    END SUBROUTINE