    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************使用PHI更新得到速度场和压力场 R-K格式********************************************!
    SUBROUTINE TIMEADVANCE3_UPDATEUP_RK
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::A1,A2,A3,B1,B2,B3,CORRECTION

    !------更新VM------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            VM(I,J)=VHAT(I,J)-ALPHA(NSUBSTEP)*DT*( PHI(I,J)-PHI(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    !------施加边界条件------!
    VM(1:IM-1:1,1       )=BCV_AB*VM(1:IM-1:1,2       )+BCV_BB
    VM(1:IM-1:1,JM      )=BCV_AT*VM(1:IM-1:1,JM-1    )+BCV_BT
    VM(0       ,:       )=BCV_AL*VM(1       ,:       )+BCV_BL
    VM(IM      ,:       )=BCV_AR*VM(IM-1    ,:       )+BCV_BR
    
    
    !------更新UM------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            UM(I,J)=UHAT(I,J)-ALPHA(NSUBSTEP)*DT*( PHI(I,J)-PHI(I-1,J) )/( XPV(I)-XPV(I-1) )
        END DO
    END DO
    !------施加边界条件------!
    UM(1       ,1:JM-1:1)=BCU_AL*UM(2       ,1:JM-1:1)+BCU_BL
    UM(IM      ,1:JM-1:1)=BCU_AR*UM(IM-1    ,1:JM-1:1)+BCU_BR
    UM(:       ,0       )=BCU_AB*UM(:       ,1       )+BCU_BB
    UM(:       ,JM      )=BCU_AT*UM(:       ,JM-1    )+BCU_BT
    
    !!------施加出口边界条件(满足边界cell的质量守恒)------!
    !IF(BCTYPE_L==2)THEN
    !    DO J=1,JM-1,1
    !        UM(1 ,J)=UM(2   ,J)+(X(2 )-X(1   ))/(Y(J+1)-Y(J))*(VM(1   ,J+1)-VM(1   ,J))
    !    END DO
    !END IF
    !IF(BCTYPE_R==2)THEN
    !    DO J=1,JM-1,1
    !        UM(IM,J)=UM(IM-1,J)-(X(IM)-X(IM-1))/(Y(J+1)-Y(J))*(VM(IM-1,J+1)-VM(IM-1,J))
    !    END DO
    !END IF
    !IF(BCTYPE_B==2)THEN
    !    DO I=1,IM-1,1
    !        VM(I,1 )=VM(I,2   )+(Y(2 )-Y(1   ))/(X(I+1)-X(I))*(UM(I+1,1   )-UM(I,1   ))
    !    END DO
    !END IF
    !IF(BCTYPE_T==2)THEN
    !    DO I=1,IM-1,1
    !        VM(I,JM)=VM(I,JM-1)-(Y(JM)-Y(JM-1))/(X(I+1)-X(I))*(UM(I+1,JM-1)-UM(I,JM-1))
    !    END DO
    !END IF
    
    !------施加出口边界条件(修正出口速度，满足全域质量守恒)------!
    IF(BCTYPE_L==2 .OR. BCTYPE_L==5)THEN
        DO J=1,JM-1,1
            UM(1 ,J)=UHAT(1 ,J)
        END DO
    END IF
    IF(BCTYPE_R==2 .OR. BCTYPE_R==5)THEN
        DO J=1,JM-1,1
            UM(IM,J)=UHAT(IM,J)
        END DO
    END IF
    IF(BCTYPE_B==2 .OR. BCTYPE_B==5)THEN
        DO I=1,IM-1,1
            VM(I,1 )=VHAT(I,1 )
        END DO
    END IF
    IF(BCTYPE_T==2 .OR. BCTYPE_T==5)THEN
        DO I=1,IM-1,1
            VM(I,JM)=VHAT(I,JM)
        END DO
    END IF
    !------施加出口边界条件(对流边条额外)------!
    IF(BCTYPE_L==5)THEN
        DO J=1,JM,1
            VM(0 ,J)=VHAT(0 ,J)
        END DO
    END IF
    IF(BCTYPE_R==5)THEN
        DO J=1,JM,1
            VM(IM,J)=VHAT(IM,J)
        END DO
    END IF
    IF(BCTYPE_B==5)THEN
        DO I=1,IM,1
            UM(I,0 )=UHAT(I,0 )
        END DO
    END IF
    IF(BCTYPE_T==5)THEN
        DO I=1,IM,1
            UM(I,JM)=UHAT(I,JM)
        END DO
    END IF
    


    !------更新P------!
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
    !------数值锚------!
    !貌似只有解析解算例需要了，因为要和解析解的压力对上，其他算例的压力早在PHI即得到修正了
    IF( DABS(BCPHI_AL*BCPHI_AR*BCPHI_AB*BCPHI_AT-1.0D0)<=CRITERIA )THEN
        I=1
        J=1
        IF(IB_LOCOMOTION==-2)THEN
            !之前错误的两种形式：根据解析解的赋值方法
            !CORRECTION=P(I,J)-0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*(T-0.5D0*DT)/Re)
            CORRECTION=P(I,J)-0.25D0*( DCOS( 2.0D0*XPV(I) )+DCOS( 2.0D0*YPU(J) ) )*DEXP(-4.0D0*T/Re)
            !另外一种形式，普适任何解析解给定方式
            !CORRECTION=SUM(P)/(IM-1)/(JM-1)

            P=P-CORRECTION
            WRITE(1234566,*) "解析解修正：",NSTEP,CORRECTION
        END IF
    END IF
    PMAX=MAXVAL(P)
    PLOC=MAXLOC(P)

    RETURN
    END

    !***************************************使用PHI更新得到速度场和压力场并计算变化量 A-B格式********************************************!
    SUBROUTINE TIMEADVANCE3_UPDATEUP_AB
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8),ALLOCATABLE::ERRORU(:,:),ERRORV(:,:)!ERRORVELO(:,:)
    REAL(KIND=8)::A1,A2,A3,B1,B2,B3,CORRECTION

    ALLOCATE( ERRORU(IM,0:JM),ERRORV(0:IM,JM) )

    !------更新V------!
    DO J=2,JM-1,1
        DO I=1,IM-1,1
            V(I,J)=VHAT(I,J)-DT*( PHI(I,J)-PHI(I,J-1) )/( YPU(J)-YPU(J-1) )
        END DO
    END DO
    !------施加边界条件------!
    V(1:IM-1:1,1       )=BCV_AB*V(1:IM-1:1,2       )+BCV_BB
    V(1:IM-1:1,JM      )=BCV_AT*V(1:IM-1:1,JM-1    )+BCV_BT
    !若下出口
    IF(BCU_AB==1.0D0 .AND. BCU_BB==0.0D0 .AND. BCV_AB==1.0D0 .AND. BCV_BB==0.0D0 .AND. BCPHI_AB==0.0D0 .AND. BCPHI_BB==0.0D0)THEN
        DO I=1,IM-1,1
            V(I,1       )=VHAT(I,1       )-DT*( PHI(I,1 )-PHI(I,0   ) )/( YPU(1 )-YPU(0   ) )
        END DO
    END IF
    !若上出口
    IF(BCU_AT==1.0D0 .AND. BCU_BT==0.0D0 .AND. BCV_AT==1.0D0 .AND. BCV_BT==0.0D0 .AND. BCPHI_AT==0.0D0 .AND. BCPHI_BT==0.0D0)THEN
        DO I=1,IM-1,1
            V(I,JM      )=VHAT(I,JM      )-DT*( PHI(I,JM)-PHI(I,JM-1) )/( YPU(JM)-YPU(JM-1) )
        END DO
    END IF
    V(0       ,:       )=BCV_AL*V(1       ,:       )+BCV_BL
    V(IM      ,:       )=BCV_AR*V(IM-1    ,:       )+BCV_BR

    !------计算变化量最大值赋值n/n-1时间层------!
    ERRORV=DABS(VN-V)
    EVMAX=MAXVAL(ERRORV)
    EVLOC=MAXLOC(ERRORV)
    VMAX=MAXVAL(V)
    VLOC=MAXLOC(V)
    VN1=VN
    VN =V

    !------更新U------!
    DO J=1,JM-1,1
        DO I=2,IM-1,1
            U(I,J)=UHAT(I,J)-DT*( PHI(I,J)-PHI(I-1,J) )/( XPV(I)-XPV(I-1) )
        END DO
    END DO
    !------施加边界条件------!
    U(1       ,1:JM-1:1)=BCU_AL*U(2       ,1:JM-1:1)+BCU_BL
    U(IM      ,1:JM-1:1)=BCU_AR*U(IM-1    ,1:JM-1:1)+BCU_BR
    !若左出口
    IF(BCU_AL==1.0D0 .AND. BCU_BL==0.0D0 .AND. BCV_AL==1.0D0 .AND. BCV_BL==0.0D0 .AND. BCPHI_AL==0.0D0 .AND. BCPHI_BL==0.0D0)THEN
        DO J=1,JM-1,1
            U(1       ,J)=UHAT(1     ,J)-DT*( PHI(1 ,J)-PHI(0   ,J) )/( XPV(1 )-XPV(0   ) )
        END DO
    END IF
    !若右出口
    IF(BCU_AR==1.0D0 .AND. BCU_BR==0.0D0 .AND. BCV_AR==1.0D0 .AND. BCV_BR==0.0D0 .AND. BCPHI_AR==0.0D0 .AND. BCPHI_BR==0.0D0)THEN
        DO J=1,JM-1,1
            U(IM      ,J)=UHAT(IM    ,J)-DT*( PHI(IM,J)-PHI(IM-1,J) )/( XPV(IM)-XPV(IM-1) )
        END DO
    END IF
    U(:       ,0       )=BCU_AB*U(:       ,1       )+BCU_BB
    U(:       ,JM      )=BCU_AT*U(:       ,JM-1    )+BCU_BT

    !------计算变化量最大值赋值n/n-1时间层------!
    ERRORU=DABS(UN-U)
    EUMAX=MAXVAL(ERRORU)
    EULOC=MAXLOC(ERRORU)
    UMAX=MAXVAL(U)
    ULOC=MAXLOC(U)
    UN1=UN
    UN =U

    !------更新P------!
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
    !------数值锚------!
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

    !------计算速度变化量最大值(无物理意义)------!
    ERRORVELOMAX=DMAX1( EUMAX,EVMAX )
    VELOMAX=MAX( UMAX,VMAX )

    !------残差文件输出------!
    WRITE(30,"( I6,(1X,F10.6),2(1X,F10.6),2(1X,F12.8),(1X,F12.8) )") NSTEP,T,UMAX,VMAX,EUMAX,EVMAX,PMAX

    WRITE(40,"( '步数：',I6,'  时间：',(1X,F10.6) )") NSTEP,T
    WRITE(40,"( 'U：',3I6 )") ULOC
    WRITE(40,"( 'V：',3I6 )") VLOC
    WRITE(40,"( 'P：',3I6 )") PLOC
    WRITE(40,"( 'ERRORU：',3I6 )") EULOC
    WRITE(40,"( 'ERRORV：',3I6 )") EVLOC

    RETURN
    END