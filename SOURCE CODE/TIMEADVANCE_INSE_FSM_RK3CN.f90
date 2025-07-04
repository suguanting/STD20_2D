    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !#                  使用FRACTIONAL STEP METHOD时间推进                #!
    !#             RUNGE-KUTTA-3-CRANK-NICOLSON离散的不可压N-S方程        #!
    !#                                                                    #!
    !######################################################################!
    SUBROUTINE TIMEADVANCE_INSE_FSM_RK3CN
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
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
        !!------求解该时间步左端项系数和移到右端项的大小------!
        !!运动边界算法中被废除
        !CALL IBM_PRIMITIVE2DERIVATIVE_TRUECARTESIAN

        !------对RK-CN离散的动量方程进行时间推进得到UHAT------!
        CALL TIMEADVANCE1_ME_RK_CN
        !------求解压力泊松方程得到PHI------!
        CALL TIMEADVANCE2_PPE
        !------使用PHI更新得到速度场和压力场------!
        CALL TIMEADVANCE3_UPDATEUP_RK

        !------计算该时间层IB右端项------!
        CALL CAL_NONFLUIDIC_CONVECT_K
        CALL CAL_NONFLUIDIC_LAPLACE_K

        !------下一循环------!
        UM2=UM1
        VM2=VM1
        UM1=UM
        VM1=VM
        !m=1虚拟时间层结束时赋值TYPEM1，IB_ITSCT_XXM1和IB_IPSVL_XXM1
        IF(NSUBSTEP==1)THEN
            TYPEUXM1=TYPEUX
            TYPEVXM1=TYPEVX
            TYPEUYM1=TYPEUY
            TYPEVYM1=TYPEVY

            IB_ITSCT_UXM1=IB_ITSCT_UX
            IB_ITSCT_VXM1=IB_ITSCT_VX
            IB_IPSVL_UXM1=IB_IPSVL_UX
            IB_IPSVL_VXM1=IB_IPSVL_VX
            IB_ITSCT_UYM1=IB_ITSCT_UY
            IB_ITSCT_VYM1=IB_ITSCT_VY
            IB_IPSVL_UYM1=IB_IPSVL_UY
            IB_IPSVL_VYM1=IB_IPSVL_VY
        END IF

    END DO

    !------循环结束，更新速度------!
    U=UM
    V=VM

    !------求解气动参数------!
    IF( MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NCLCT) ) )==0 )THEN
        IF( BOUNDARY_EXISTENCE_1==1 )THEN
            CALL CAL_CLCT_SRFC_INTGRTN_EXACT(QUADRIC_GEOMETRICAL_1,&
                QUADRIC_KINETIC_1,1)
            !!CALL CAL_CLCT_SRFC_INTGRTN_ALL(QUADRIC_GEOMETRICAL_1,&
            !    QUADRIC_KINETIC_1,1)
            !CALL CAL_CLCT_SRFC_INTGRTN_SOME(QUADRIC_GEOMETRICAL_1,&
            !    QUADRIC_KINETIC_1,1)
        END IF
        IF( BOUNDARY_EXISTENCE_2==1 )THEN
            CALL CAL_CLCT_SRFC_INTGRTN_EXACT(QUADRIC_GEOMETRICAL_2,&
                QUADRIC_KINETIC_2,2)
            !CALL CAL_CLCT_SRFC_INTGRTN_ALL(QUADRIC_GEOMETRICAL_2,&
            !    QUADRIC_KINETIC_2,2)
            !CALL CAL_CLCT_SRFC_INTGRTN_SOME(QUADRIC_GEOMETRICAL_2,&
            !    QUADRIC_KINETIC_2,2)
        END IF
        CALL CAL_CLCT_CONTROL_VOLUME
    END IF

    !------计算变化量最大值赋值n时间层------!
    ERRORV=DABS(VN-V)
    EVMAX=MAXVAL(ERRORV)
    EVLOC=MAXLOC(ERRORV)
    VMAX=MAXVAL(V)
    VLOC=MAXLOC(V)
    VN =V
    !------计算变化量最大值赋值n时间层------!
    ERRORU=DABS(UN-U)
    EUMAX=MAXVAL(ERRORU)
    EULOC=MAXLOC(ERRORU)
    UMAX=MAXVAL(U)
    ULOC=MAXLOC(U)
    UN =U


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
    END SUBROUTINE