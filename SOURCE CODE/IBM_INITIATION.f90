    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !****************************IB相关初始化：储存上一层数据，初始化新一层数据******************************!
    SUBROUTINE IBM_INITIATION
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    IF(NSTEP==NSTART)THEN
        !计算刚开始,没有上一层数据，需人为设置或求解
        IF(CASE_TYPE==1)THEN
            !人为设置
            TYPEUXN=10
            TYPEVXN=10
            TYPEUYN=10
            TYPEVYN=10

            IB_ITSCT_UXN=0.0D0
            IB_ITSCT_VXN=0.0D0
            IB_IPSVL_UXN=0.0D0
            IB_IPSVL_VXN=0.0D0
            IB_ITSCT_UYN=0.0D0
            IB_ITSCT_VYN=0.0D0
            IB_IPSVL_UYN=0.0D0
            IB_IPSVL_VYN=0.0D0

            IB_IPSVL_UXVN=0.0D0
            IB_IPSVL_VXUN=0.0D0
            IB_IPSVL_UYVN=0.0D0
            IB_IPSVL_VYUN=0.0D0

        ELSE IF(CASE_TYPE==2)THEN
            !使用上一时间层边界信息求解
            CALL TYPE_IB_REINITIATE
            T=DT*DBLE(NSTEP-1)
            CALL CAL_QUADRIC_2D
            IF( BOUNDARY_EXISTENCE_1==1)CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_1,QUADRIC_KINETIC_1)
            IF( BOUNDARY_EXISTENCE_2==1)CALL INTERSECTION_QUADRIC(QUADRIC_GEOMETRICAL_2,QUADRIC_KINETIC_2)
            T=DT*DBLE(NSTEP)

            TYPEUXN=TYPEUX
            TYPEVXN=TYPEVX
            TYPEUYN=TYPEUY
            TYPEVYN=TYPEVY

            IB_ITSCT_UXN=IB_ITSCT_UX
            IB_ITSCT_VXN=IB_ITSCT_VX
            IB_IPSVL_UXN=IB_IPSVL_UX
            IB_IPSVL_VXN=IB_IPSVL_VX
            IB_ITSCT_UYN=IB_ITSCT_UY
            IB_ITSCT_VYN=IB_ITSCT_VY
            IB_IPSVL_UYN=IB_IPSVL_UY
            IB_IPSVL_VYN=IB_IPSVL_VY

            IB_IPSVL_UXVN=IB_IPSVL_UXV
            IB_IPSVL_VXUN=IB_IPSVL_VXU
            IB_IPSVL_UYVN=IB_IPSVL_UYV
            IB_IPSVL_VYUN=IB_IPSVL_VYU
        END IF
    ELSE
        !直接储存上一层数据
        TYPEUXN=TYPEUX
        TYPEVXN=TYPEVX
        TYPEUYN=TYPEUY
        TYPEVYN=TYPEVY

        IB_ITSCT_UXN=IB_ITSCT_UX
        IB_ITSCT_VXN=IB_ITSCT_VX
        IB_IPSVL_UXN=IB_IPSVL_UX
        IB_IPSVL_VXN=IB_IPSVL_VX
        IB_ITSCT_UYN=IB_ITSCT_UY
        IB_ITSCT_VYN=IB_ITSCT_VY
        IB_IPSVL_UYN=IB_IPSVL_UY
        IB_IPSVL_VYN=IB_IPSVL_VY

        IB_IPSVL_UXVN=IB_IPSVL_UXV
        IB_IPSVL_VXUN=IB_IPSVL_VXU
        IB_IPSVL_UYVN=IB_IPSVL_UYV
        IB_IPSVL_VYUN=IB_IPSVL_VYU
    END IF

    !储存上一层粘性项
    !CALL CAL_VISCOUS_TERM_N
    !!CALL CAL_VISCOUS_TERM_N_2
    !CALL CAL_CONVECT_TERM_N
    !初始化新一层数据
    CALL TYPE_IB_REINITIATE

    RETURN
    END SUBROUTINE


    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !*****************************************初始化新一层数据**************************************************!
    !TYPE置10及-10,
    !IB_R/A/B/C置0
    !IB_ITSCT/IPSVL置0
    !IBN1_IPSVL置0
    SUBROUTINE TYPE_IB_REINITIATE
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    !绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    REAL(KIND=8)::QUADRIC_VALUE

    TYPEUX=10
    IB_RUX=0.0D0
    IB_AUX=0.0D0
    IB_CUX=0.0D0
    TYPEUY=10
    IB_RUY=0.0D0
    IB_BUY=0.0D0
    IB_CUY=0.0D0
    TYPEVX=10
    IB_RVX=0.0D0
    IB_AVX=0.0D0
    IB_CVX=0.0D0
    TYPEVY=10
    IB_RVY=0.0D0
    IB_BVY=0.0D0
    IB_CVY=0.0D0

    IB_ITSCT_UX=0.0D0
    IB_ITSCT_VX=0.0D0
    IB_IPSVL_UX=0.0D0
    IB_IPSVL_VX=0.0D0
    IB_ITSCT_UY=0.0D0
    IB_ITSCT_VY=0.0D0
    IB_IPSVL_UY=0.0D0
    IB_IPSVL_VY=0.0D0

    IB_IPSVL_UXV=0.0D0
    IB_IPSVL_VXU=0.0D0
    IB_IPSVL_UYV=0.0D0
    IB_IPSVL_VYU=0.0D0

    IBN1_IPSVL_UXN=0.0D0
    IBN1_IPSVL_VXN=0.0D0
    IBN1_IPSVL_UYN=0.0D0
    IBN1_IPSVL_VYN=0.0D0

    !-----------边界1-----------!
    IF( BOUNDARY_EXISTENCE_1==1 )THEN
        COX2=QUADRIC_GEOMETRICAL_1(7)
        COY2=QUADRIC_GEOMETRICAL_1(8)
        COXY=QUADRIC_GEOMETRICAL_1(9)
        COX=QUADRIC_GEOMETRICAL_1(10)
        COY=QUADRIC_GEOMETRICAL_1(11)
        COM=QUADRIC_GEOMETRICAL_1(13)

        !------U------!
        DO J=1,JM-1,1
            DO I=2,IM-1,1
                QUADRIC_VALUE=COX2*X(I)**2.0D0+COY2*YPU(J)**2.0D0+COXY*X(I)*YPU(J)+COX*X(I)+COY*YPU(J)+COM
                IF( QUADRIC_VALUE<0.0D0 )THEN
                    TYPEUX(I,J)=-10
                    TYPEUY(I,J)=-10
                END IF
            END DO
        END DO

        !------V------!
        DO J=2,JM-1,1
            DO I=1,IM-1,1
                QUADRIC_VALUE=COX2*XPV(I)**2.0D0+COY2*Y(J)**2.0D0+COXY*XPV(I)*Y(J)+COX*XPV(I)+COY*Y(J)+COM
                IF( QUADRIC_VALUE<0.0D0 )THEN
                    TYPEVX(I,J)=-10
                    TYPEVY(I,J)=-10
                END IF
            END DO
        END DO

    END IF

    !-----------边界2-----------!
    IF( BOUNDARY_EXISTENCE_2==1 )THEN
        COX2=QUADRIC_GEOMETRICAL_2(7)
        COY2=QUADRIC_GEOMETRICAL_2(8)
        COXY=QUADRIC_GEOMETRICAL_2(9)
        COX=QUADRIC_GEOMETRICAL_2(10)
        COY=QUADRIC_GEOMETRICAL_2(11)
        COM=QUADRIC_GEOMETRICAL_2(13)

        !------U------!
        DO J=1,JM-1,1
            DO I=2,IM-1,1
                QUADRIC_VALUE=COX2*X(I)**2.0D0+COY2*YPU(J)**2.0D0+COXY*X(I)*YPU(J)+COX*X(I)+COY*YPU(J)+COM
                IF( QUADRIC_VALUE<0.0D0 )THEN
                    TYPEUX(I,J)=-10
                    TYPEUY(I,J)=-10
                END IF
            END DO
        END DO

        !------V------!
        DO J=2,JM-1,1
            DO I=1,IM-1,1
                QUADRIC_VALUE=COX2*XPV(I)**2.0D0+COY2*Y(J)**2.0D0+COXY*XPV(I)*Y(J)+COX*XPV(I)+COY*Y(J)+COM
                IF( QUADRIC_VALUE<0.0D0 )THEN
                    TYPEVX(I,J)=-10
                    TYPEVY(I,J)=-10
                END IF
            END DO
        END DO

    END IF

    RETURN
    END SUBROUTINE

    !!******************************************由TYPEN和UN求解VISCOUS_N*****************************************************!
    !SUBROUTINE CAL_VISCOUS_TERM_N
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !
    !REAL(KIND=8)::DB,DE,IB_B,IB_I,IB_IE
    !
    !VISCOUS_VXN=0.0D0
    !VISCOUS_VYN=0.0D0
    !VISCOUS_UXN=0.0D0
    !VISCOUS_UYN=0.0D0
    !
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !
    !    IF(TYPEVXN(I,J)==1)THEN
    !        DE=XPV(I+1)-XPV(I)
    !        DB=XPV(I)-IB_ITSCT_VX(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_VXN(I,J)=IB_B*IB_IPSVL_VX(I,J)+IB_I*VN(I,J)+IB_IE*VN(I+1,J)
    !    ELSE IF(TYPEVXN(I,J)==-1)THEN
    !        DE=-XPV(I-1)+XPV(I)
    !        DB=-XPV(I)+IB_ITSCT_VX(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_VXN(I,J)=IB_B*IB_IPSVL_VX(I,J)+IB_I*VN(I,J)+IB_IE*VN(I-1,J)
    !    ELSE IF(TYPEVXN(I,J)==0)THEN
    !        VISCOUS_VXN(I,J)=0.0D0
    !    ELSE IF(TYPEVXN(I,J)==-10)THEN
    !        VISCOUS_VXN(I,J)=0.0D0
    !    END IF
    !
    !    IF(TYPEVYN(I,J)==1)THEN
    !        DE=Y(J+1)-Y(J)
    !        DB=Y(J)-IB_ITSCT_VY(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_VYN(I,J)=IB_B*IB_IPSVL_VY(I,J)+IB_I*VN(I,J)+IB_IE*VN(I,J+1)
    !    ELSE IF(TYPEVYN(I,J)==-1)THEN
    !        DE=-Y(J-1)+Y(J)
    !        DB=-Y(J)+IB_ITSCT_VY(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_VYN(I,J)=IB_B*IB_IPSVL_VY(I,J)+IB_I*VN(I,J)+IB_IE*VN(I,J-1)
    !    ELSE IF(TYPEVYN(I,J)==0)THEN
    !        VISCOUS_VYN(I,J)=0.0D0
    !    ELSE IF(TYPEVYN(I,J)==-10)THEN
    !        VISCOUS_VYN(I,J)=0.0D0
    !    END IF
    !
    !    IF(TYPEUXN(I,J)==1)THEN
    !        DE=X(I+1)-X(I)
    !        DB=X(I)-IB_ITSCT_UX(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_UXN(I,J)=IB_B*IB_IPSVL_UX(I,J)+IB_I*UN(I,J)+IB_IE*UN(I+1,J)
    !    ELSE IF(TYPEUXN(I,J)==-1)THEN
    !        DE=-X(I-1)+X(I)
    !        DB=-X(I)+IB_ITSCT_UX(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_UXN(I,J)=IB_B*IB_IPSVL_UX(I,J)+IB_I*UN(I,J)+IB_IE*UN(I-1,J)
    !    ELSE IF(TYPEUXN(I,J)==0)THEN
    !        VISCOUS_UXN(I,J)=0.0D0
    !    ELSE IF(TYPEUXN(I,J)==-10)THEN
    !        VISCOUS_UXN(I,J)=0.0D0
    !    END IF
    !
    !    IF(TYPEUYN(I,J)==1)THEN
    !        DE=YPU(J+1)-YPU(J)
    !        DB=YPU(J)-IB_ITSCT_UY(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_UYN(I,J)=IB_B*IB_IPSVL_UY(I,J)+IB_I*UN(I,J)+IB_IE*UN(I,J+1)
    !    ELSE IF(TYPEUYN(I,J)==-1)THEN
    !        DE=-YPU(J-1)+YPU(J)
    !        DB=-YPU(J)+IB_ITSCT_UY(I,J)
    !        IB_B = 2.0D0/( DB*(DB+DE) )
    !        IB_I =-2.0D0/( DB*DE )
    !        IB_IE= 2.0D0/( DE*(DB+DE) )
    !        VISCOUS_UYN(I,J)=IB_B*IB_IPSVL_UY(I,J)+IB_I*UN(I,J)+IB_IE*UN(I,J-1)
    !    ELSE IF(TYPEUYN(I,J)==0)THEN
    !        VISCOUS_UYN(I,J)=0.0D0
    !    ELSE IF(TYPEUYN(I,J)==-10)THEN
    !        VISCOUS_UYN(I,J)=0.0D0
    !    END IF
    !
    !    END DO
    !END DO
    !
    !RETURN
    !END SUBROUTINE
    !
    !!******************************************由TYPEN和UN求解VISCOUS_N*****************************************************!
    !SUBROUTINE CAL_VISCOUS_TERM_N_2
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !
    !REAL(KIND=8)::DB,DE,IB_B,IB_I,IB_IE
    !
    !VISCOUS_VXN=0.0D0
    !VISCOUS_VYN=0.0D0
    !VISCOUS_UXN=0.0D0
    !VISCOUS_UYN=0.0D0
    !
    !RETURN
    !END SUBROUTINE
    !
    !!******************************************由TYPEN和UN求解CONVECT_N*****************************************************!
    !SUBROUTINE CAL_CONVECT_TERM_N
    !USE DECLARATION
    !USE IMMERSED_BOUNDARY
    !IMPLICIT NONE
    !
    !REAL(KIND=8)::D1,D2,IB_1,IB_0,IB_2
    !REAL(KIND=8)::UR,UL,UU,UD,VR,VL,VU,VD!各方向合成速度，n时间层
    !
    !CONVECT_VXN=0.0D0
    !CONVECT_VYN=0.0D0
    !CONVECT_UXN=0.0D0
    !CONVECT_UYN=0.0D0
    !
    !DO J=1,JM-1,1
    !    DO I=1,IM-1,1
    !
    !        IF(TYPEVXN(I,J)==1)THEN
    !            CALL LINEAR_INTERPOLATION( VN(I+1,J),VR,VN(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
    !            CALL LINEAR_INTERPOLATION( UN(I+1,J),UR,UN(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CONVECT_VXN(I,J)=( UR*VR-IB_IPSVL_VX(I,J)*IB_IPSVL_VXU(I,J) )/( X(I+1)-IB_ITSCT_VX(I,J) )
    !        ELSE IF(TYPEVXN(I,J)==-1)THEN
    !            CALL LINEAR_INTERPOLATION( VN(I  ,J),VL,VN(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
    !            CALL LINEAR_INTERPOLATION( UN(I  ,J),UL,UN(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CONVECT_VXN(I,J)=( UL*VL-IB_IPSVL_VX(I,J)*IB_IPSVL_VXU(I,J) )/( X(I)-IB_ITSCT_VX(I,J) )
    !        ELSE IF(TYPEVXN(I,J)==0)THEN
    !            CONVECT_VXN(I,J)=0.0D0
    !        ELSE IF(TYPEVXN(I,J)==-10)THEN
    !            CONVECT_VXN(I,J)=0.0D0
    !        END IF
    !
    !        IF(TYPEVYN(I,J)==1)THEN
    !            VU=0.5D0*( VN(I,J+1)+VN(I,J) )
    !            CONVECT_VYN(I,J)=( VU**2.0D0-IB_IPSVL_VY(I,J)**2.0D0 )/( YPU(J)-IB_ITSCT_VY(I,J) )
    !        ELSE IF(TYPEVYN(I,J)==-1)THEN
    !            VD=0.5D0*( VN(I,J-1)+VN(I,J) )
    !            CONVECT_VYN(I,J)=( VD**2.0D0-IB_IPSVL_VY(I,J)**2.0D0 )/( YPU(J-1)-IB_ITSCT_VY(I,J) )
    !        ELSE IF(TYPEVYN(I,J)==0)THEN
    !            CONVECT_VYN(I,J)=0.0D0
    !        ELSE IF(TYPEVYN(I,J)==-10)THEN
    !            CONVECT_VYN(I,J)=0.0D0
    !        END IF
    !
    !
    !        IF(TYPEUXN(I,J)==1)THEN
    !            UR=0.5D0*( UN(I+1,J)+UN(I,J) )
    !            CONVECT_UXN(I,J)=( UR**2.0D0-IB_IPSVL_UX(I,J)**2.0D0 )/( XPV(I)-IB_ITSCT_UX(I,J) )
    !        ELSE IF(TYPEUXN(I,J)==-1)THEN
    !            UL=0.5D0*( UN(I-1,J)+UN(I,J) )
    !            CONVECT_UXN(I,J)=( UL**2.0D0-IB_IPSVL_UX(I,J)**2.0D0 )/( XPV(I-1)-IB_ITSCT_UX(I,J) )
    !        ELSE IF(TYPEUXN(I,J)==0)THEN
    !            CONVECT_UXN(I,J)=0.0D0
    !        ELSE IF(TYPEUXN(I,J)==-10)THEN
    !            CONVECT_UXN(I,J)=0.0D0
    !        END IF
    !
    !        IF(TYPEUYN(I,J)==1)THEN
    !            CALL LINEAR_INTERPOLATION( UN(I,J+1),UU,UN(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J+1),VU,VN(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
    !            CONVECT_UYN(I,J)=( UU*VU-IB_IPSVL_UY(I,J)*IB_IPSVL_UYV(I,J) )/( Y(J+1)-IB_ITSCT_UY(I,J) )
    !        ELSE IF(TYPEUYN(I,J)==-1)THEN
    !            CALL LINEAR_INTERPOLATION( UN(I,J  ),UD,UN(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
    !            CALL LINEAR_INTERPOLATION( VN(I,J  ),VD,VN(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
    !            CONVECT_UYN(I,J)=( UD*VD-IB_IPSVL_UY(I,J)*IB_IPSVL_UYV(I,J) )/( Y(J)-IB_ITSCT_UY(I,J) )
    !        ELSE IF(TYPEUYN(I,J)==0)THEN
    !            CONVECT_UYN(I,J)=0.0D0
    !        ELSE IF(TYPEUYN(I,J)==-10)THEN
    !            CONVECT_UYN(I,J)=0.0D0
    !        END IF
    !
    !    END DO
    !END DO
    !
    !RETURN
    !END SUBROUTINE