    !########################################################################!
    !#                                                                      #!
    !#                              步骤子函数                              #!
    !#                                                                      #!
    !#                       纯笛卡尔方法相关的子程序                       #!
    !#         计算使用上时间层（物理或虚拟）信息，在时间推进前计算         #!
    !########################################################################!

    !!****************************进行一个方向的速度探测******************************!
    !SUBROUTINE VELOCITY_PROBING(VELO_E,VELO_IB,VELO_B,VELO_PROBE,POSI_E,POSI_IB,POSI_B,POSI_PROBE)
    !IMPLICIT NONE
    !
    !REAL(KIND=8)::VELO_E,VELO_IB,VELO_B,VELO_PROBE
    !REAL(KIND=8)::POSI_E,POSI_IB,POSI_B,POSI_PROBE
    !
    !
    !IF( DABS(POSI_IB-POSI_B) >= DABS(POSI_B-POSI_PROBE) )THEN
    !    CALL LINEAR_INTERPOLATION(VELO_IB,VELO_PROBE,VELO_B,POSI_IB,POSI_PROBE,POSI_B)
    !ELSE
    !    CALL LINEAR_INTERPOLATION(VELO_IB,VELO_PROBE,VELO_E,POSI_IB,2.0D0*POSI_B-POSI_PROBE,POSI_E)
    !    VELO_PROBE=2.0D0*VELO_B-VELO_PROBE
    !END IF
    !
    !RETURN
    !END SUBROUTINE

    !****************************进行一个方向的速度探测******************************!
    SUBROUTINE VELOCITY_PROBING(VELO_E,VELO_I,VELO_B,VELO_S,POSI_E,POSI_I,POSI_B,POSI_S)
    IMPLICIT NONE

    REAL(KIND=8)::VELO_E,VELO_I,VELO_B,VELO_S
    REAL(KIND=8)::POSI_E,POSI_I,POSI_B,POSI_S
    REAL(KIND=8)::DB!IB点到边界的距离
    REAL(KIND=8)::DE!IB点到流域外点的距离
    REAL(KIND=8)::DS!边界到固域内点的距离
    INTEGER::PROBING_METHOD!边界到固域内点的距离

    DS=DABS(POSI_B-POSI_S)
    DB=DABS(POSI_I-POSI_B)
    DE=DABS(POSI_E-POSI_I)

    !允许部分情况下采用单点probing时等于1
    PROBING_METHOD=1

    IF( DS <= DB .AND. PROBING_METHOD==1)THEN
        VELO_S=(1.0D0+DS/DB)*VELO_B-DS/DB*VELO_I
    ELSE
        VELO_S=2.0D0*VELO_B-(DE+DB-DS)/DE*VELO_I-(DS-DB)/DE*VELO_E
    END IF

    RETURN
    END SUBROUTINE

    !*********计算K时间层(对应下一时间层的K-1时间层)，非流体内场点的拉普拉斯项*********!
    SUBROUTINE CAL_NONFLUIDIC_LAPLACE_K
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    REAL(KIND=8)::U_PROBE,V_PROBE
    REAL(KIND=8)::A1,A2,A3,B1,B2,B3

    LAPLACE_VXM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    LAPLACE_VYM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    LAPLACE_UXM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    LAPLACE_UYM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)

    DO J=1,JM-1,1
        DO I=1,IM-1,1


            IF( IABS(TYPEVX(I,J))==1 )THEN!TYPEVX(I,J)==1,-1

                A1=2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I-1) ) )
                A2=2.0D0/( ( XPV(I  )-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
                A3=2.0D0/( ( XPV(I+1)-XPV(I-1) )*( XPV(I+1)-XPV(I  ) ) )
                IF( TYPEVX(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(VM(I+1,J),VM(I,J),IB_IPSVL_VX(I,J),V_PROBE,XPV(I+1),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                    LAPLACE_VXM1(I,J)= A1*V_PROBE - A2*VM(I,J) + A3*VM(I+1,J)
                ELSE IF( TYPEVX(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I-1,J),VM(I,J),IB_IPSVL_VX(I,J),V_PROBE,XPV(I-1),XPV(I),IB_ITSCT_VX(I,J),XPV(I+1))
                    LAPLACE_VXM1(I,J)= A1*VM(I-1,J) - A2*VM(I,J) + A3*V_PROBE
                END IF

            END IF

            IF( IABS(TYPEVY(I,J))==1 )THEN!TYPEVY(I,J)==1,-1

                B1=2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J-1) ) )
                B2=2.0D0/( ( Y(J  )-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
                B3=2.0D0/( ( Y(J+1)-Y(J-1) )*( Y(J+1)-Y(J  ) ) )
                IF( TYPEVY(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(VM(I,J+1),VM(I,J),IB_IPSVL_VY(I,J),V_PROBE,Y(J+1),Y(J),IB_ITSCT_VY(I,J),Y(J-1))
                    LAPLACE_VYM1(I,J)= B1*V_PROBE - B2*VM(I,J) + B3*VM(I,J+1)
                ELSE IF( TYPEVY(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I,J-1),VM(I,J),IB_IPSVL_VY(I,J),V_PROBE,Y(J-1),Y(J),IB_ITSCT_VY(I,J),Y(J+1))
                    LAPLACE_VYM1(I,J)= B1*VM(I,J-1) - B2*VM(I,J) + B3*V_PROBE
                END IF

            END IF

            IF( IABS(TYPEUX(I,J))==1 )THEN!TYPEUX(I,J)==1,-1

                A1=2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I-1) ) )
                A2=2.0D0/( ( X(I  )-X(I-1) )*( X(I+1)-X(I  ) ) )
                A3=2.0D0/( ( X(I+1)-X(I-1) )*( X(I+1)-X(I  ) ) )
                IF( TYPEUX(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(UM(I+1,J),UM(I,J),IB_IPSVL_UX(I,J),U_PROBE,X(I+1),X(I),IB_ITSCT_UX(I,J),X(I-1))
                    LAPLACE_UXM1(I,J)= A1*U_PROBE - A2*UM(I,J) + A3*UM(I+1,J)
                ELSE IF( TYPEUX(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I-1,J),UM(I,J),IB_IPSVL_UX(I,J),U_PROBE,X(I-1),X(I),IB_ITSCT_UX(I,J),X(I+1))
                    LAPLACE_UXM1(I,J)= A1*UM(I-1,J) - A2*UM(I,J) + A3*U_PROBE
                END IF

            END IF

            IF( IABS(TYPEUY(I,J))==1 )THEN!TYPEUY(I,J)==1,-1

                B1=2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J-1) ) )
                B2=2.0D0/( ( YPU(J  )-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
                B3=2.0D0/( ( YPU(J+1)-YPU(J-1) )*( YPU(J+1)-YPU(J  ) ) )
                IF( TYPEUY(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(UM(I,J+1),UM(I,J),IB_IPSVL_UY(I,J),U_PROBE,YPU(J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                    LAPLACE_UYM1(I,J)= B1*U_PROBE - B2*UM(I,J) + B3*UM(I,J+1)
                ELSE IF( TYPEUY(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I,J-1),UM(I,J),IB_IPSVL_UY(I,J),U_PROBE,YPU(J-1),YPU(J),IB_ITSCT_UY(I,J),YPU(J+1))
                    LAPLACE_UYM1(I,J)= B1*UM(I,J-1) - B2*UM(I,J) + B3*U_PROBE
                END IF

            END IF

        END DO
    END DO



    RETURN
    END SUBROUTINE


    !*********计算K时间层(对应下一时间层的K-1时间层)，非流体内场点的对流项*********!
    SUBROUTINE CAL_NONFLUIDIC_CONVECT_K
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    USE QUADRIC_PARAMETER
    IMPLICIT NONE

    REAL(KIND=8)::U_PROBE,V_PROBE
    REAL(KIND=8)::D1,D2
    REAL(KIND=8)::URK,ULK,VRK,VLK,XRK,XLK
    REAL(KIND=8)::UTK,UBK,VTK,VBK,YTK,YBK
    REAL(KIND=8)::CVCMARK!特殊处理的标记


    CONVECT_VXM2=CONVECT_VXM1
    CONVECT_VYM2=CONVECT_VYM1
    CONVECT_UXM2=CONVECT_UXM1
    CONVECT_UYM2=CONVECT_UYM1
    !置零这样在TIMEMARCH中直接不进行CONVECT_***2的计算了
    IF(NSUBSTEP==3)THEN
        CONVECT_VXM2=0.0D0
        CONVECT_VYM2=0.0D0
        CONVECT_UXM2=0.0D0
        CONVECT_UYM2=0.0D0
    END IF

    CONVECT_VXM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    CONVECT_VYM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    CONVECT_UXM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    CONVECT_UYM1=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)


    DO J=1,JM-1,1
        DO I=1,IM-1,1

            !------V: V^2/Y------!
            IF( IABS(TYPEVY(I,J))==1 )THEN!TYPEVY(I,J)==1,-1

                D1=Y(J-1)-Y(J)
                D2=Y(J+1)-Y(J)
                IF( TYPEVY(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(VM(I,J+1),VM(I,J),IB_IPSVL_VY(I,J),V_PROBE,Y(J+1),Y(J),IB_ITSCT_VY(I,J),Y(J-1))
                    CONVECT_VYM1(I,J)= ( D2*D2*V_PROBE*V_PROBE+(D1*D1-D2*D2)*VM(I,J)*VM(I,J)-D1*D1*VM(I,J+1)*VM(I,J+1) )/( D1*D2*(D2-D1) )
                ELSE IF( TYPEVY(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I,J-1),VM(I,J),IB_IPSVL_VY(I,J),V_PROBE,Y(J-1),Y(J),IB_ITSCT_VY(I,J),Y(J+1))
                    CONVECT_VYM1(I,J)= ( D2*D2*VM(I,J-1)*VM(I,J-1)+(D1*D1-D2*D2)*VM(I,J)*VM(I,J)-D1*D1*V_PROBE*V_PROBE )/( D1*D2*(D2-D1) )
                END IF

            END IF

            !------V: UV/X------!URK,ULK,VRK,VLK,XRK,XLK
            IF(TYPEVX(I,J)/=-10 .AND. TYPEVX(I,J)/=0)THEN
                CVCMARK=0
                !URK,ULK,VRK,VLK,XRK,XLK正常数值
                CALL LINEAR_INTERPOLATION( VM(I  ,J),VLK,VM(I-1,J),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( VM(I+1,J),VRK,VM(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
                CALL LINEAR_INTERPOLATION( UM(I  ,J),ULK,UM(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( UM(I+1,J),URK,UM(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                XRK=X(I+1)
                XLK=X(I  )
                !VLK
                IF( TYPEVX(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(VM(I+1,J),VM(I,J),IB_IPSVL_VX(I,J),V_PROBE,XPV(I+1),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                    CALL LINEAR_INTERPOLATION( VM(I  ,J),VLK,V_PROBE,XPV(I  ),X(I  ),XPV(I-1) )
                    CVCMARK=1
                END IF
                !VRK
                IF( TYPEVX(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I-1,J),VM(I,J),IB_IPSVL_VX(I,J),V_PROBE,XPV(I-1),XPV(I),IB_ITSCT_VX(I,J),XPV(I+1))
                    CALL LINEAR_INTERPOLATION( V_PROBE,VRK,VM(I  ,J),XPV(I+1),X(I+1),XPV(I  ) )
                    CVCMARK=1
                END IF
                !ULK
                IF( TYPEUY(I,J)==1 .AND. TYPEUY(I,J-1)==-10 )THEN
                    CALL VELOCITY_PROBING(UM(I,J+1),UM(I,J),IB_IPSVL_UY(I,J),U_PROBE,YPU(J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                    CALL LINEAR_INTERPOLATION( UM(I  ,J),ULK,U_PROBE,YPU(J  ),Y(J  ),YPU(J-1) )
                    CVCMARK=1
                ELSE IF( TYPEUY(I,J)==-10 .AND. TYPEUY(I,J-1)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I,J-2),UM(I,J-1),IB_IPSVL_UY(I,J-1),U_PROBE,YPU(J-2),YPU(J-1),IB_ITSCT_UY(I,J-1),YPU(J))
                    CALL LINEAR_INTERPOLATION( U_PROBE,ULK,UM(I  ,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                    CVCMARK=1
                END IF
                !URK
                IF( TYPEUY(I+1,J)==1 .AND. TYPEUY(I+1,J-1)==-10 )THEN
                    CALL VELOCITY_PROBING(UM(I+1,J+1),UM(I+1,J),IB_IPSVL_UY(I+1,J),U_PROBE,YPU(J+1),YPU(J),IB_ITSCT_UY(I+1,J),YPU(J-1))
                    CALL LINEAR_INTERPOLATION( UM(I+1,J),URK,U_PROBE,YPU(J  ),Y(J  ),YPU(J-1) )
                    CVCMARK=1
                ELSE IF( TYPEUY(I+1,J)==-10 .AND. TYPEUY(I+1,J-1)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I+1,J-2),UM(I+1,J-1),IB_IPSVL_UY(I+1,J-1),U_PROBE,YPU(J-2),YPU(J-1),IB_ITSCT_UY(I+1,J-1),YPU(J))
                    CALL LINEAR_INTERPOLATION( U_PROBE,URK,UM(I+1,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                    CVCMARK=1
                END IF
                !ULK,VLK,XLK
                IF( TYPEUY(I,J)==-10 .AND. TYPEUY(I,J-1)==-10 )THEN
                    XLK=IB_ITSCT_VX(I,J)
                    ULK=IB_IPSVL_VXU(I,J)
                    VLK=IB_IPSVL_VX(I,J)
                    CVCMARK=1
                END IF
                !URK,VRK,XRK
                IF( TYPEUY(I+1,J)==-10 .AND. TYPEUY(I+1,J-1)==-10 )THEN
                    XRK=IB_ITSCT_VX(I,J)
                    URK=IB_IPSVL_VXU(I,J)
                    VRK=IB_IPSVL_VX(I,J)
                    CVCMARK=1
                END IF

                !最终计算
                IF(CVCMARK==1)THEN
                    CONVECT_VXM1(I,J)= ( URK*VRK-ULK*VLK )/( XRK-XLK )
                END IF
            END IF

            !------U: U^2/X------!
            IF( IABS(TYPEUX(I,J))==1 )THEN!TYPEUX(I,J)==1,-1

                D1=X(I-1)-X(I)
                D2=X(I+1)-X(I)
                IF( TYPEUX(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(UM(I+1,J),UM(I,J),IB_IPSVL_UX(I,J),U_PROBE,X(I+1),X(I),IB_ITSCT_UX(I,J),X(I-1))
                    CONVECT_UXM1(I,J)= ( D2*D2*U_PROBE*U_PROBE+(D1*D1-D2*D2)*UM(I,J)*UM(I,J)-D1*D1*UM(I+1,J)*UM(I+1,J) )/( D1*D2*(D2-D1) )
                ELSE IF( TYPEUX(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I-1,J),UM(I,J),IB_IPSVL_UX(I,J),U_PROBE,X(I-1),X(I),IB_ITSCT_UX(I,J),X(I+1))
                    CONVECT_UXM1(I,J)= ( D2*D2*UM(I-1,J)*UM(I-1,J)+(D1*D1-D2*D2)*UM(I,J)*UM(I,J)-D1*D1*U_PROBE*U_PROBE )/( D1*D2*(D2-D1) )
                END IF

            END IF


            !------U: UV/Y------!UTK,UBK,VTK,VBK,YTK,YBK
            IF(TYPEUY(I,J)/=-10 .AND. TYPEUY(I,J)/=0)THEN
                CVCMARK=0
                !UTK,UBK,VTK,VBK,YTK,YBK正常数值
                CALL LINEAR_INTERPOLATION( UM(I,J+1),UTK,UM(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
                CALL LINEAR_INTERPOLATION( UM(I,J  ),UBK,UM(I,J-1),YPU(J  ),Y(J  ),YPU(J-1) )
                CALL LINEAR_INTERPOLATION( VM(I,J+1),VTK,VM(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
                CALL LINEAR_INTERPOLATION( VM(I,J  ),VBK,VM(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
                YTK=Y(J+1)
                YBK=Y(J  )
                !UBK
                IF( TYPEUY(I,J)==1 )THEN
                    CALL VELOCITY_PROBING(UM(I,J+1),UM(I,J),IB_IPSVL_UY(I,J),U_PROBE,YPU(J+1),YPU(J),IB_ITSCT_UY(I,J),YPU(J-1))
                    CALL LINEAR_INTERPOLATION( UM(I,J  ),UBK,U_PROBE,YPU(J  ),Y(J  ),YPU(J-1) )
                    CVCMARK=1
                END IF
                !UTK
                IF( TYPEUY(I,J)==-1 )THEN
                    CALL VELOCITY_PROBING(UM(I,J-1),UM(I,J),IB_IPSVL_UY(I,J),U_PROBE,YPU(J-1),YPU(J),IB_ITSCT_UY(I,J),YPU(J+1))
                    CALL LINEAR_INTERPOLATION( U_PROBE,UTK,UM(I,J  ),YPU(J+1),Y(J+1),YPU(J  ) )
                    CVCMARK=1
                END IF
                !VBK
                IF( TYPEVX(I,J)==1 .AND. TYPEVX(I-1,J)==-10 )THEN
                    CALL VELOCITY_PROBING(VM(I+1,J),VM(I,J),IB_IPSVL_VX(I,J),V_PROBE,XPV(I+1),XPV(I),IB_ITSCT_VX(I,J),XPV(I-1))
                    CALL LINEAR_INTERPOLATION( VM(I,J  ),VBK,V_PROBE,XPV(I  ),X(I  ),XPV(I-1) )
                    CVCMARK=1
                ELSE IF( TYPEVX(I,J)==-10 .AND. TYPEVX(I-1,J)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I-2,J),VM(I-1,J),IB_IPSVL_VX(I-1,J),V_PROBE,XPV(I-2),XPV(I-1),IB_ITSCT_VX(I-1,J),XPV(I))
                    CALL LINEAR_INTERPOLATION( V_PROBE,VBK,VM(I-1,J  ),XPV(I  ),X(I  ),XPV(I-1) )
                    CVCMARK=1
                END IF
                !VTK
                IF( TYPEVX(I,J+1)==1 .AND. TYPEVX(I-1,J+1)==-10 )THEN
                    CALL VELOCITY_PROBING(VM(I+1,J+1),VM(I,J+1),IB_IPSVL_VX(I,J+1),V_PROBE,XPV(I+1),XPV(I),IB_ITSCT_VX(I,J+1),XPV(I-1))
                    CALL LINEAR_INTERPOLATION( VM(I,J+1),VTK,V_PROBE,XPV(I  ),X(I  ),XPV(I-1) )
                    CVCMARK=1
                ELSE IF( TYPEVX(I,J+1)==-10 .AND. TYPEVX(I-1,J+1)==-1 )THEN
                    CALL VELOCITY_PROBING(VM(I-2,J+1),VM(I-1,J+1),IB_IPSVL_VX(I-1,J+1),V_PROBE,XPV(I-2),XPV(I-1),IB_ITSCT_VX(I-1,J+1),XPV(I))
                    CALL LINEAR_INTERPOLATION( V_PROBE,VTK,VM(I-1,J+1),XPV(I  ),X(I  ),XPV(I-1) )
                    CVCMARK=1
                END IF

                !UBK,VBK,YBK
                IF( TYPEVX(I,J)==-10 .AND. TYPEVX(I-1,J)==-10 )THEN
                    YBK=IB_ITSCT_UY(I,J)
                    UBK=IB_IPSVL_UY(I,J)
                    VBK=IB_IPSVL_UYV(I,J)
                    CVCMARK=1
                END IF
                !UTK,VTK,YTK
                IF( TYPEVX(I,J+1)==-10 .AND. TYPEVX(I-1,J+1)==-10 )THEN
                    YTK=IB_ITSCT_UY(I,J)
                    UTK=IB_IPSVL_UY(I,J)
                    VTK=IB_IPSVL_UYV(I,J)
                    CVCMARK=1
                END IF

                !最终计算
                IF(CVCMARK==1)THEN
                    CONVECT_UYM1(I,J)=( UTK*VTK-UBK*VBK )/( YTK-YBK )
                END IF
            END IF


        END DO
    END DO

    RETURN
    END SUBROUTINE


    !*********计算K时间层非流体内场点的左端系数和移到右端项的大小*********!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!DELTA_VELO_BOUNDARY的取值现在严格按照静边界给定，适用于动边界的之后再想办法!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE IBM_PRIMITIVE2DERIVATIVE_TRUECARTESIAN
    USE IMMERSED_BOUNDARY
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::DB!IB点到边界的距离
    REAL(KIND=8)::DE!IB点到流域外点的距离
    REAL(KIND=8)::DS!边界到固域内点的距离
    REAL(KIND=8)::DELTA_VELO_BOUNDARY!边界点速度差值
    INTEGER::TSIGN!坐标符号

    REAL(KIND=8)::AS,AI,AE,BS,BI,BE!正常离散格式左端项系数
    REAL(KIND=8)::CSB,CSI,CSE!VELOCITY_PROBING系数

    !------V------!
    DO J=1,JM,1
        DO I=0,IM,1
            !------X------!
            IF(TYPEVX(I,J)==0 .OR. TYPEVX(I,J)==-10)THEN
                IF( NSUBSTEP==1)THEN
                    IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    IB_RVX(I,J)=0.0D0
                END IF
            END IF
            IF( IABS(TYPEVX(I,J))==1 )THEN
                !1
                TSIGN=TYPEVX(I,J)
                DB=DBLE(TSIGN)*(XPV(I)-IB_ITSCT_VX(I,J))
                DE=DBLE(TSIGN)*(XPV(I+TSIGN)-XPV(I))
                DS=DBLE(TSIGN)*(IB_ITSCT_VX(I,J)-XPV(I-TSIGN))
                !2
                AS=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/(DS+DB+DE)
                AI=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/DE
                AE=-ALPHA(NSUBSTEP)*DT/Re/DE/(DS+DB+DE)
                !3
                IF(DS<=DB)THEN
                    CSB=1.0D0+DS/DB
                    CSI=-DS/DB
                    CSE=0.0D0
                ELSE
                    CSB=2.0D0
                    CSI=-(DE+DB-DS)/DE
                    CSE=-(DS-DB)/DE
                END IF
                !4
                IF(CASE_TYPE==1 .AND. NSTEP==1 .AND. NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_VX(I,J)-V_FREESTREAM
                ELSE IF( NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_VX(I,J)-IB_IPSVL_VXM1(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    DELTA_VELO_BOUNDARY=0.0D0
                END IF
                IB_RVX(I,J)=-AS*CSB*DELTA_VELO_BOUNDARY
                IB_CVX(I,J)=-AI+AS*CSI
                IB_AVX(I,J)=AE+AS*CSE

            END IF

            !------Y------!
            IF(TYPEVY(I,J)==0 .OR. TYPEVY(I,J)==-10)THEN
                IF( NSUBSTEP==1)THEN
                    IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    IB_RVY(I,J)=0.0D0
                END IF
            END IF
            IF( IABS(TYPEVY(I,J))==1 )THEN
                !1
                TSIGN=TYPEVY(I,J)
                DB=DBLE(TSIGN)*(Y(J)-IB_ITSCT_VY(I,J))
                DE=DBLE(TSIGN)*(Y(J+TSIGN)-Y(J))
                DS=DBLE(TSIGN)*(IB_ITSCT_VY(I,J)-Y(J-TSIGN))
                !2
                AS=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/(DS+DB+DE)
                AI=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/DE
                AE=-ALPHA(NSUBSTEP)*DT/Re/DE/(DS+DB+DE)
                !3
                IF(DS<=DB)THEN
                    CSB=1.0D0+DS/DB
                    CSI=-DS/DB
                    CSE=0.0D0
                ELSE
                    CSB=2.0D0
                    CSI=-(DE+DB-DS)/DE
                    CSE=-(DS-DB)/DE
                END IF
                !4
                IF(CASE_TYPE==1 .AND. NSTEP==1 .AND. NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_VY(I,J)-V_FREESTREAM
                ELSE IF( NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_VY(I,J)-IB_IPSVL_VYM1(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    DELTA_VELO_BOUNDARY=0.0D0
                END IF
                IB_RVY(I,J)=-AS*CSB*DELTA_VELO_BOUNDARY
                IB_CVY(I,J)=-AI+AS*CSI
                IB_BVY(I,J)=AE+AS*CSE

            END IF

        END DO
    END DO

    !------U------!
    DO J=0,JM,1
        DO I=1,IM,1
            !------X------!
            IF(TYPEUX(I,J)==0 .OR. TYPEUX(I,J)==-10)THEN
                IF( NSUBSTEP==1)THEN
                    IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    IB_RUX(I,J)=0.0D0
                END IF
            END IF
            IF( IABS(TYPEUX(I,J))==1 )THEN
                !1
                TSIGN=TYPEUX(I,J)
                DB=DBLE(TSIGN)*(X(I)-IB_ITSCT_UX(I,J))
                DE=DBLE(TSIGN)*(X(I+TSIGN)-X(I))
                DS=DBLE(TSIGN)*(IB_ITSCT_UX(I,J)-X(I-TSIGN))
                !2
                AS=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/(DS+DB+DE)
                AI=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/DE
                AE=-ALPHA(NSUBSTEP)*DT/Re/DE/(DS+DB+DE)
                !3
                IF(DS<=DB)THEN
                    CSB=1.0D0+DS/DB
                    CSI=-DS/DB
                    CSE=0.0D0
                ELSE
                    CSB=2.0D0
                    CSI=-(DE+DB-DS)/DE
                    CSE=-(DS-DB)/DE
                END IF
                !4
                IF(CASE_TYPE==1 .AND. NSTEP==1 .AND. NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_UX(I,J)-U_FREESTREAM
                ELSE IF( NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_UX(I,J)-IB_IPSVL_UXM1(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    DELTA_VELO_BOUNDARY=0.0D0
                END IF
                IB_RUX(I,J)=-AS*CSB*DELTA_VELO_BOUNDARY
                IB_CUX(I,J)=-AI+AS*CSI
                IB_AUX(I,J)=AE+AS*CSE

            END IF

            !------Y------!
            IF(TYPEUY(I,J)==0 .OR. TYPEUY(I,J)==-10)THEN
                IF( NSUBSTEP==1)THEN
                    IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    IB_RUY(I,J)=0.0D0
                END IF
            END IF
            IF( IABS(TYPEUY(I,J))==1 )THEN
                !1
                TSIGN=TYPEUY(I,J)
                DB=DBLE(TSIGN)*(YPU(J)-IB_ITSCT_UY(I,J))
                DE=DBLE(TSIGN)*(YPU(J+TSIGN)-YPU(J))
                DS=DBLE(TSIGN)*(IB_ITSCT_UY(I,J)-YPU(J-TSIGN))
                !2
                AS=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/(DS+DB+DE)
                AI=-ALPHA(NSUBSTEP)*DT/Re/(DS+DB)/DE
                AE=-ALPHA(NSUBSTEP)*DT/Re/DE/(DS+DB+DE)
                !3
                IF(DS<=DB)THEN
                    CSB=1.0D0+DS/DB
                    CSI=-DS/DB
                    CSE=0.0D0
                ELSE
                    CSB=2.0D0
                    CSI=-(DE+DB-DS)/DE
                    CSE=-(DS-DB)/DE
                END IF
                !4
                IF(CASE_TYPE==1 .AND. NSTEP==1 .AND. NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_UY(I,J)-U_FREESTREAM
                ELSE IF( NSUBSTEP==1)THEN
                    DELTA_VELO_BOUNDARY=IB_IPSVL_UY(I,J)-IB_IPSVL_UYM1(I,J)
                ELSE IF( NSUBSTEP/=1)THEN
                    DELTA_VELO_BOUNDARY=0.0D0
                END IF
                IB_RUY(I,J)=-AS*CSB*DELTA_VELO_BOUNDARY
                IB_CUY(I,J)=-AI+AS*CSI
                IB_BUY(I,J)=AE+AS*CSE

            END IF

        END DO
    END DO


    RETURN
    END SUBROUTINE

