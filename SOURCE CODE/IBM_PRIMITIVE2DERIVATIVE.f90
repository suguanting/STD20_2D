    SUBROUTINE IBM_PRIMITIVE2DERIVATIVE
    USE IMMERSED_BOUNDARY
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::DB,DE

    !------V------!
    DO J=1,JM,1
        DO I=0,IM,1
            !------X------!
            IF(TYPEVX(I,J)==0)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            ELSE IF(TYPEVX(I,J)==1)THEN
                DE=XPV(I+1)-XPV(I)
                DB=XPV(I)-IB_ITSCT_VX(I,J)
                IB_CVX(I,J)= DT/Re/(DB*DE)
                IB_AVX(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VX(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VX(I,J)-IBN1_IPSVL_VXN(I,J))
                END IF
            ELSE IF(TYPEVX(I,J)==-1)THEN
                DE=-XPV(I-1)+XPV(I)
                DB=-XPV(I)+IB_ITSCT_VX(I,J)
                IB_CVX(I,J)= DT/Re/(DB*DE)
                IB_AVX(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VX(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VX(I,J)-IBN1_IPSVL_VXN(I,J))
                END IF
            ELSE IF(TYPEVX(I,J)==-10)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            END IF

            !------Y------!
            IF(TYPEVY(I,J)==0)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            ELSE IF(TYPEVY(I,J)==1)THEN
                DE=Y(J+1)-Y(J)
                DB=Y(J)-IB_ITSCT_VY(I,J)
                IB_CVY(I,J)= DT/Re/(DB*DE)
                IB_BVY(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VY(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VY(I,J)-IBN1_IPSVL_VYN(I,J))
                END IF
            ELSE IF(TYPEVY(I,J)==-1)THEN
                DE=-Y(J-1)+Y(J)
                DB=-Y(J)+IB_ITSCT_VY(I,J)
                IB_CVY(I,J)= DT/Re/(DB*DE)
                IB_BVY(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VY(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_VY(I,J)-IBN1_IPSVL_VYN(I,J))
                END IF
            ELSE IF(TYPEVY(I,J)==-10)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            END IF

        END DO
    END DO

    !------U------!
    DO J=0,JM,1
        DO I=1,IM,1
            !------X------!
            IF(TYPEUX(I,J)==0)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            ELSE IF(TYPEUX(I,J)==1)THEN
                DE=X(I+1)-X(I)
                DB=X(I)-IB_ITSCT_UX(I,J)
                IB_CUX(I,J)= DT/Re/(DB*DE)
                IB_AUX(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UX(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UX(I,J)-IBN1_IPSVL_UXN(I,J))
                END IF
            ELSE IF(TYPEUX(I,J)==-1)THEN
                DE=-X(I-1)+X(I)
                DB=-X(I)+IB_ITSCT_UX(I,J)
                IB_CUX(I,J)= DT/Re/(DB*DE)
                IB_AUX(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UX(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUX(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UX(I,J)-IBN1_IPSVL_UXN(I,J))
                END IF
            ELSE IF(TYPEUX(I,J)==-10)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            END IF

            !------Y------!
            IF(TYPEUY(I,J)==0)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            ELSE IF(TYPEUY(I,J)==1)THEN
                DE=YPU(J+1)-YPU(J)
                DB=YPU(J)-IB_ITSCT_UY(I,J)
                IB_CUY(I,J)= DT/Re/(DB*DE)
                IB_BUY(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UY(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UY(I,J)-IBN1_IPSVL_UYN(I,J))
                END IF
            ELSE IF(TYPEUY(I,J)==-1)THEN
                DE=-YPU(J-1)+YPU(J)
                DB=-YPU(J)+IB_ITSCT_UY(I,J)
                IB_CUY(I,J)= DT/Re/(DB*DE)
                IB_BUY(I,J)=-DT/Re/(DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UY(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUY(I,J)=DT/Re/(DB*(DB+DE))*(IB_IPSVL_UY(I,J)-IBN1_IPSVL_UYN(I,J))
                END IF
            ELSE IF(TYPEUY(I,J)==-10)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            END IF

        END DO
    END DO


    RETURN
    END SUBROUTINE


    SUBROUTINE IBM_PRIMITIVE2DERIVATIVE_NO_PECULARITY_POINT
    USE IMMERSED_BOUNDARY
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::DB,DE

    !------V------!
    DO J=1,JM,1
        DO I=0,IM,1
            !------X------!
            IF(TYPEVX(I,J)==0)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            ELSE IF(TYPEVX(I,J)==1)THEN
                DE=XPV(I+1)-XPV(I)
                DB=XPV(I)-IB_ITSCT_VX(I,J)
                IB_CVX(I,J)= DT/Re/(DE*DE)
                IB_AVX(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VX(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VX(I,J)-IBN1_IPSVL_VXN(I,J))
                END IF
            ELSE IF(TYPEVX(I,J)==-1)THEN
                DE=-XPV(I-1)+XPV(I)
                DB=-XPV(I)+IB_ITSCT_VX(I,J)
                IB_CVX(I,J)= DT/Re/(DE*DE)
                IB_AVX(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VX(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VX(I,J)-IBN1_IPSVL_VXN(I,J))
                END IF
            ELSE IF(TYPEVX(I,J)==-10)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            END IF

            !------Y------!
            IF(TYPEVY(I,J)==0)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            ELSE IF(TYPEVY(I,J)==1)THEN
                DE=Y(J+1)-Y(J)
                DB=Y(J)-IB_ITSCT_VY(I,J)
                IB_CVY(I,J)= DT/Re/(DE*DE)
                IB_BVY(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VY(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VY(I,J)-IBN1_IPSVL_VYN(I,J))
                END IF
            ELSE IF(TYPEVY(I,J)==-1)THEN
                DE=-Y(J-1)+Y(J)
                DB=-Y(J)+IB_ITSCT_VY(I,J)
                IB_CVY(I,J)= DT/Re/(DE*DE)
                IB_BVY(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VY(I,J)-V_FREESTREAM)
                ELSE
                    IB_RVY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_VY(I,J)-IBN1_IPSVL_VYN(I,J))
                END IF
            ELSE IF(TYPEVY(I,J)==-10)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            END IF

        END DO
    END DO

    !------U------!
    DO J=0,JM,1
        DO I=1,IM,1
            !------X------!
            IF(TYPEUX(I,J)==0)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            ELSE IF(TYPEUX(I,J)==1)THEN
                DE=X(I+1)-X(I)
                DB=X(I)-IB_ITSCT_UX(I,J)
                IB_CUX(I,J)= DT/Re/(DE*DE)
                IB_AUX(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UX(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UX(I,J)-IBN1_IPSVL_UXN(I,J))
                END IF
            ELSE IF(TYPEUX(I,J)==-1)THEN
                DE=-X(I-1)+X(I)
                DB=-X(I)+IB_ITSCT_UX(I,J)
                IB_CUX(I,J)= DT/Re/(DE*DE)
                IB_AUX(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UX(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUX(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UX(I,J)-IBN1_IPSVL_UXN(I,J))
                END IF
            ELSE IF(TYPEUX(I,J)==-10)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            END IF

            !------Y------!
            IF(TYPEUY(I,J)==0)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            ELSE IF(TYPEUY(I,J)==1)THEN
                DE=YPU(J+1)-YPU(J)
                DB=YPU(J)-IB_ITSCT_UY(I,J)
                IB_CUY(I,J)= DT/Re/(DE*DE)
                IB_BUY(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UY(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UY(I,J)-IBN1_IPSVL_UYN(I,J))
                END IF
            ELSE IF(TYPEUY(I,J)==-1)THEN
                DE=-YPU(J-1)+YPU(J)
                DB=-YPU(J)+IB_ITSCT_UY(I,J)
                IB_CUY(I,J)= DT/Re/(DE*DE)
                IB_BUY(I,J)=-DT*DB/Re/(DE*DE*(DB+DE))
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UY(I,J)-U_FREESTREAM)
                ELSE
                    IB_RUY(I,J)=DT/Re/(DE*(DB+DE))*(IB_IPSVL_UY(I,J)-IBN1_IPSVL_UYN(I,J))
                END IF
            ELSE IF(TYPEUY(I,J)==-10)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            END IF

        END DO
    END DO


    RETURN
    END SUBROUTINE


    SUBROUTINE IBM_PRIMITIVE2DERIVATIVE_SECOND_ORDER
    USE IMMERSED_BOUNDARY
    USE DECLARATION
    IMPLICIT NONE

    REAL(KIND=8)::DB,DE,DBI
    REAL(KIND=8)::CA,CB,CC
    REAL(KIND=8)::DENOM

    !------V------!
    DO J=1,JM,1
        DO I=0,IM,1
            !------X------!
            IF(TYPEVX(I,J)==0)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            ELSE IF(TYPEVX(I,J)==1)THEN
                DE=XPV(I+1)-XPV(I)
                DB=IB_ITSCT_VX(I,J)-XPV(I)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CVX(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_AVX(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=0.0D0
                ELSE
                    IB_RVX(I,J)=0.0D0
                END IF
            ELSE IF(TYPEVX(I,J)==-1)THEN
                DE=XPV(I-1)-XPV(I)
                DB=IB_ITSCT_VX(I,J)-XPV(I)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CVX(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_AVX(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVX(I,J)=0.0D0
                ELSE
                    IB_RVX(I,J)=0.0D0
                END IF
            ELSE IF(TYPEVX(I,J)==-10)THEN
                IB_RVX(I,J)=IB_IPSVL_VX(I,J)-VN(I,J)
            END IF

            !------Y------!
            IF(TYPEVY(I,J)==0)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            ELSE IF(TYPEVY(I,J)==1)THEN
                DE=Y(J+1)-Y(J)
                DB=IB_ITSCT_VY(I,J)-Y(J)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CVY(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_BVY(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=0.0D0
                ELSE
                    IB_RVY(I,J)=0.0D0
                END IF
            ELSE IF(TYPEVY(I,J)==-1)THEN
                DE=Y(J-1)-Y(J)
                DB=IB_ITSCT_VY(I,J)-Y(J)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CVY(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_BVY(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RVY(I,J)=0.0D0
                ELSE
                    IB_RVY(I,J)=0.0D0
                END IF
            ELSE IF(TYPEVY(I,J)==-10)THEN
                IB_RVY(I,J)=IB_IPSVL_VY(I,J)-VN(I,J)
            END IF

        END DO
    END DO

    !------U------!
    DO J=0,JM,1
        DO I=1,IM,1
            !------X------!
            IF(TYPEUX(I,J)==0)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            ELSE IF(TYPEUX(I,J)==1)THEN
                DE=X(I+1)-X(I)
                DB=IB_ITSCT_UX(I,J)-X(I)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CUX(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_AUX(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=0.0D0
                ELSE
                    IB_RUX(I,J)=0.0D0
                END IF
            ELSE IF(TYPEUX(I,J)==-1)THEN
                DE=X(I-1)-X(I)
                DB=IB_ITSCT_UX(I,J)-X(I)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CUX(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_AUX(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUX(I,J)=0.0D0
                ELSE
                    IB_RUX(I,J)=0.0D0
                END IF
            ELSE IF(TYPEUX(I,J)==-10)THEN
                IB_RUX(I,J)=IB_IPSVL_UX(I,J)-UN(I,J)
            END IF

            !------Y------!
            IF(TYPEUY(I,J)==0)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            ELSE IF(TYPEUY(I,J)==1)THEN
                DE=YPU(J+1)-YPU(J)
                DB=IB_ITSCT_UY(I,J)-YPU(J)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CUY(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_BUY(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=0.0D0
                ELSE
                    IB_RUY(I,J)=0.0D0
                END IF
            ELSE IF(TYPEUY(I,J)==-1)THEN
                DE=YPU(J-1)-YPU(J)
                DB=IB_ITSCT_UY(I,J)-YPU(J)
                DBI=DB-DE
                CA=(DBI*DBI-DB*DB)/DE
                CB=(DE*DE-DBI*DBI)/DB
                CC=(DB*DB-DE*DE)/DBI
                DENOM=DE*DE*CA+DB*DB*CB+DBI*DBI*CC
                IB_CUY(I,J)= DT/Re*( CA+CB+CC )/DENOM
                IB_BUY(I,J)=-DT/Re*CA/DENOM
                IF(CASE_TYPE==1 .AND. NSTEP==1)THEN
                    IB_RUY(I,J)=0.0D0
                ELSE
                    IB_RUY(I,J)=0.0D0
                END IF
            ELSE IF(TYPEUY(I,J)==-10)THEN
                IB_RUY(I,J)=IB_IPSVL_UY(I,J)-UN(I,J)
            END IF

        END DO
    END DO


    RETURN
    END SUBROUTINE