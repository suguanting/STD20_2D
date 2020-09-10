    !######################################################################!
    !#                                                                    #!
    !#                              �����Ӻ���                            #!
    !#                                                                    #!
    !######################################################################!

    !***********************************************�ж��Ƿ�������ͬʱ�����Ҫ������Ϣ************************************************************!
    SUBROUTINE CONVERGENCE_JUDGE
    USE DECLARATION
    IMPLICIT NONE

    IF( NSTEP/=NSTART .AND. (ERRORVELOMAX<=CRITERIA .OR. VELOMAX>=1.0D-4/CRITERIA) )THEN
        CONVERGENCE=1

        CALL OUTPUT_PLT_1_STAGGERED
        !CALL OUTPUT_FULL_STAGGERED

    END IF
    
    !IF( IB_LOCOMOTION==8 .AND. ERRORVELOMAX<=CRITERIA )THEN
    !        CONVERGENCE=0
    !
    !END IF

    RETURN
    END SUBROUTINE
