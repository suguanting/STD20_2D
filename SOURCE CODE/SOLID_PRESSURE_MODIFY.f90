    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************消除固体内部异常压力方法1********************************************!
    SUBROUTINE SOLID_PRESSURE_MODIFY_METHOD1
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE

    INTEGER     ,ALLOCATABLE::TYPEP(:,:),TYPEP_T1(:,:),TYPEP_T2(:,:)!10-FLOW
    INTEGER     ,ALLOCATABLE::I_10MIN(:),I_10MAX(:)
    INTEGER     ,ALLOCATABLE::J_10MIN(:),J_10MAX(:)
    INTEGER::IJPMAX(2),IJPMIN(2)
    INTEGER::I_PMAX,J_PMAX,I_PMIN,J_PMIN
    INTEGER::MODIFY_COUNT
    INTEGER::MODIFY_DISTANCE

    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( TYPEP(IM-1,JM-1),TYPEP_T1(IM-1,JM-1),TYPEP_T2(IM-1,JM-1) )
    ALLOCATE( I_10MIN(JM-1),I_10MAX(JM-1) )
    ALLOCATE( J_10MIN(IM-1),J_10MAX(IM-1) )

    TYPEP=10
    TYPEP_T1=10
    TYPEP_T2=10

    I_10MIN=0
    I_10MAX=0
    J_10MIN=0
    J_10MAX=0

    MODIFY_DISTANCE=IDNINT(0.05D0/DX3)

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            IF( TYPEUX(I,J)==-10 .AND. TYPEUY(I,J)==-10 .AND. TYPEUX(I+1,J)==-10 .AND. TYPEUY(I+1,J)==-10 .AND. TYPEVX(I,J)==-10 .AND. TYPEVY(I,J)==-10 .AND. TYPEVX(I,J+1)==-10 .AND. TYPEVY(I,J+1)==-10 )THEN
                TYPEP_T1(I,J)=-10
            END IF
        END DO
    END DO

    DO J=2,JM-2,1
        DO I=2,IM-2,1
            IF( TYPEP_T1(I-1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MIN(J)=I
            END IF
            IF( TYPEP_T1(I+1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MAX(J)=I
            END IF
            IF( TYPEP_T1(I,J-1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MIN(I)=J
            END IF
            IF( TYPEP_T1(I,J+1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MAX(I)=J
            END IF
        END DO
    END DO

    DO I=2,IM-2,1
        DO J=2,JM-2,1
            IF( TYPEP_T1(I,J)==-10 .AND. I-I_10MIN(J)>=MODIFY_DISTANCE .AND. I-I_10MAX(J)<=-MODIFY_DISTANCE .AND. J-J_10MIN(I)>=MODIFY_DISTANCE .AND. J-J_10MAX(I)<=-MODIFY_DISTANCE )THEN
                TYPEP_T2(I,J)=-10
            END IF
        END DO
    END DO

    TYPEP=TYPEP_T2

    IJPMAX=MAXLOC( P )
    IJPMIN=MINLOC( P )
    I_PMAX=IJPMAX( 1 )
    J_PMAX=IJPMAX( 2 )
    I_PMIN=IJPMIN( 1 )
    J_PMIN=IJPMIN( 2 )
    MODIFY_COUNT=0
    IF( TYPEP(I_PMAX,J_PMAX)==-10 .OR. TYPEP(I_PMIN,J_PMIN)==-10 )THEN
        DO I=1,IM-1,1
            DO J=1,JM-1,1
                IF( TYPEP(I,J)==-10 )THEN
                    P(I,J)=0.0D0
                    MODIFY_COUNT=MODIFY_COUNT+1
                END IF
            END DO
        END DO
        WRITE(*,*)"已修改:",MODIFY_COUNT

    END IF

    RETURN
    END SUBROUTINE

    
    !***************************************消除固体内部异常压力方法(x方向运动用)********************************************!
    SUBROUTINE SOLID_PRESSURE_MODIFY_METHOD3X
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE

    INTEGER     ,ALLOCATABLE::TYPEP(:,:),TYPEP_T1(:,:),TYPEP_T2(:,:)!10-FLOW
    INTEGER     ,ALLOCATABLE::I_10MIN(:),I_10MAX(:)
    INTEGER     ,ALLOCATABLE::J_10MIN(:),J_10MAX(:)
    INTEGER::IJPMAX(2),IJPMIN(2)
    INTEGER::I_PMAX,J_PMAX,I_PMIN,J_PMIN
    INTEGER::I_MODIFYTEMP(1),I_MODIFY
    INTEGER::MODIFY_COUNT
    INTEGER::MODIFY_DISTANCE

    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( TYPEP(IM-1,JM-1),TYPEP_T1(IM-1,JM-1),TYPEP_T2(IM-1,JM-1) )
    ALLOCATE( I_10MIN(JM-1),I_10MAX(JM-1) )
    ALLOCATE( J_10MIN(IM-1),J_10MAX(IM-1) )

    TYPEP=10
    TYPEP_T1=10
    TYPEP_T2=10

    I_10MIN=0
    I_10MAX=0
    J_10MIN=0
    J_10MAX=0

    MODIFY_DISTANCE=IDNINT(0.1D0/DX3)

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            IF( TYPEUX(I,J)==-10 .AND. TYPEUY(I,J)==-10 .AND. TYPEUX(I+1,J)==-10 .AND. TYPEUY(I+1,J)==-10 .AND. TYPEVX(I,J)==-10 .AND. TYPEVY(I,J)==-10 .AND. TYPEVX(I,J+1)==-10 .AND. TYPEVY(I,J+1)==-10 )THEN
                TYPEP_T1(I,J)=-10
            END IF
        END DO
    END DO

    DO J=2,JM-2,1
        DO I=2,IM-2,1
            IF( TYPEP_T1(I-1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MIN(J)=I
            END IF
            IF( TYPEP_T1(I+1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MAX(J)=I
            END IF
            IF( TYPEP_T1(I,J-1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MIN(I)=J
            END IF
            IF( TYPEP_T1(I,J+1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MAX(I)=J
            END IF
        END DO
    END DO

    !J_10MAX=J_10MAX-J_10MIN
    !I_MODIFYTEMP=MAXLOC(J_10MAX)
    !I_MODIFY=I_MODIFYTEMP(1)

    DO I=2,IM-2,1
        DO J=2,JM-2,1
            IF( TYPEP_T1(I,J)==-10 .AND. I-I_10MIN(J)>=MODIFY_DISTANCE .AND. I-I_10MAX(J)<=-MODIFY_DISTANCE .AND. J-J_10MIN(I)>=1 .AND. J-J_10MAX(I)<=-1 )THEN
                !IF( TYPEP_T1(I,J)==-10 .AND. (I-I_10MAX(J)<=-MODIFY_DISTANCE .OR. I<=I_MODIFY) )THEN
                TYPEP_T2(I,J)=-10
            END IF
        END DO
    END DO

    !TYPEP=TYPEP_T1
    TYPEP=TYPEP_T2

    IJPMAX=MAXLOC( P )
    IJPMIN=MINLOC( P )
    I_PMAX=IJPMAX( 1 )
    J_PMAX=IJPMAX( 2 )
    I_PMIN=IJPMIN( 1 )
    J_PMIN=IJPMIN( 2 )
    MODIFY_COUNT=0
    IF( TYPEP(I_PMAX,J_PMAX)==-10 .OR. TYPEP(I_PMIN,J_PMIN)==-10 )THEN
        DO I=1,IM-1,1
            DO J=1,JM-1,1
                IF( TYPEP_T2(I,J)==-10 )THEN
                    P(I,J)=0.0D0
                    MODIFY_COUNT=MODIFY_COUNT+1
                END IF
            END DO
        END DO
        WRITE(*,*)"已修改:",MODIFY_COUNT

    END IF

    IF(MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPLT) ) )==0 .OR. NSTEP==NSTART )THEN
        WRITE(CHAR_STEP,'(I6.6)') NSTEP
        REYNOLDS=IDNINT(Re)
        WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

        OPEN(UNIT=10,FILE='TYPEP'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
        WRITE(10,*) 'TITLE="NONAME"'
        WRITE(10,*) 'VARIABLES="X","Y","TYPEP_T1","TYPEP_T2","TYPEP"'!
        WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

        DO J=1,JM-1,1
            DO I=1,IM-1,1
                WRITE(10,*) XPV(I),YPU(J),TYPEP_T1(I,J),TYPEP_T2(I,J),TYPEP(I,J)
            END DO
        END DO

        CLOSE(10)
    END IF

    RETURN
    END SUBROUTINE
    
    !***************************************消除固体内部异常压力方法(Y方向运动用)********************************************!
    SUBROUTINE SOLID_PRESSURE_MODIFY_METHOD3Y
    USE DECLARATION
    USE IMMERSED_BOUNDARY
    IMPLICIT NONE

    INTEGER     ,ALLOCATABLE::TYPEP(:,:),TYPEP_T1(:,:),TYPEP_T2(:,:)!10-FLOW
    INTEGER     ,ALLOCATABLE::I_10MIN(:),I_10MAX(:)
    INTEGER     ,ALLOCATABLE::J_10MIN(:),J_10MAX(:)
    INTEGER::IJPMAX(2),IJPMIN(2)
    INTEGER::I_PMAX,J_PMAX,I_PMIN,J_PMIN
    INTEGER::I_MODIFYTEMP(1),I_MODIFY
    INTEGER::MODIFY_COUNT
    INTEGER::MODIFY_DISTANCE

    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER REYNOLDS

    ALLOCATE( TYPEP(IM-1,JM-1),TYPEP_T1(IM-1,JM-1),TYPEP_T2(IM-1,JM-1) )
    ALLOCATE( I_10MIN(JM-1),I_10MAX(JM-1) )
    ALLOCATE( J_10MIN(IM-1),J_10MAX(IM-1) )

    TYPEP=10
    TYPEP_T1=10
    TYPEP_T2=10

    I_10MIN=0
    I_10MAX=0
    J_10MIN=0
    J_10MAX=0

    MODIFY_DISTANCE=IDNINT(0.1D0/DX3)

    DO J=1,JM-1,1
        DO I=1,IM-1,1
            IF( TYPEUX(I,J)==-10 .AND. TYPEUY(I,J)==-10 .AND. TYPEUX(I+1,J)==-10 .AND. TYPEUY(I+1,J)==-10 .AND. TYPEVX(I,J)==-10 .AND. TYPEVY(I,J)==-10 .AND. TYPEVX(I,J+1)==-10 .AND. TYPEVY(I,J+1)==-10 )THEN
                TYPEP_T1(I,J)=-10
            END IF
        END DO
    END DO

    DO J=2,JM-2,1
        DO I=2,IM-2,1
            IF( TYPEP_T1(I-1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MIN(J)=I
            END IF
            IF( TYPEP_T1(I+1,J)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                I_10MAX(J)=I
            END IF
            IF( TYPEP_T1(I,J-1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MIN(I)=J
            END IF
            IF( TYPEP_T1(I,J+1)==10 .AND. TYPEP_T1(I,J)==-10 )THEN
                J_10MAX(I)=J
            END IF
        END DO
    END DO

    !J_10MAX=J_10MAX-J_10MIN
    !I_MODIFYTEMP=MAXLOC(J_10MAX)
    !I_MODIFY=I_MODIFYTEMP(1)

    DO I=2,IM-2,1
        DO J=2,JM-2,1
            IF( TYPEP_T1(I,J)==-10 .AND. I-I_10MIN(J)>=1 .AND. I-I_10MAX(J)<=-1 .AND. J-J_10MIN(I)>=MODIFY_DISTANCE .AND. J-J_10MAX(I)<=-MODIFY_DISTANCE )THEN
                !IF( TYPEP_T1(I,J)==-10 .AND. (I-I_10MAX(J)<=-MODIFY_DISTANCE .OR. I<=I_MODIFY) )THEN
                TYPEP_T2(I,J)=-10
            END IF
        END DO
    END DO

    !TYPEP=TYPEP_T1
    TYPEP=TYPEP_T2

    IJPMAX=MAXLOC( P )
    IJPMIN=MINLOC( P )
    I_PMAX=IJPMAX( 1 )
    J_PMAX=IJPMAX( 2 )
    I_PMIN=IJPMIN( 1 )
    J_PMIN=IJPMIN( 2 )
    MODIFY_COUNT=0
    IF( TYPEP(I_PMAX,J_PMAX)==-10 .OR. TYPEP(I_PMIN,J_PMIN)==-10 )THEN
        DO I=1,IM-1,1
            DO J=1,JM-1,1
                IF( TYPEP_T2(I,J)==-10 )THEN
                    P(I,J)=0.0D0
                    MODIFY_COUNT=MODIFY_COUNT+1
                END IF
            END DO
        END DO
        WRITE(*,*)"已修改:",MODIFY_COUNT

    END IF

    IF(MOD( NSTEP,IDNINT( DBLE(NCYCLE)/DBLE(NPLT) ) )==0 .OR. NSTEP==NSTART )THEN
        WRITE(CHAR_STEP,'(I6.6)') NSTEP
        REYNOLDS=IDNINT(Re)
        WRITE(CHAR_REYNOLDS,'(I5.5)') REYNOLDS

        OPEN(UNIT=10,FILE='TYPEP'//TRIM(CHAR_REYNOLDS)//'N'//TRIM(CHAR_STEP)//'.PLT')!,FORM='UNFORMATTED'
        WRITE(10,*) 'TITLE="NONAME"'
        WRITE(10,*) 'VARIABLES="X","Y","TYPEP_T1","TYPEP_T2","TYPEP"'!
        WRITE(10,*) 'ZONE T="NONAME", I=',IM-1,', J=',JM-1,', F=POINT'

        DO J=1,JM-1,1
            DO I=1,IM-1,1
                WRITE(10,*) XPV(I),YPU(J),TYPEP_T1(I,J),TYPEP_T2(I,J),TYPEP(I,J)
            END DO
        END DO

        CLOSE(10)
    END IF

    RETURN
    END SUBROUTINE