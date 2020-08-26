    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !******************************************从文件中读取二维流场(全功能版)*****************************************************!
    SUBROUTINE READFROMFILE_2D_STAGGERED(FILENAME)
    USE DECLARATION
    IMPLICIT NONE
    CHARACTER(LEN=80)::FILENAME
    REAL(KIND=8)::XREAD,YREAD,OMEGAR,PHIR,DIVHR
    LOGICAL ALIVE
    INTEGER::ERROR
    !原有流场直接读取量
    REAL(KIND=8),ALLOCATABLE::DIVR(:,:),UPR(:,:),VPR(:,:),XPVR(:),YPUR(:),PR(:,:)
    !原有流场计算量
    REAL(KIND=8),ALLOCATABLE::DIVC(:,:),UR(:,:),VR(:,:),XR(:),YR(:)

    REAL(KIND=8)::ERRORDIV
    REAL(KIND=8)::MODIFY1,MODIFY2

    INTEGER::IMR,JMR
    INTEGER::RESTART_TYPE
    REAL(KIND=8)::LEFTR,RIGHR,BOTTR,TOPPR
    REAL(KIND=8),ALLOCATABLE::DISTANCE_XU(:),DISTANCE_YV(:),DISTANCE_XPV(:),DISTANCE_YPU(:)
    INTEGER,ALLOCATABLE::INDEX_XU(:),INDEX_YV(:),INDEX_XPV(:),INDEX_YPU(:)
    REAL(KIND=8)::DX,DY,VX,VY


    INQUIRE(FILE=FILENAME, EXIST=ALIVE)
    IF(.NOT. ALIVE)THEN
        WRITE(*,*)TRIM(FILENAME)," DOESN'T EXIST"
        STOP
    END IF

    WRITE(*,*)"输入数字确定续算情况：1.网格设置相同 2.计算域相同，网格密度不同 3.计算域和网格密度均不同,计算域不大于原有设置"
    READ(*,*)RESTART_TYPE

    IF(RESTART_TYPE==1)THEN
        !------直接读取------!
        ALLOCATE( DIVC(IM-1,JM-1),DIVR(IM-1,JM-1),UPR(IM-1,0:JM-1),VPR(0:IM-1,JM-1) )

        !读取报错该怎么写没想好
        OPEN(50, FILE=FILENAME)
        !DO WHILE(.TRUE.)
        READ(50,*)
        READ(50,*)
        READ(50,*)
        DO J=1,JM-1,1
            DO I=1,IM-1,1
                READ(50,*,IOSTAT=ERROR)XREAD,YREAD,UPR(I,J),VPR(I,J),P(I,J),PHI(I,J),OMEGAR,DIVR(I,J),DIVHR
                IF(ERROR/=0)EXIT
                !WRITE(*,*)XREAD,YREAD,ZREAD,U(I,J,1),V(I,J,1),W(I,J,1),P(I,J,1),OMEGAXR,OMEGAYR,OMEGAZR,OMEGAR
            END DO
        END DO
        !IF(ERROR/=0)EXIT
        !END DO
        READ(50,*)
        DO J=1,JM-1,1
            READ(50,*,IOSTAT=ERROR) U(1,J)
            IF(ERROR/=0)EXIT
        END DO
        READ(50,*)
        DO I=1,IM-1,1
            READ(50,*,IOSTAT=ERROR) V(I,1)
            IF(ERROR/=0)EXIT
        END DO  

        !------U------!
        DO J=1,JM-1,1 
            DO I=2,IM,1
                U(I,J)=2.0D0*UPR(I-1,J)-U(I-1,J)
            END DO
        END DO
        U(:       ,0       )=BCU_AB*U(:       ,1       )+BCU_BB
        U(:       ,JM      )=BCU_AT*U(:       ,JM-1    )+BCU_BT
        !------V------!
        DO I=1,IM-1,1
            DO J=2,JM,1
                V(I,J)=2.0D0*VPR(I,J-1)-V(I,J-1)
            END DO
        END DO
        V(0       ,:       )=BCV_AL*V(0       ,:       )+BCV_BL
        V(IM      ,:       )=BCV_AR*V(IM-1    ,:       )+BCV_BR

        !------DIV------!
        DO J=1,JM-1,1
            DO I=1,IM-1,1
                DIVC(I,J)=( U(I+1,J)-U(I,J) )/( X(I+1)-X(I) )+( V(I,J+1)-V(I,J) )/( Y(J+1)-Y(J) )
            END DO
        END DO
        ERRORDIV=MAXVAL(DABS(DIVC-DIVR))
        WRITE(*,*) ERRORDIV,MAXLOC(DABS(DIVC-DIVR))

    ELSE
        !-----------读取数据-----------!
        WRITE(*,*)"输入原算例IM-1"
        READ(*,*)IMR
        WRITE(*,*)"输入原算例JM-1"
        READ(*,*)JMR

        ALLOCATE( DIVC(IMR,JMR),DIVR(IMR,JMR),UPR(IMR,JMR),VPR(IMR,JMR),PR(0:IMR+1,0:JMR+1) )
        ALLOCATE( XR(IMR+1),YR(JMR+1),XPVR(0:IMR+1),YPUR(0:JMR+1),UR(IMR+1,0:JMR+1),VR(0:IMR+1,JMR+1) )

        IF(RESTART_TYPE==2)THEN
            LEFTR=LEFT
            RIGHR=RIGH
            BOTTR=BOTT
            TOPPR=TOPP
        ELSE IF(RESTART_TYPE==3)THEN
            WRITE(*,*)"输入原算例LEFT"
            READ(*,*)LEFTR
            WRITE(*,*)"输入原算例RIGH"
            READ(*,*)RIGHR
            WRITE(*,*)"输入原算例BOTT"
            READ(*,*)BOTTR
            WRITE(*,*)"输入原算例TOPP"
            READ(*,*)TOPPR
        END IF

        !读取报错该怎么写没想好
        OPEN(50, FILE=FILENAME)
        !DO WHILE(.TRUE.)
        READ(50,*)
        READ(50,*)
        READ(50,*)
        DO J=1,JMR,1
            DO I=1,IMR,1
                READ(50,*,IOSTAT=ERROR)XPVR(I),YPUR(J),UPR(I,J),VPR(I,J),PR(I,J),PHIR,OMEGAR,DIVR(I,J),DIVHR
                IF(ERROR/=0)EXIT
                !WRITE(*,*)XREAD,YREAD,ZREAD,UR(I,J,1),VR(I,J,1),W(I,J,1),P(I,J,1),OMEGAXR,OMEGAYR,OMEGAZR,OMEGAR
            END DO
        END DO
        !IF(ERROR/=0)EXIT
        !END DO
        READ(50,*)
        DO J=1,JMR,1
            READ(50,*,IOSTAT=ERROR) UR(1,J)
            IF(ERROR/=0)EXIT
        END DO
        READ(50,*)
        DO I=1,JMR,1
            READ(50,*,IOSTAT=ERROR) VR(I,1)
            IF(ERROR/=0)EXIT
        END DO  
        !------补充XPVR,YPUR------!
        XPVR(0    )=2.0D0*LEFTR-XPVR(1)
        XPVR(IMR+1)=2.0D0*RIGHR-XPVR(IMR)
        YPUR(0    )=2.0D0*BOTTR-YPUR(1)
        YPUR(JMR+1)=2.0D0*TOPPR-YPUR(JMR)
        !------补充PR------!
        PR(0       ,1:JMR:1)=PR(1       ,1:JMR:1)
        PR(IMR+1   ,1:JMR:1)=0.0D0
        PR(:       ,0      )=PR(:       ,1      )
        PR(:       ,JMR+1  )=PR(:       ,JMR    )

        !------UR------!
        DO J=1,JMR,1 
            DO I=2,IMR+1,1
                UR(I,J)=2.0D0*UPR(I-1,J)-UR(I-1,J)
            END DO
        END DO
        UR(:       ,0       )=UR(:       ,1       )
        UR(:       ,JMR+1   )=UR(:       ,JMR     )
        !------VR------!
        DO I=1,IMR,1 
            DO J=2,JMR+1,1
                VR(I,J)=2.0D0*VPR(I,J-1)-VR(I,J-1)
            END DO
        END DO
        VR(0       ,:       )=V_FREESTREAM
        VR(IMR+1   ,:       )=VR(IMR     ,:      )

        !------XR------!
        XR(1)=LEFTR
        DO I=2,IMR,1
            XR(I)=2.0D0*XPVR(I-1)-XR(I-1)
        END DO
        XR(IMR+1)=RIGHR

        !------YR------!
        YR(1)=BOTTR
        DO J=2,JMR,1
            YR(J)=2.0D0*YPUR(J-1)-YR(J-1)
        END DO
        YR(JMR+1)=TOPPR

        !------DIV------!
        DO J=1,JMR,1
            DO I=1,IMR,1
                DIVC(I,J)=( UR(I+1,J)-UR(I,J) )/( XR(I+1)-XR(I) )+( VR(I,J+1)-VR(I,J) )/( YR(J+1)-YR(J) )
            END DO
        END DO
        ERRORDIV=MAXVAL(DABS(DIVC-DIVR))
        WRITE(*,*) ERRORDIV,MAXLOC(DABS(DIVC-DIVR))

        !-----------数据插值----------!
        ALLOCATE( DISTANCE_XU(IMR+1),DISTANCE_YV(JMR+1),DISTANCE_XPV(0:IMR+1),DISTANCE_YPU(0:JMR+1) )
        !ALLOCATE( INDEX_XU(IM),INDEX_YV(JM),INDEX_XPV(0:IM),INDEX_YPU(0:JM) )
        ALLOCATE( INDEX_XU(2:IM-1),INDEX_YV(2:JM-1),INDEX_XPV(1:IM-1),INDEX_YPU(1:JM-1) )

        !------INDEX------!
        DO I=2,IM-1,1
            DISTANCE_XU=XR(:)-X(I)
            INDEX_XU(I)=MINLOC( DABS(DISTANCE_XU),1 )+LBOUND(DISTANCE_XU,1)-1
            IF( XR(INDEX_XU(I))>X(I) ) INDEX_XU(I)=INDEX_XU(I)-1
        END DO
        DO I=1,IM-1,1
            DISTANCE_XPV=XPVR(:)-XPV(I)
            INDEX_XPV(I)=MINLOC( DABS(DISTANCE_XPV),1 )+LBOUND(DISTANCE_XPV,1)-1
            IF( XPVR(INDEX_XPV(I))>XPV(I) ) INDEX_XPV(I)=INDEX_XPV(I)-1
        END DO
        DO J=2,JM-1,1
            DISTANCE_YV=YR(:)-Y(J)
            INDEX_YV(J)=MINLOC( DABS(DISTANCE_YV),1 )+LBOUND(DISTANCE_YV,1)-1
            IF( YR(INDEX_YV(J))>Y(J) ) INDEX_YV(J)=INDEX_YV(J)-1
        END DO
        DO J=1,JM-1,1
            DISTANCE_YPU=YPUR(:)-YPU(J)
            INDEX_YPU(J)=MINLOC( DABS(DISTANCE_YPU),1 )+LBOUND(DISTANCE_YPU,1)-1
            IF( YPUR(INDEX_YPU(J))>YPU(J) ) INDEX_YPU(J)=INDEX_YPU(J)-1
        END DO
        !------V------!
        DO J=2,JM-1,1
            DO I=1,IM-1,1
                DX=XPVR(INDEX_XPV(I)+1)-XPVR(INDEX_XPV(I))
                DY=YR  (INDEX_YV (J)+1)-YR  (INDEX_YV (J))

                VX=XPV(I)-XPVR(INDEX_XPV(I))
                VY=Y  (J)-YR  (INDEX_YV (J))

                V(I,J)=&
                VR(INDEX_XPV(I)  ,INDEX_YV(J)  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
                VR(INDEX_XPV(I)+1,INDEX_YV(J)  )*       VX/DX *(1.0D0-VY/DY)+&
                VR(INDEX_XPV(I)  ,INDEX_YV(J)+1)*(1.0D0-VX/DX)*       VY/DY +&
                VR(INDEX_XPV(I)+1,INDEX_YV(J)+1)*       VX/DX *       VY/DY  
            END DO
        END DO
        DO J=2,JM-1,1
            V(0       ,J)=2.0D0*V_FREESTREAM-V(1       ,J)
        END DO
        V(IM      ,2:JM-1:1)=V(IM-1    ,1:JM-1:1)
        V(:       ,1       )=V_FREESTREAM!V(:       ,2       )
        V(:       ,JM      )=V_FREESTREAM!V(:       ,JM-1    )
        !------U------!
        DO J=1,JM-1,1
            DO I=2,IM-1,1
                DX=XR  (INDEX_XU (I)+1)-XR  (INDEX_XU (I))
                DY=YPUR(INDEX_YPU(J)+1)-YPUR(INDEX_YPU(J))
                VX=X  (I)-XR  (INDEX_XU (I))
                VY=YPU(J)-YPUR(INDEX_YPU(J))
                U(I,J)=&
                UR(INDEX_XU(I)  ,INDEX_YPU(J)  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
                UR(INDEX_XU(I)+1,INDEX_YPU(J)  )*       VX/DX *(1.0D0-VY/DY)+&
                UR(INDEX_XU(I)  ,INDEX_YPU(J)+1)*(1.0D0-VX/DX)*       VY/DY +&
                UR(INDEX_XU(I)+1,INDEX_YPU(J)+1)*       VX/DX *       VY/DY
            END DO
        END DO
        U(1       ,1:JM-1:1)=U_FREESTREAM
        DO J=1,JM-1,1
            U(IM      ,J)=U(IM-1    ,J)-( X(IM)-X(IM-1) )*( V(IM-1,J+1)-V(IM-1,J) )/( Y(J+1)-Y(J) )
        END DO
        U(:       ,0       )=U(:       ,1       )
        U(:       ,JM      )=U(:       ,JM-1    )
        !------P------!
        DO J=1,JM-1,1
            DO I=1,IM-1,1
                DX=XPVR(INDEX_XPV(I)+1)-XPVR(INDEX_XPV(I))
                DY=YPUR(INDEX_YPU(J)+1)-YPUR(INDEX_YPU(J))

                VX=XPV(I)-XPVR(INDEX_XPV(I))
                VY=YPU(J)-YPUR(INDEX_YPU(J))

                P(I,J)=&
                PR(INDEX_XPV(I)  ,INDEX_YPU(J)  )*(1.0D0-VX/DX)*(1.0D0-VY/DY)+&
                PR(INDEX_XPV(I)+1,INDEX_YPU(J)  )*       VX/DX *(1.0D0-VY/DY)+&
                PR(INDEX_XPV(I)  ,INDEX_YPU(J)+1)*(1.0D0-VX/DX)*       VY/DY +&
                PR(INDEX_XPV(I)+1,INDEX_YPU(J)+1)*       VX/DX *       VY/DY  
            END DO
        END DO




    END IF

    CLOSE(50)

    CALL OUTPUT_PLT_1_STAGGERED

    RETURN
    END