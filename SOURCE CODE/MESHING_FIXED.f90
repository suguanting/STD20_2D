    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !***************************************************划分固定网格点******************************************************!
    SUBROUTINE MESHING_FIXED
    USE DECLARATION
    IMPLICIT NONE

    !绘制一级网格
    IF(CONTINUOUS_MESH==0)THEN
        DO I=1,IM,1
            IF(I<=ILM1)THEN
                X(I)=LEFT+DBLE(I-1   )*DX1
            ELSE IF(I<=ILM2)THEN
                X(I)=LEM1+DBLE(I-ILM1)*DX21
            ELSE IF(I<=ILM3)THEN
                X(I)=LEM2+DBLE(I-ILM2)*DX22
            ELSE IF(I<=IL)THEN
                X(I)=LEM3+DBLE(I-ILM3)*DX23
            ELSE IF(I<=IR)THEN
                X(I)=LEIN+DBLE(I-IL  )*DX3
            ELSE IF(I<=IRM3)THEN
                X(I)=RIIN+DBLE(I-IR  )*DX23
            ELSE IF(I<=IRM2)THEN
                X(I)=RIM3+DBLE(I-IRM3)*DX22
            ELSE IF(I<=IRM1)THEN
                X(I)=RIM2+DBLE(I-IRM2)*DX21
            ELSE
                X(I)=RIM1+DBLE(I-IRM1)*DX1
            END IF
        END DO

        DO J=1,JM,1
            IF(J<=JBM1)THEN
                Y(J)=BOTT+DBLE(J-1   )*DX1
            ELSE IF(J<=JBM2)THEN
                Y(J)=BOM1+DBLE(J-JBM1)*DX21
            ELSE IF(J<=JBM3)THEN
                Y(J)=BOM2+DBLE(J-JBM2)*DX22
            ELSE IF(J<=JB)THEN
                Y(J)=BOM3+DBLE(J-JBM3)*DX23
            ELSE IF(J<=JT)THEN
                Y(J)=BOIN+DBLE(J-JB  )*DX3
            ELSE IF(J<=JTM3)THEN
                Y(J)=TOIN+DBLE(J-JT  )*DX23
            ELSE IF(J<=JTM2)THEN
                Y(J)=TOM3+DBLE(J-JTM3)*DX22
            ELSE IF(J<=JTM1)THEN
                Y(J)=TOM2+DBLE(J-JTM2)*DX21
            ELSE
                Y(J)=TOM1+DBLE(J-JTM1)*DX1
            END IF
        END DO
    ELSE IF(CONTINUOUS_MESH==1)THEN
        DO I=1,IM,1
            IF(I<IL)THEN
                X(I)=LEIN+HL*( BL+1.0D0-(BL-1.0D0)*( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(IL-I)/DBLE(NODEL-1) ) )/( ( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(IL-I)/DBLE(NODEL-1) )+1.0D0 )
            ELSE IF(I<=IR)THEN
                X(I)=LEIN+DBLE(I-IL  )*DX3
            ELSE
                X(I)=RIIN+HR*( BR+1.0D0-(BR-1.0D0)*( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(I-IR)/DBLE(NODER-1) ) )/( ( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(I-IR)/DBLE(NODER-1) )+1.0D0 )
            END IF
        END DO

        DO J=1,JM,1
            IF(J<JB)THEN
                Y(J)=BOIN+HB*( BB+1.0D0-(BB-1.0D0)*( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(JB-J)/DBLE(NODEB-1) ) )/( ( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(JB-J)/DBLE(NODEB-1) )+1.0D0 )
            ELSE IF(J<=JT)THEN
                Y(J)=BOIN+DBLE(J-JB  )*DX3
            ELSE
                Y(J)=TOIN+HT*( BT+1.0D0-(BT-1.0D0)*( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(J-JT)/DBLE(NODET-1) ) )/( ( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(J-JT)/DBLE(NODET-1) )+1.0D0 )
            END IF
        END DO
        !
        !WRITE(*,*) X(2,1)-X(1,1)
        !WRITE(*,*) X(IM,1)-X(IM-1,1)
        !WRITE(*,*) X(IL,1)-X(IL-1,1)
        !WRITE(*,*) X(IR+1,1)-X(IR,1)
        !
        !WRITE(*,*) Y(1,2)-Y(1,1)
        !WRITE(*,*) Y(1,JB)-Y(1,JB-1)
        !WRITE(*,*) Y(1,JT+1)-Y(1,JT)
        !WRITE(*,*) Y(1,JM)-Y(1,JM-1)
    END IF
    
    WRITE(*,*)'网格绘制误差情况'
    WRITE(*,"('X(IM)-RIGH:'(1X,F21.18))") X(IM)-RIGH
    WRITE(*,"('X(1 )-LEFT:'(1X,F21.18))") X(1)-LEFT
    WRITE(*,"('X(IR)-RIIN:'(1X,F21.18))") X(IR)-RIIN
    WRITE(*,"('Y(JM)-TOPP:'(1X,F21.18))") Y(JM)-TOPP
    WRITE(*,"('Y(1 )-BOTT:'(1X,F21.18))") Y(1)-BOTT
    WRITE(*,"('Y(JT)-TOIN:'(1X,F21.18))") Y(JT)-TOIN

    !绘制二级网格
    DO I=1,IM-1,1
        XPV(I)=( X(I)+X(I+1) )/2.0D0
    END DO
    XPV(0)=X(1)-( X(2)-X(1) )/2.0D0
    XPV(IM)=X(IM)+( X(IM)-X(IM-1) )/2.0D0

    DO J=1,JM-1,1
        YPU(J)=( Y(J)+Y(J+1) )/2.0D0
    END DO
    YPU(0)=Y(1)-( Y(2)-Y(1) )/2.0D0
    YPU(JM)=Y(JM)+( Y(JM)-Y(JM-1) )/2.0D0

    RETURN
    END SUBROUTINE