    !######################################################################!
    !#                                                                    #!
    !#                              功能子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !******************************************计算网格规模*****************************************************!
    SUBROUTINE CAL_GRID_SIZE
    USE DECLARATION
    IMPLICIT NONE
    REAL(KIND=8)::X1,Y1    

    IF(CONTINUOUS_MESH==0)THEN
        !非连续网格
        ILM1=1   +IDNINT((LEM1-LEFT)/DX1 )
        ILM2=ILM1+IDNINT((LEM2-LEM1)/DX21)
        ILM3=ILM2+IDNINT((LEM3-LEM2)/DX22)
        IL  =ILM3+IDNINT((LEIN-LEM3)/DX23)
        IR  =IL  +IDNINT((RIIN-LEIN)/DX3 )
        IRM3=IR  +IDNINT((RIM3-RIIN)/DX23)
        IRM2=IRM3+IDNINT((RIM2-RIM3)/DX22)
        IRM1=IRM2+IDNINT((RIM1-RIM2)/DX21)
        IM  =IRM1+IDNINT((RIGH-RIM1)/DX1 )

        JBM1=1   +IDNINT((BOM1-BOTT)/DX1 )
        JBM2=JBM1+IDNINT((BOM2-BOM1)/DX21)
        JBM3=JBM2+IDNINT((BOM3-BOM2)/DX22)
        JB  =JBM3+IDNINT((BOIN-BOM3)/DX23)
        JT  =JB  +IDNINT((TOIN-BOIN)/DX3 )
        JTM3=JT  +IDNINT((TOM3-TOIN)/DX23)
        JTM2=JTM3+IDNINT((TOM2-TOM3)/DX22)
        JTM1=JTM2+IDNINT((TOM1-TOM2)/DX21)
        JM  =JTM1+IDNINT((TOPP-TOM1)/DX1 )
    ELSE IF(CONTINUOUS_MESH==1)THEN
        !连续网格
        HL=LEFT-LEIN
        HR=RIGH-RIIN
        HB=BOTT-BOIN
        HT=TOPP-TOIN

        !NODEL=63
        !NODER=88
        !NODEB=55
        !NODET=55

        NODEL=1
        DO WHILE( DABS(HL)>CRITERIA )
            NODEL=NODEL+10
            X1=LEIN+HL*( BL+1.0D0-(BL-1.0D0)*( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEL-1) ) )/( ( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEL-1) )+1.0D0 )
            IF( LEIN-X1 <=DX3 )EXIT
        END DO
        DO WHILE( DABS(HL)>CRITERIA )
            NODEL=NODEL-1
            X1=LEIN+HL*( BL+1.0D0-(BL-1.0D0)*( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEL-1) ) )/( ( (BL+1.0D0)/(BL-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEL-1) )+1.0D0 )
            IF( LEIN-X1 >=DX3 )EXIT
        END DO

        NODER=1
        DO WHILE( DABS(HR)>CRITERIA )
            NODER=NODER+10
            X1=RIIN+HR*( BR+1.0D0-(BR-1.0D0)*( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODER-1) ) )/( ( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODER-1) )+1.0D0 )
            IF( X1-RIIN <=DX3 )EXIT
        END DO
        DO WHILE( DABS(HR)>CRITERIA )
            NODER=NODER-1
            X1=RIIN+HR*( BR+1.0D0-(BR-1.0D0)*( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODER-1) ) )/( ( (BR+1.0D0)/(BR-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODER-1) )+1.0D0 )
            IF( X1-RIIN >=DX3 )EXIT
        END DO

        NODEB=1
        DO WHILE( DABS(HB)>CRITERIA )
            NODEB=NODEB+10
            Y1=BOIN+HB*( BB+1.0D0-(BB-1.0D0)*( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEB-1) ) )/( ( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEB-1) )+1.0D0 )
            IF( BOIN-Y1 <=DX3 )EXIT
        END DO
        DO WHILE( DABS(HB)>CRITERIA )
            NODEB=NODEB-1
            Y1=BOIN+HB*( BB+1.0D0-(BB-1.0D0)*( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEB-1) ) )/( ( (BB+1.0D0)/(BB-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODEB-1) )+1.0D0 )
            IF( BOIN-Y1 >=DX3 )EXIT
        END DO

        NODET=1
        DO WHILE( DABS(HT)>CRITERIA )
            NODET=NODET+10
            Y1=TOIN+HT*( BT+1.0D0-(BT-1.0D0)*( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODET-1) ) )/( ( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODET-1) )+1.0D0 )
            IF( Y1-TOIN <=DX3 )EXIT
        END DO
        DO WHILE( DABS(HT)>CRITERIA )
            NODET=NODET-1
            Y1=TOIN+HT*( BT+1.0D0-(BT-1.0D0)*( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODET-1) ) )/( ( (BT+1.0D0)/(BT-1.0D0) )**( 1.0D0-DBLE(1)/DBLE(NODET-1) )+1.0D0 )
            IF( Y1-TOIN >=DX3 )EXIT
        END DO
        WRITE(*,*)NODEL,NODER,NODEB,NODET

        IL=NODEL
        IR=IL+IDNINT((RIIN-LEIN)/DX3 )
        IM=IR+NODER-1

        JB=NODEB
        JT=JB+IDNINT((TOIN-BOIN)/DX3 )
        JM=JT+NODET-1
    END IF

    IIM=IR-IL+1
    JIM=JT-JB+1
    WRITE(*,*)IM,JM

    RETURN
    END