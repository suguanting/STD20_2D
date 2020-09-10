    !2020/09/10
    !!$$$!
    !REAL(KIND=8),ALLOCATABLE::TRIAD(:),TRIBD(:),TRICD(:),TRIDD(:)
    !!$$$!
    !!$$$!
    !IF( ALLOCATED(TRIAD) )THEN
    !    DEALLOCATE(TRIAD,TRIBD,TRICD,TRIDD)
    !END IF
    !ALLOCATE( TRIAD(IM-1),TRIBD(IM-1),TRICD(IM-1),TRIDD(IM-1) )
    !TRIAD=0.0D0
    !TRIBD=0.0D0
    !TRICD=0.0D0
    !TRIDD=0.0D0
    !!$$$!
    !!$$$!
                !IF(I==174 .AND. J==158)THEN
                !    WRITE(*,*)I
                !END IF
                !!$$$!
    !!$$$!
            !DO I=2,IM-2,1
            !    IF(TYPEVX(I,J)==10)THEN
            !        TRIAD(I)=A1(I,J)
            !        TRIBD(I)=1.0D0-A2(I,J)
            !        TRICD(I)=A3(I,J)
            !        TRIDD(I)=RV(I,J)
            !    ELSE IF(TYPEVX(I,J)==0)THEN
            !        TRIAD(I)=0.0D0
            !        TRIBD(I)=1.0D0
            !        TRICD(I)=0.0D0
            !        TRIDD(I)=IB_RVX(I,J)
            !    ELSE IF(TYPEVX(I,J)==1)THEN
            !        TRIAD(I)=0.0D0
            !        TRIBD(I)=1.0D0+IB_CVX(I,J)
            !        TRICD(I)=IB_AVX(I,J)
            !        TRIDD(I)=RV(I,J)+IB_RVX(I,J)
            !    ELSE IF(TYPEVX(I,J)==-1)THEN
            !        TRIAD(I)=IB_AVX(I,J)
            !        TRIBD(I)=1.0D0+IB_CVX(I,J)
            !        TRICD(I)=0.0D0
            !        TRIDD(I)=RV(I,J)+IB_RVX(I,J)
            !    ELSE IF(TYPEVX(I,J)==-10)THEN
            !        TRIAD(I)=0.0D0
            !        TRIBD(I)=1.0D0
            !        TRICD(I)=0.0D0
            !        TRIDD(I)=IB_RVX(I,J)
            !    END IF
            !    !$$$!
            !    IF(I==174 .AND. J==158)THEN
            !        WRITE(*,*)I
            !    END IF
            !    !$$$!
            !END DO
            !
            !WRITE(34567,*)'J',J,':',MAXVAL(DABS(TRIAD(2:IM-2)-TRIA(2:IM-2)))
            !WRITE(34567,*)'J',J,':',MAXVAL(DABS(TRIBD(2:IM-2)-TRIB(2:IM-2)))
            !WRITE(34567,*)'J',J,':',MAXVAL(DABS(TRICD(2:IM-2)-TRIC(2:IM-2)))
            !WRITE(34567,*)'J',J,':',MAXVAL(DABS(TRIDD(2:IM-2)-TRID(2:IM-2)))
            !IF(DABS(Y(J))<=0.5)THEN
            !    WRITE(34567,*)'有固体点'
            !END IF
            !
            !!$$$!
    !!$$$!
    !IF( ALLOCATED(TRIAD) )THEN
    !    DEALLOCATE(TRIAD,TRIBD,TRICD,TRIDD)
    !END IF
    !ALLOCATE( TRIAD(2:JM-1),TRIBD(2:JM-1),TRICD(2:JM-1),TRIDD(2:JM-1) )
    !TRIAD=0.0D0
    !TRIBD=0.0D0
    !TRICD=0.0D0
    !TRIDD=0.0D0
    !!$$$!
    !!$$$!
            !DO J=3,JM-2,1
            !    IF(TYPEVY(I,J)==10)THEN
            !        TRIA(J)=B1(I,J)
            !        TRIB(J)=1.0D0-B2(I,J)
            !        TRIC(J)=B3(I,J)
            !        TRID(J)=DVP(I,J)
            !    ELSE IF(TYPEVY(I,J)==0)THEN
            !        TRIA(J)=0.0D0
            !        TRIB(J)=1.0D0
            !        TRIC(J)=0.0D0
            !        TRID(J)=IB_RVY(I,J)
            !    ELSE IF(TYPEVY(I,J)==1)THEN
            !        TRIA(J)=0.0D0
            !        TRIB(J)=1.0D0+IB_CVY(I,J)
            !        TRIC(J)=IB_BVY(I,J)
            !        TRID(J)=DVP(I,J)+IB_RVY(I,J)
            !    ELSE IF(TYPEVY(I,J)==-1)THEN
            !        TRIA(J)=IB_BVY(I,J)
            !        TRIB(J)=1.0D0+IB_CVY(I,J)
            !        TRIC(J)=0.0D0
            !        TRID(J)=DVP(I,J)+IB_RVY(I,J)
            !    ELSE IF(TYPEVY(I,J)==-10)THEN
            !        TRIA(J)=0.0D0
            !        TRIB(J)=1.0D0
            !        TRIC(J)=0.0D0
            !        TRID(J)=IB_RVY(I,J)
            !    END IF
            !END DO
            !
            !WRITE(34561,*)'I',I,':',MAXVAL(DABS(TRIAD(3:JM-2)-TRIA(3:JM-2)))
            !WRITE(34561,*)'I',I,':',MAXVAL(DABS(TRIBD(3:JM-2)-TRIB(3:JM-2)))
            !WRITE(34561,*)'I',I,':',MAXVAL(DABS(TRICD(3:JM-2)-TRIC(3:JM-2)))
            !WRITE(34561,*)'I',I,':',MAXVAL(DABS(TRIDD(3:JM-2)-TRID(3:JM-2)))
            !
            !!$$$!
    !!$$$!
    !IF( ALLOCATED(TRIAD) )THEN
    !    DEALLOCATE(TRIAD,TRIBD,TRICD,TRIDD)
    !END IF
    !ALLOCATE( TRIAD(2:IM-1),TRIBD(2:IM-1),TRICD(2:IM-1),TRIDD(2:IM-1) )
    !TRIAD=0.0D0
    !TRIBD=0.0D0
    !TRICD=0.0D0
    !TRIDD=0.0D0
    !!$$$!
    !!$$$!
    !IF( ALLOCATED(TRIAD) )THEN
    !    DEALLOCATE(TRIAD,TRIBD,TRICD,TRIDD)
    !END IF
    !ALLOCATE( TRIAD(JM-1),TRIBD(JM-1),TRICD(JM-1),TRIDD(JM-1) )
    !TRIAD=0.0D0
    !TRIBD=0.0D0
    !TRICD=0.0D0
    !TRIDD=0.0D0
    !!$$$!
    !$$$!
    !CALL OUTPUT_IB_STAGGERED
    !$$$!
    
        !2020/09/10正规化气动力求解
    !!2
    !CV_LEFT=-0.7D0
    !CV_RIGH= 0.7D0
    !CV_BOTT=-0.7D0
    !CV_TOPP= 0.7D0
    !
    !!求解最近坐标
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!力初始化
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!求解各项面力
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!求解非定常
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!合力
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!计算域绝对坐标系下的力
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30002,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL

    !!3
    !CV_LEFT=-0.6D0
    !CV_RIGH= 0.9D0
    !CV_BOTT=-0.6D0
    !CV_TOPP= 0.9D0
    !
    !!求解最近坐标
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!力初始化
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!求解各项面力
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!求解非定常
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!合力
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!计算域绝对坐标系下的力
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30003,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !
    !!4
    !CV_LEFT=-0.8D0
    !CV_RIGH= 0.8D0
    !CV_BOTT=-0.8D0
    !CV_TOPP= 0.8D0
    !
    !!求解最近坐标
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!力初始化
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!求解各项面力
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!求解非定常
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!合力
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!计算域绝对坐标系下的力
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30004,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !
    !!5
    !CV_LEFT=-0.9D0
    !CV_RIGH= 0.9D0
    !CV_BOTT=-0.9D0
    !CV_TOPP= 0.9D0
    !
    !!求解最近坐标
    !DISTANCE_X=X(:)-CV_LEFT
    !INDEX_XL=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XL)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVL SEARCH ERROR"
    !DISTANCE_X=X(:)-CV_RIGH
    !INDEX_XR=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    !IF( DABS(X(INDEX_XR)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVR SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_BOTT
    !INDEX_YB=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YB)-CV_LEFT)>CRITERIA ) WRITE(*,*)"CVB SEARCH ERROR"
    !DISTANCE_Y=Y(:)-CV_TOPP
    !INDEX_YT=MINLOC( DABS(DISTANCE_Y),1 )+LBOUND(DISTANCE_Y,1)-1
    !IF( DABS(Y(INDEX_YT)-CV_RIGH)>CRITERIA ) WRITE(*,*)"CVT SEARCH ERROR"
    !
    !!力初始化
    !CXC_TOTAL=0.0D0
    !CYC_TOTAL=0.0D0
    !
    !CFTX=0.0D0
    !CFTY=0.0D0
    !
    !CFCLX=0.0D0
    !CFCLY=0.0D0
    !CFCRX=0.0D0
    !CFCRY=0.0D0
    !CFCBX=0.0D0
    !CFCBY=0.0D0
    !CFCTX=0.0D0
    !CFCTY=0.0D0
    !
    !CFPL=0.0D0
    !CFPR=0.0D0
    !CFPB=0.0D0
    !CFPT=0.0D0
    !
    !CFVLX=0.0D0
    !CFVLY=0.0D0
    !CFVRX=0.0D0
    !CFVRY=0.0D0
    !CFVBX=0.0D0
    !CFVBY=0.0D0
    !CFVTX=0.0D0
    !CFVTY=0.0D0
    !
    !CFCX=0.0D0
    !CFCY=0.0D0
    !CFPX=0.0D0
    !CFPY=0.0D0
    !CFVX=0.0D0
    !CFVY=0.0D0
    !
    !!求解各项面力
    !I=INDEX_XL
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCLX=CFCLX+DX3*U(I,J)*U(I,J)
    !    CFCLY=CFCLY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPL =CFPL +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVLX=CFVLX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVLY=CFVLY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !I=INDEX_XR
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    CFCRX=CFCRX+DX3*U(I,J)*U(I,J)
    !    CFCRY=CFCRY+DX3*U(I,J)*0.25D0*(V(I,J)+V(I,J+1)+V(I-1,J)+V(I-1,J+1))
    !    CFPR =CFPR +DX3*0.5D0*(P(I,J)+P(I-1,J))
    !    CFVRX=CFVRX+DX3/Re*2.0D0*(U(I+1,J)-U(I-1,J))/(X(I+1)-X(I-1))
    !    CFVRY=CFVRY+DX3/Re*( (U(I,J+1)-U(I,J-1))/(YPU(J+1)-YPU(J-1)) &
    !        & + (V(I,J)+V(I,J+1)-V(I-1,J)-V(I-1,J+1))/(XPV(I)-XPV(I-1))/2.0D0 )
    !END DO
    !J=INDEX_YB
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCBX=CFCBX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCBY=CFCBY+DX3*V(I,J)*V(I,J)
    !    CFPB =CFPB +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVBX=CFVBX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVBY=CFVBY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !J=INDEX_YT
    !DO I=INDEX_XL,INDEX_XR-1,1
    !    CFCTX=CFCTX+DX3*V(I,J)*0.25D0*(U(I,J)+U(I,J-1)+U(I+1,J)+U(I+1,J-1))
    !    CFCTY=CFCTY+DX3*V(I,J)*V(I,J)
    !    CFPT =CFPT +DX3*0.5D0*(P(I,J)+P(I,J-1))
    !    CFVTX=CFVTX+DX3/Re*( (V(I+1,J)-V(I-1,J))/(XPV(I+1)-XPV(I-1)) &
    !        & + (U(I,J)+U(I+1,J)-U(I,J-1)-U(I+1,J-1))/(YPU(J)-YPU(J-1))/2.0D0 )
    !    CFVTY=CFVTY+DX3/Re*2.0D0*(V(I,J+1)-V(I,J-1))/(Y(J+1)-Y(J-1))
    !END DO
    !!求解非定常
    !DO J=INDEX_YB,INDEX_YT-1,1
    !    DO I=INDEX_XL,INDEX_XR-1,1
    !        CFTX=CFTX+DX3*DX3*(U(I+1,J)+U(I,J)-UN(I+1,J)-UN(I,J))/2.0D0/DT
    !        CFTY=CFTY+DX3*DX3*(V(I,J+1)+V(I,J)-VN(I,J+1)-VN(I,J))/2.0D0/DT
    !    END DO
    !END DO
    !
    !!合力
    !CFCX=-CFCLX+CFCRX-CFCBX+CFCTX
    !CFCY=-CFCLY+CFCRY-CFCBY+CFCTY
    !CFPX=-CFPL+CFPR
    !CFPY=-CFPB+CFPT
    !CFVX=-CFVLX+CFVRX-CFVBX+CFVTX
    !CFVY=-CFVLY+CFVRY-CFVBY+CFVTY
    !
    !!计算域绝对坐标系下的力
    !CXC_TOTAL=-2.0D0*(CFTX+CFCX+CFPX-CFVX)
    !CYC_TOTAL=-2.0D0*(CFTY+CFCY+CFPY-CFVY)
    !
    !WRITE(30005,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL

    !6
    
!***************************************************求解升力推力（控制体方法,简化）******************************************************!
    SUBROUTINE CAL_CLCT_CONTROL_VOLUME_SIMPLIFIED
    USE DECLARATION
    IMPLICIT NONE

    !---------输出相关---------!
    CHARACTER(LEN=6)CHAR_STEP,CHAR_REYNOLDS
    INTEGER::REYNOLDS

    !控制体区域
    REAL(KIND=8)::CV_LEFT,CV_RIGH,CV_BOTT,CV_TOPP

    !寻址
    REAL(KIND=8),ALLOCATABLE::DISTANCE_X(:)
    INTEGER::INDEX_XL_U,INDEX_XR_U
    INTEGER::INDEX_XL_P,INDEX_XR_P

    !界面信息
    REAL(KIND=8),ALLOCATABLE::U_LEFT_BOUNDARY(:),U_RIGH_BOUNDARY(:)
    REAL(KIND=8),ALLOCATABLE::P_LEFT_BOUNDARY(:),P_RIGH_BOUNDARY(:)
    REAL(KIND=8),ALLOCATABLE::BOUNDARY_AREA(:)

    !------力系数------!
    REAL(KIND=8)::CXC_TOTAL,CYC_TOTAL!计算域绝对坐标系下的力


    ALLOCATE( BOUNDARY_AREA(JM-1) )

    ALLOCATE( DISTANCE_X(IM) )
    ALLOCATE( U_LEFT_BOUNDARY(JM-1),U_RIGH_BOUNDARY(JM-1) )
    ALLOCATE( P_LEFT_BOUNDARY(JM-1),P_RIGH_BOUNDARY(JM-1) )


    CV_LEFT=-1.0D0
    CV_RIGH= 1.0D0
    CV_BOTT=BOTT
    CV_TOPP=TOPP


    !求解最近坐标
    DISTANCE_X=X(:)-CV_LEFT
    INDEX_XL_U=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( X(INDEX_XL_U)>CV_LEFT ) INDEX_XL_U=INDEX_XL_U-1
    DISTANCE_X=X(:)-CV_RIGH
    INDEX_XR_U=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( X(INDEX_XR_U)>CV_RIGH ) INDEX_XR_U=INDEX_XR_U-1
    !求解最近坐标
    DISTANCE_X=XPV(1:IM)-CV_LEFT
    INDEX_XL_P=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( XPV(INDEX_XL_P)>CV_LEFT ) INDEX_XL_P=INDEX_XL_P-1
    DISTANCE_X=XPV(1:IM)-CV_RIGH
    INDEX_XR_P=MINLOC( DABS(DISTANCE_X),1 )+LBOUND(DISTANCE_X,1)-1
    IF( XPV(INDEX_XR_P)>CV_RIGH ) INDEX_XR_P=INDEX_XR_P-1

    !插值速度和压力
    DO J=1,JM-1,1
        CALL LINEAR_INTERPOLATION(U(INDEX_XL_U,J),U_LEFT_BOUNDARY(J),U(INDEX_XL_U+1,J),X  (INDEX_XL_U),CV_LEFT,X  (INDEX_XL_U+1))
        CALL LINEAR_INTERPOLATION(U(INDEX_XR_U,J),U_RIGH_BOUNDARY(J),U(INDEX_XR_U+1,J),X  (INDEX_XR_U),CV_RIGH,X  (INDEX_XR_U+1))
        CALL LINEAR_INTERPOLATION(P(INDEX_XL_P,J),P_LEFT_BOUNDARY(J),P(INDEX_XL_P+1,J),XPV(INDEX_XL_P),CV_LEFT,XPV(INDEX_XL_P+1))
        CALL LINEAR_INTERPOLATION(P(INDEX_XR_P,J),P_RIGH_BOUNDARY(J),P(INDEX_XR_P+1,J),XPV(INDEX_XR_P),CV_RIGH,XPV(INDEX_XR_P+1))
    END DO

    !------求边界大小------!
    DO J=1,JM-1,1
        BOUNDARY_AREA(J)=Y(J+1)-Y(J)
    END DO

    CXC_TOTAL=-2.0D0*SUM((U_RIGH_BOUNDARY*U_RIGH_BOUNDARY+P_RIGH_BOUNDARY)*BOUNDARY_AREA-(U_LEFT_BOUNDARY*U_LEFT_BOUNDARY+P_LEFT_BOUNDARY)*BOUNDARY_AREA)


    WRITE(30001,"( I6,2(1X,F9.5)                       )")NSTEP,CXC_TOTAL,CYC_TOTAL
    !WRITE(30001,*)NSTEP,CXC_TOTAL,CYC_TOTAL

    RETURN
    END SUBROUTINE