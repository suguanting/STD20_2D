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
            !    WRITE(34567,*)'ÓÐ¹ÌÌåµã'
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