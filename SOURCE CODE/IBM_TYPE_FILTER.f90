    !######################################################################!
    !#                                                                    #!
    !#                              步骤子函数                            #!
    !#                                                                    #!
    !######################################################################!

    !*************************************权宜之计：去除相邻两向外点******************************************!
    SUBROUTINE IBM_TYPE_FILTER
    USE IMMERSED_BOUNDARY
    USE DECLARATION

    !------U------!
    DO J=0,JM-1,1
        DO I=1,IM-1,1
            !------X------!
            IF( TYPEUX(I,J)==-1 .AND. TYPEUX(I,J)*TYPEUX(I+1,J)==-1 )THEN
                TYPEUX(I  ,J)=10
                TYPEUX(I+1,J)=10
            END IF
            !------Y------!
            IF( TYPEUY(I,J)==-1 .AND. TYPEUY(I,J)*TYPEUY(I,J+1)==-1 )THEN
                TYPEUY(I,J  )=10
                TYPEUY(I,J+1)=10
            END IF

        END DO
    END DO

    !------V------!
    DO J=1,JM-1,1
        DO I=0,IM-1,1
            !------X------!
            IF( TYPEVX(I,J)==-1 .AND. TYPEVX(I,J)*TYPEVX(I+1,J)==-1 )THEN
                TYPEVX(I  ,J)=10
                TYPEVX(I+1,J)=10
            END IF
            !------Y------!
            IF( TYPEVY(I,J)==-1 .AND. TYPEVY(I,J)*TYPEVY(I,J+1)==-1 )THEN
                TYPEVY(I,J  )=10
                TYPEVY(I,J+1)=10
            END IF

        END DO
    END DO

    RETURN
    END SUBROUTINE

