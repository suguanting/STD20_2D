    !######################################################################!
    !#                                                                    #!
    !#                              定义模块                              #!
    !#                                                                    #!
    !######################################################################!
    MODULE DECLARATION

    REAL(KIND=8),PARAMETER::PI=3.1415926535897932384626433832795D0
    !------算例种类------!
    INTEGER::CASE_TYPE
    !------任务种类------!
    INTEGER::TASK_TYPE
    !------内边界运动情况------!
    INTEGER::IB_LOCOMOTION
    !------粘性项计算方法------!
    INTEGER::VISCOUS_TERM_METHOD
    !------网格种类------!
    INTEGER::CONTINUOUS_MESH
    !------动边界形状------!
    INTEGER::IB_SHAPE
    !------连续网格控制量------!
    REAL(KIND=8)::HL,HR,HB,HT
    REAL(KIND=8)::BL,BR,BB,BT
    INTEGER::NODEL,NODER,NODEB,NODET
    !------计算域尺度量（边界只出现在网格最密处）------!
    REAL(KIND=8)::LEFT=-7.5D0,LEM1=-4.0D0,LEM2=-2.0D0,LEM3=-1.6D0,LEIN=-1.7D0!-0.7D0!-1.2D0
    REAL(KIND=8)::RIGH=15.0D0,RIM1= 4.0D0,RIM2= 2.0D0,RIM3= 1.6D0,RIIN= 2.3D0! 0.7D0! 1.2D0
    REAL(KIND=8)::BOTT=-7.5D0,BOM1=-5.0D0,BOM2=-4.0D0,BOM3=-3.0D0,BOIN=-3.0D0!-1.0D0!-2.5D0
    REAL(KIND=8)::TOPP= 7.5D0,TOM1= 5.0D0,TOM2= 4.0D0,TOM3= 3.0D0,TOIN= 3.0D0! 1.0D0! 2.5D0
    !------网格密度量------!
    REAL(KIND=8)::DX1 !外层
    REAL(KIND=8)::DX21!中外层
    REAL(KIND=8)::DX22!中中层
    REAL(KIND=8)::DX23!中内层
    REAL(KIND=8)::DX3 !内层
    !------各种无量纲参数和准则数------!
    REAL(KIND=8)::Re
    !------时间信息------!
    REAL(KIND=8)::T,DT
    !------计时器------!
    REAL(KIND=8)::WALLCLOCK1,WALLCLOCK2,WALLCLOCK
    REAL(KIND=8)::TARRAY(2)
    !------迭代控制------!(权宜之计，在判断交点时可能需要准则放大)(这个备注什么鬼？？)
    INTEGER::NSTEP,NSTART,NMAX,NCYCLE,NDURATION,NSUBSTEP
    REAL(KIND=8)::ERRORVELOMAX,EUMAX,EVMAX,VELOMAX,UMAX,VMAX,PMAX
    INTEGER::PLOC(2),ULOC(2),EULOC(2),VLOC(2),EVLOC(2),WLOC(2),EWLOC(2)
    REAL(KIND=8)::CRITERIA=1.0D-6
    INTEGER::CONVERGENCE=0
    !------续算文件名------!
    CHARACTER(LEN=80)::FILENAME_RESTART
    !------调用输出次数控制------!
    INTEGER::NPROBE,NCLCT,NIB
    REAL(KIND=8)::NPLT
    !------边界条件系数------!BC_A*U+BC_B
    REAL(KIND=8)::BCU_AL,BCU_AR,BCU_AB,BCU_AT
    REAL(KIND=8)::BCU_BL,BCU_BR,BCU_BB,BCU_BT
    REAL(KIND=8)::BCU_CL,BCU_CR,BCU_CB,BCU_CT
    REAL(KIND=8)::BCV_AL,BCV_AR,BCV_AB,BCV_AT
    REAL(KIND=8)::BCV_BL,BCV_BR,BCV_BB,BCV_BT
    REAL(KIND=8)::BCV_CL,BCV_CR,BCV_CB,BCV_CT
    REAL(KIND=8)::BCPHI_AL,BCPHI_AR,BCPHI_AB,BCPHI_AT
    REAL(KIND=8)::BCPHI_BL,BCPHI_BR,BCPHI_BB,BCPHI_BT
    INTEGER::BCTYPE_L,BCTYPE_R,BCTYPE_B,BCTYPE_T
    !------格点数------!
    INTEGER::IM,JM
    INTEGER::ILM1,ILM2,ILM3,IL,IR,IRM3,IRM2,IRM1
    INTEGER::JBM1,JBM2,JBM3,JB,JT,JTM3,JTM2,JTM1
    INTEGER::IIM,JIM
    INTEGER::I,J
    !------固定格点信息------!
    REAL(KIND=8),ALLOCATABLE::X(:),Y(:)
    REAL(KIND=8),ALLOCATABLE::XPV(:),YPU(:)
    REAL(KIND=8),ALLOCATABLE::U(:,:),V(:,:),P(:,:),PHI(:,:),UHAT(:,:),VHAT(:,:)
    !------该时间层(N)和上一时间层(N-1)的速度，A-B格式用，计算残差用------!
    REAL(KIND=8),ALLOCATABLE::UN(:,:),VN(:,:)
    REAL(KIND=8),ALLOCATABLE::UN1(:,:),VN1(:,:)
    !------该子时间层(K)和上一子时间层(K-1)的速度，上二子时间层(K-2)的速度，R-K格式用------!
    REAL(KIND=8),ALLOCATABLE::UK(:,:),VK(:,:)
    REAL(KIND=8),ALLOCATABLE::UK1(:,:),VK1(:,:)
    REAL(KIND=8),ALLOCATABLE::UK2(:,:),VK2(:,:)
    !------临时计数器------!
    INTEGER::N,M,MT,MB!N为循环计数器
    !INTEGER::MM,MM1!对应组按限制条件筛选时的临时计数器
    !------升阻力系数（其他算例时需更改）------!
    REAL(KIND=8)::CX1,CY1,CX2,CY2
    !REAL(KIND=8),ALLOCATABLE::ZREL_STAGE(:)
    !------储存变量名的字符变量------!
    !CHARACTER(LEN=10)::VARNAME
    !------探针坐标值------!
    REAL(KIND=8)::PROBE_X1,PROBE_X2,PROBE_X3,PROBE_X4
    REAL(KIND=8)::PROBE_Y1,PROBE_Y2,PROBE_Y3,PROBE_Y4
    INTEGER::PROBE_IPV1,PROBE_IPV2,PROBE_IPV3,PROBE_IPV4
    INTEGER::PROBE_JPU1,PROBE_JPU2,PROBE_JPU3,PROBE_JPU4
    INTEGER::PROBE_IU1,PROBE_IU2,PROBE_IU3,PROBE_IU4
    INTEGER::PROBE_JV1,PROBE_JV2,PROBE_JV3,PROBE_JV4
    !------RK系数------!GAMMA写为GAMA
    REAL(KIND=8)::GAMA(3),RHO(3),ALPHA(3)

    !$!
    !------旋转的特殊参数------!
    !背景速度场的旋转角，影响：速度初场和速度边界条件
    REAL(KIND=8)::FREESTREAM_TILT!为零代表沿计算域X+，逆时针旋转为正，不参与旋转矩阵的计算
    !由计算域绝对坐标系X+指向上拍方向，影响：PSI-确定拍动平面角；MAT_ABS2FLP-由计算域绝对坐标系到扑翼坐标系，用于扑翼气动阻力和升力计算
    !出现在CAL_CLCT_ELLIPTIC POSE_VELO_QUADRIC_2D_PERIODIC/STATIC/INTERMITTENT_FORE/INTERMITTENT_HIND
    REAL(KIND=8)::ABSX_UPSTROKE_ANGLE!逆时针方向为正
    !飞行角，由真实世界X+指向飞行方向，影响：计算速度在X、Y方向分量，用于有效功率计算
    !出现在CAL_CLCT_ELLIPTIC（CAL_CLCT_CIRCULAR）
    REAL(KIND=8)::TRUX_FLIGHT_ANGLE!逆时针方向为正
    !由计算域绝对坐标系X+指向真实世界X+，影响：MAT_ABS2TRU-由计算域绝对坐标系到真实世界坐标系，用于真实气动力计算
    !出现在CAL_CLCT_ELLIPTIC
    REAL(KIND=8)::ABSX_TRUX_ANGLE!逆时针方向为正
    !前进比
    REAL(KIND=8)::J_FORW
    !飞行速度与特征速度比值
    REAL(KIND=8)::VELO_RATIO!不一定等于前进比
    REAL(KIND=8)::U_FREESTREAM,V_FREESTREAM

    END MODULE
    !**********************************************浸入式边界模块所需变量*************************************************!
    MODULE IMMERSED_BOUNDARY

    !------DERIVATIVE VARIABLES------!
    !------XSWEEP------!
    REAL(KIND=8),ALLOCATABLE::IB_RUX(:,:),IB_RVX(:,:)!aiVi
    REAL(KIND=8),ALLOCATABLE::IB_AUX(:,:),IB_AVX(:,:)!bi
    REAL(KIND=8),ALLOCATABLE::IB_CUX(:,:),IB_CVX(:,:)!
    !------YSWEEP------!
    REAL(KIND=8),ALLOCATABLE::IB_RUY(:,:),IB_RVY(:,:)!ajVj
    REAL(KIND=8),ALLOCATABLE::IB_BUY(:,:),IB_BVY(:,:)!bj
    REAL(KIND=8),ALLOCATABLE::IB_CUY(:,:),IB_CVY(:,:)!

    !------PRIMITIVE VARIABLES------!
    !------XSWEEP------!
    INTEGER     ,ALLOCATABLE::TYPEUX(:,:),TYPEVX(:,:)!10-FLOW
    REAL(KIND=8),ALLOCATABLE::IB_ITSCT_UX(:,:),IB_ITSCT_VX(:,:)!INTERSECTION LOCATION
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UX(:,:),IB_IPSVL_VX(:,:)!IMPOSED VELOCITY
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UXV(:,:),IB_IPSVL_VXU(:,:)!IMPOSED VELOCITY
    !------YSWEEP------!
    INTEGER     ,ALLOCATABLE::TYPEUY(:,:),TYPEVY(:,:)!10-FLOW
    REAL(KIND=8),ALLOCATABLE::IB_ITSCT_UY(:,:),IB_ITSCT_VY(:,:)!INTERSECTION LOCATION
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UY(:,:),IB_IPSVL_VY(:,:)!IMPOSED VELOCITY
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UYV(:,:),IB_IPSVL_VYU(:,:)!IMPOSED VELOCITY

    !------上一时间层------!
    REAL(KIND=8),ALLOCATABLE::LAPLACE_VXK1(:,:),LAPLACE_VYK1(:,:)!
    REAL(KIND=8),ALLOCATABLE::LAPLACE_UXK1(:,:),LAPLACE_UYK1(:,:)!
    REAL(KIND=8),ALLOCATABLE::CONVECT_VXK1(:,:),CONVECT_VYK1(:,:)!
    REAL(KIND=8),ALLOCATABLE::CONVECT_UXK1(:,:),CONVECT_UYK1(:,:)!
    REAL(KIND=8),ALLOCATABLE::CONVECT_VXK2(:,:),CONVECT_VYK2(:,:)!
    REAL(KIND=8),ALLOCATABLE::CONVECT_UXK2(:,:),CONVECT_UYK2(:,:)!

    INTEGER     ,ALLOCATABLE::TYPEUXN(:,:),TYPEVXN(:,:)!10:FLOW,-10:SOLID
    INTEGER     ,ALLOCATABLE::TYPEUYN(:,:),TYPEVYN(:,:)!10:FLOW,-10:SOLID
    REAL(KIND=8),ALLOCATABLE::IB_ITSCT_UXN(:,:),IB_ITSCT_VXN(:,:)!INTERSECTION LOCATION
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UXN(:,:),IB_IPSVL_VXN(:,:)!IMPOSED VELOCITY
    REAL(KIND=8),ALLOCATABLE::IB_ITSCT_UYN(:,:),IB_ITSCT_VYN(:,:)!INTERSECTION LOCATION
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UYN(:,:),IB_IPSVL_VYN(:,:)!IMPOSED VELOCITY

    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UXVN(:,:),IB_IPSVL_VXUN(:,:)!IMPOSED VELOCITY 只起到储存数据的作用，不实际用到
    REAL(KIND=8),ALLOCATABLE::IB_IPSVL_UYVN(:,:),IB_IPSVL_VYUN(:,:)!IMPOSED VELOCITY 只起到储存数据的作用，不实际用到
    !------本时间层在上一时间层的速度（由插值等方法得到）------!
    REAL(KIND=8),ALLOCATABLE::IBN1_IPSVL_UXN(:,:),IBN1_IPSVL_VXN(:,:)
    REAL(KIND=8),ALLOCATABLE::IBN1_IPSVL_UYN(:,:),IBN1_IPSVL_VYN(:,:)

    END MODULE
    !**********************************************边界运动信息*************************************************!
    !增加边界时，需在此处增加几何和运动参数定义，然后CAL_QUADRIC、TYPE_IB_REINITIATE、CAL_CLCT_CIRCULAR中复制代码，然后IBM_QUADRIC中多调用一次INTERSECTION函数,IBM_INITIATION也需要
    MODULE QUADRIC_PARAMETER

    REAL(KIND=8)::XREL,YREL
    !展向位置
    REAL(KIND=8)::SPAN

    INTEGER::BOUNDARY_EXISTENCE_1,BOUNDARY_EXISTENCE_2
    !------边界1------!
    REAL(KIND=8)::QUADRIC_GEOMETRICAL_1(36),QUADRIC_KINETIC_1(9)
    !------边界2------!
    REAL(KIND=8)::QUADRIC_GEOMETRICAL_2(36),QUADRIC_KINETIC_2(9)

    !------上一时间层的边界信息------!
    !------边界1------!
    REAL(KIND=8)::QUADRIC_GEOMETRICAL_N1(36),QUADRIC_KINETIC_N1(9)
    !------边界2------!
    REAL(KIND=8)::QUADRIC_GEOMETRICAL_N2(36),QUADRIC_KINETIC_N2(9)

    END MODULE

    !******************************************计算QUADIRC边界运动情况时的临时数据*********************************************!
    !非主程序模组，只应用于CAL_QUADIRC相关子函数中
    MODULE CAL_QUADRIC_DECLARATION

    !三个欧拉角和一个拍动平面角
    REAL(KIND=8)::PSI!ψ拍动平面夹角，此类算例中应为0°
    REAL(KIND=8)::PHIW!ϕw拍动角
    REAL(KIND=8)::PSIW!ψw翻转角
    REAL(KIND=8)::PSI0!ψ0初始翻转角
    REAL(KIND=8)::THETAW!θw偏离角/偏移角
    REAL(KIND=8)::PHASE_DIFFERENCE!前后翅相位差,谁大谁提前
    REAL(KIND=8)::PHASE_INITIATION!前后翅相位差,谁大谁提前
    REAL(KIND=8)::TAU_INIT!由起始相位角推得起始时间

    !时间和角度定量
    REAL(KIND=8)::TAUC!周期时长
    REAL(KIND=8)::TAU!周期内时刻
    REAL(KIND=8)::TAU_R1,TAU_R2!两个翻转开始的时刻
    REAL(KIND=8)::DTAUR!单次翻转时长
    REAL(KIND=8)::PHIM!最大拍动角
    REAL(KIND=8)::PHIA!拍动角平均值
    REAL(KIND=8)::PSIM!最大攻角

    !坐标转换矩阵
    REAL(KIND=8)::TRANMAT(2,2)
    REAL(KIND=8)::T11,T12,T21,T22
    !坐标转换逆矩阵
    REAL(KIND=8)::TRANMAT_INVERSE(2,2)

    !平动速度（相对+绝对）
    REAL(KIND=8)::VELO_TRAN_R(2),VELO_TRAN_A(2)
    !转动中心(相对于旋转平面的坐标+绝对坐标）
    REAL(KIND=8)::CEN_P(2),CEN(2)
    REAL(KIND=8)::CEN_TRANSLATION(2)!各边界转动中心平移量
    REAL(KIND=8)::CEN_DEVIATION(2)!各边界转动中心在弦向和拍动向相对二次图形几何图形中心偏移量
    !转动速度
    REAL(KIND=8)::VELO_ANGL

    !平移后绝对坐标系下平面的数学表达式系数
    REAL(KIND=8)::COX2,COY2,COXY,COX,COY,COM!,COZ2,COXZ,COYZ,COZ
    REAL(KIND=8)::LAXIS,SAXIS

    !------是否闭合，1闭合，0非闭合------!
    INTEGER::CLOSED

    !------非无限平面四个角绝对坐标系限定条件------!
    !上表面
    REAL(KIND=8)::XMAXT,XMINT
    REAL(KIND=8)::YMAXT,YMINT
    !下表面
    REAL(KIND=8)::XMAXB,XMINB
    REAL(KIND=8)::YMAXB,YMINB
    !------中间数------!
    !单位坐标转换矩阵（拍动平面角和翻转角）
    REAL(KIND=8)::MATP(2,2),MATW(2,2)
    !相对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::RCOX2,RCOY2,RCOXY,RCOX,RCOY,RCOM!,RCOZ2,RCOXZ,RCOYZ,RCOZ
    !未平移绝对坐标系下二次曲面的数学表达式系数
    REAL(KIND=8)::PCOX2,PCOY2,PCOXY,PCOX,PCOY,PCOM!,PCOZ2,PCOXZ,PCOYZ,PCOZ
    !------非无限平面四个角相对坐标系限定条件------!
    !上下表面X,Z共同使用，Y分别用一个，为COM/C
    REAL(KIND=8)::RXTRL,RXLED
    REAL(KIND=8)::RYTOP,RYBOT
    !------非无限平面四个角绝对坐标值------!
    REAL(KIND=8)::X1,X2
    REAL(KIND=8)::Y1,Y2

    !------间歇性飞行相关------!
    REAL(KIND=8)::DTAUH,DTAUI,DTAUL
    REAL(KIND=8)::TSTART,TSTOP,TRESTART,TEND
    REAL(KIND=8)::PHII,PSII,TAUI
    !间歇飞行起始所处周期数
    REAL(KIND=8)::INTERMITTENT_INITIATION

    !------HEAVING&PLUNGING相关------!
    REAL(KIND=8)::K,H,KH

    END MODULE