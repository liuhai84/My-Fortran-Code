module TimeClass
!
!            时间类
!
  implicit none
  !type definition
  type , public :: Time
    character(len = 8)  :: data , time , zone
    integer(kind = 8)   :: values(8)
    real(kind = 8)      :: Start , Finish
   contains
   !Bound procedures
    procedure , private :: GetDataTime   => GetDataTimeSub    !获取日期时间
    procedure , public :: GetStartTime  => GetStartTimeSub   !程序开始时间
    procedure , public :: GetFinishTime => GetFinishTimeSub  !程序结束时间
    procedure , private :: GetUseTime    => GetUseTimeSub
    procedure , public :: to_string_fn
  end type Time
  contains
  !Now add methods
   character(len = 10) function to_string_fn(this)
    !
    ! Represent the data as a string : MM/DD/YYYY
    !
    implicit none
    !Declare calling arguments
    class(Time) :: this     !date object
      !Declare local variables
    character(len = 2)  :: dd       !Day
    character(len = 2)  :: mm       !month
    character(len = 4)  :: yy       !yy
    !get data
    call this%GetDataTime()
    !Get components
    write(dd , '(I2.2)') this%values(3)
    write(mm , '(I2.2)') this%values(2)
    write(yy , '(I4)')   this%values(1)
    !return string
    to_string_fn = yy // '/' // mm // '/' // dd

   end function to_string_fn

   subroutine GetDataTimeSub(this)
    !----------------------------------------------------------------------
    !--------                  获取日期时间                    ------------
    !----------------------------------------------------------------------
     implicit none
     class(Time) , intent(inout) :: this
     !**********************************************************************
     call date_and_time(this%data , this%time , this%zone , this%values)

     !**********************************************************************
   end subroutine GetDataTimeSub
   !*************************************************************************
   subroutine GetStartTimeSub(this)
    !----------------------------------------------------------------------
    !--------                  程序开始时间                    ------------
    !----------------------------------------------------------------------
    implicit none
    class(Time) , intent(inout) :: this
    !**********************************************************************
    call cpu_time(this%Start)
    !***********************************************************************
   end subroutine GetStartTimeSub
    !*************************************************************************
   subroutine GetFinishTimeSub(this)
    !----------------------------------------------------------------------
    !--------                  程序结束时间                    ------------
    !----------------------------------------------------------------------
    implicit none
    class(Time) , intent(inout) :: this
    !**********************************************************************
    call cpu_time(this%Finish)
    call this%GetUseTime()
    !***********************************************************************
   end subroutine GetFinishTimeSub

   subroutine GetUseTimeSub(this)
    !----------------------------------------------------------------------
    !--------                  程序结束时间                    ------------
    !----------------------------------------------------------------------
    implicit none
    class(Time) , intent(inout) :: this
    write(*,'(1x,A20 , F5.3 , A1)') 'Execution of time = ' , this%Finish - this%Start , 's'
   end subroutine
end module TimeClass

module ConfigClass
  !
  ! 参数类，读取参数
  !
  implicit none
  !type definition
  type , public :: Config
   ! private
   real(kind = 8)               :: Mu , K                      !水文地质参数，Mu for 弹性给水度，K for 导水系数
   real(kind = 8)               :: deltaX , deltaZ , deltaT    !deltaX for △x , deltaZ for △Z , delta for △t
   integer(kind = 8)            :: Nx , Nz                     !单元格划分，Nx for x_direction , Nz for z_direction
  end type Config
  !****************************************************************************************************************
  interface Config
     module procedure CreateConfig
  end interface Config
  !****************************************************************************************************************
  contains
  type(Config) function CreateConfig(cfgFilePath)
   implicit none
   character(len = *) , intent(in) :: cfgFilePath
   integer(kind = 8)               :: cfgFileID
   logical                         :: aLive
   !**************************************************************************************************************
   inquire(file = cfgFilePath , exist = aLive)
   if(aLive) then
    open( newunit=cfgFileID, file=trim(cfgFilePath), status='old', action='read' )
    read(cfgFileID , *) CreateConfig%Mu
    read(cfgFileID , *) CreateConfig%K
    read(cfgFileID , *) CreateConfig%deltaX
    read(cfgFileID , *) CreateConfig%deltaZ
    read(cfgFileID , *) CreateConfig%deltaT
    read(cfgFileID , *) CreateConfig%Nx
    read(cfgFileID , *) CreateConfig%Nz
    close(cfgFileID)
   else
    write(*,*) trim(cfgFilePath) , " doesn't exit."
    stop
   end if
  end function CreateConfig
  !***************************************************************************************************************
end module ConfigClass
!****************************************************************
module Steady_flow_two_dimensional_implicit_difference_method_Class
  !
  !   二维隐式差分法
  !
  use ConfigClass
  use TimeClass
  implicit none
  !Type definition
  type , public :: mSolver
  !Instance variables
   type(Config)                  :: mConfig
   type(Time)                    :: mTime
   real(kind = 8) , allocatable  :: Head(:,:) , Head0(:,:)          !水头值，Head for Final , head0 for Initial
   real(kind = 8)                :: DfMax , Time = 0.0 , Err = 0.00001
  contains
   !!Bound procedures
   procedure , public :: GetData           => GetDataSub           !获得参数
   procedure , public :: GetAllocatedArray => GetAllocatedArraySub !分配数组空间
   procedure , public :: Lambda            => GetLambdaFun         !计算Lambda值
   procedure , public :: DeleteArray       => DeleteArraySub       !释放数组内存
   procedure , public :: GetInitialValue   => GetInitialValueSub   !赋初值
   procedure , public :: GetBoundaryValue  => GetBoundaryValueSub  !获得边界值
   procedure , public :: GetVersion        => GetVersionSub        !程序信息
   procedure , public :: WriteData         => WriteDataSub         !写入计算结果
   procedure , public :: GetHeadToHead0    => GetHeadToHead0Sub    !本次计算结果为下一次迭代的初值
   procedure , public :: GetHeadValue      => GetHeadValueSub      !计算水头值
   procedure , public :: CalculationDfMax  => CalculationDfMaxSub  !计算迭代误差
   procedure , public :: IsStop            => IsStopSub            !误差达到，程序结束
   procedure , public :: Run               => RunSub               !运行程序
  end type mSolver
  !Now add methods
  contains
  !*************************************************************************
   subroutine GetDataSub(this , cfgFilePath)
    !----------------------------------------------------------------------
    !--------                  获取参数                        ------------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(out)        :: this
    character(len = *) , intent(inout) :: cfgFilePath
    !**********************************************************************
     this%mConfig = config(cfgFilePath)  !获取参数
    !**********************************************************************
    call this%GetAllocatedArray()
   end subroutine GetDataSub
 !*************************************************************************
   subroutine GetAllocatedArraySub(this)
    !----------------------------------------------------------------------
    !--------------           分配数组空间                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout) :: this
    !**********************************************************************
    allocate(this%Head(this%mConfig%Nx , this%mConfig%Nz) , this%Head0(this%mConfig%Nx , this%mConfig%Nz))
    !**********************************************************************
   end subroutine GetAllocatedArraySub

 !*************************************************************************
   real(kind = 8) function GetLambdaFun(this)
     implicit none
     class(mSolver) , intent(inout) :: this
     GetLambdaFun = this%mConfig%K * this%mConfig%deltaT / (this%mConfig%Mu * this%mConfig%deltaX * this%mConfig%deltaX)
   end function GetLambdaFun

   subroutine DeleteArraySub(this)
    !----------------------------------------------------------------------
    !--------------           释放数组内存                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout) :: this
    integer(kind = 8)             :: istat
    !**********************************************************************
    write(*,*) 'In Finalizer...'

    deallocate(this%Head0 , stat = istat)
    deallocate(this%Head , stat = istat)

    !**********************************************************************
   end subroutine DeleteArraySub
 !*************************************************************************
   subroutine GetInitialValueSub(this)
    !----------------------------------------------------------------------
    !--------------           获取水头初值                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    integer(kind = 8)                :: i , j
    !**********************************************************************
    !every where
    do i = 2 , this%mConfig%Nx - 1
      do j = 2 , this%mConfig%Nz - 1
        this%Head0(i , j) = 0.0
      end do
    end do
    !**********************************************************************
    ! up point 2 to 7
    do i = 2 , 7
      this%Head(i , 1)  = 5.0
      this%Head0(i , 1) = 5.0
    end do
    !**********************************************************************
    !up point 9 to end
    do i = 9 , this%mConfig%Nx - 1
      this%Head(i , 1) = 0.0
    end do
    !**********************************************************************
   end subroutine GetInitialValueSub
 !*************************************************************************
   subroutine GetBoundaryValueSub(this)
    !----------------------------------------------------------------------
    !--------------           获取边界初值                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    integer(kind = 8)                :: i , j
    !**********************************************************************
    !left and right
    do j = 2 , this%mConfig%Nz - 1
      this%Head0(1 , j) = this%Head0(3 , j)
      this%Head0(this%mConfig%Nx , j) = this%Head0(this%mConfig%Nx - 2 , j)
    end do
    !**********************************************************************
    ! bottom
    do i = 2 , this%mConfig%Nx - 1
      this%Head0(i , this%mConfig%Nz) = this%Head0(i , this%mConfig%Nz - 2)
    end do
    !**********************************************************************
   end subroutine GetBoundaryValueSub
 !*************************************************************************
  subroutine GetVersionSub(this)
    !----------------------------------------------------------------------
    !--------------           程序信息输入                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout) :: this
    !******************************************************************************************
    write(*,*)'*******************************************************************************'
    write(*,*)'------------        求解稳定流问题的二维隐式差分法计算机程序       ------------'
    write(*,*)'------------       Code By:LH  Data:',this%mTime%to_string_fn(),'  Version:1.0.0','      ------------'
    write(*,*)'*******************************************************************************'
 !*************************************************************************
   end subroutine GetVersionSub
!*************************************************************************
   subroutine WriteDataSub(this)
    !----------------------------------------------------------------------
    !--------------           计算结果写入                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    integer(kind = 8)                :: i , j
    !**********************************************************************
    call this%GetVersion()
    write(* , 36) this%Time
    36 format(1x,'Time = ' , F6.2)
    write(* , 37) ((this%Head(i , j),i = 2 , this%mConfig%Nx - 1) , j = 1 , this%mConfig%Nz - 1)
    37 format(1x , 'Head(Ni , Nj)=',/,(13F6.2, 1x))
    write(*,*) '*******************************************************************************'
    !**********************************************************************
   end subroutine WriteDataSub
 !*************************************************************************
   subroutine GetHeadToHead0Sub(this)
    !----------------------------------------------------------------------
    !--------------           终值变为初值                      -----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    integer(kind = 8)                :: i , j
    !**********************************************************************
    do i = 2 , this%mConfig%Nx - 1
      do j = 2 , this%mConfig%Nz - 1
        this%Head0(i , j) = this%Head(i , j)
      end do
    end do
    this%Head0(8 , 1) = this%Head(8 , 1)
    !**********************************************************************
   end subroutine GetHeadToHead0Sub
 !*************************************************************************
   subroutine GetHeadValueSub(this)
    !----------------------------------------------------------------------
    !--------------           计算区域水头值                     ----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    integer(kind = 8)                :: i , j
    !**********************************************************************
    this%Time = this%Time + this%mConfig%deltaT
    do i = 2 , this%mConfig%Nx - 1
      do j = 2 , this%mConfig%Nz - 1
        this%Head(i , j) = (this%Head0(i - 1 , j) + this%Head0(i + 1 , j) + &
                           this%Head0(i , j - 1) + this%Head0(i , j + 1)) / 4.0
      end do
    end do
    this%Head(8 , 1) = (this%Head0(7 , 1) + this%Head0(9 , 1) +             &
                        2 * this%Head0(8 , 2)) / 4.0
    !**********************************************************************
   end subroutine GetHeadValueSub
 !*************************************************************************
   subroutine CalculationDfMaxSub(this)
    !----------------------------------------------------------------------
    !--------------          计算结束误差                       ----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout)    :: this
    real(kind = 8)                   :: AbsErr
    integer(kind = 8)                :: i , j
    !**********************************************************************
    this%DfMax = 0.0
    do i = 2 , this%mConfig%Nx - 1
      do j = 1 , this%mConfig%Nz - 1
        AbsErr = abs(this%Head(i , j) - this%Head0(i , j))
        if(AbsErr > this%DfMax) this%DfMax = AbsErr
      end do
    end do
  !  if(this%DfMax < this%Err) exit
    !**********************************************************************
   end subroutine CalculationDfMaxSub
    !**********************************************************************`
   subroutine IsStopSub(this)
    !----------------------------------------------------------------------
    !--------------           水头计算结束                       ----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout) :: this
    !**********************************************************************
    if(this%DfMax < this%Err) then
      call this%WriteData()
      call this%DeleteArray()
      call this%mTime%GetFinishTime()
     ! call this%mTime%GetUseTime()
      stop
    end if
   end subroutine IsStopSub
!**************************************************************************
   subroutine RunSub(this , cfgFilePath)
    !----------------------------------------------------------------------
    !--------------           水头计算结束                       ----------
    !----------------------------------------------------------------------
    implicit none
    class(mSolver) , intent(inout) :: this
    character(len = 200) :: cfgFilePath
    call this%GetData(cfgFilePath)
    call this%GetInitialValue()
  !**************************************************************
    do
      call this%GetBoundaryValue()
      call this%GetHeadValue()
      call this%CalculationDfMax()
       if(this%DfMax < this%Err) exit
      call this%GetHeadToHead0()
     end do
     call this%WriteData()
     call this%DeleteArray()
   end subroutine RunSub

end module Steady_flow_two_dimensional_implicit_difference_method_Class

program main
  use Steady_flow_two_dimensional_implicit_difference_method_Class
  implicit none
  !***************************************************************
  type(mSolver)       :: SteadyFlowForTwoDim
  character(len=200) :: configFile = "data.txt"
  !***************************************************************
  !begin of program
  call SteadyFlowForTwoDim%mTime%GetStartTime()   !获取初始运行时间
  !***************************************************************
  call SteadyFlowForTwoDim%Run(configFile)
  !***************************************************************
  call SteadyFlowForTwoDim%mTime%GetFinishTime()  !获取程序运行时间
  ! end of program
end program main
