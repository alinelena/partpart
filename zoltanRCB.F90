module zoltanRCB
  use mpi
  use zoltan
  implicit none

  private

  type,public ::  particlesType 
    integer                          :: nGParticles, nLParticles
    integer(ZOLTAN_INT), allocatable :: GIDs(:)
    real(ZOLTAN_DOUBLE), allocatable :: x(:), y(:), z(:)
    character(len=6),    allocatable :: labels(:)  
   end type particlesType
   type(particlesType),public :: particles

   type,public :: zoltType
     LOGICAL                      :: changes 
     INTEGER(Zoltan_INT)          :: numGidEntries, numLidEntries
     INTEGER(Zoltan_INT)          :: numImport, numExport
     INTEGER(Zoltan_INT), POINTER :: importGlobalGids(:), exportGlobalGids(:)
     INTEGER(Zoltan_INT), POINTER :: importLocalGids(:), exportLocalGids(:)
     INTEGER(Zoltan_INT), POINTER :: importProcs(:), exportProcs(:)
     INTEGER(Zoltan_INT), POINTER :: importToPart(:), exportToPart(:)
   end type zoltType
   type(zoltType), public :: zolt

   public :: partitionParticlesWithRCB
   public :: zoltanCleanup

  contains

  subroutine partitionParticlesWithRCB()

    type(Zoltan_Struct), pointer :: zz
    integer(ZOLTAN_INT)          :: ierr

    nullify(zz)

    zz=>Zoltan_Create(MPI_COMM_WORLD)

    ierr = Zoltan_Set_Param(zz, "LB_METHOD", "RCB")
    ierr = Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE")

    ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, zoltNumObjs)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, zoltGetObjs)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, zoltNumGeom)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_GEOM_FN_TYPE, zoltGeom)
!    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE, zoltObjSize)
!    ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE, zoltObjPack)
!    ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE, zoltObjUnPack)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Use Zoltan to partition the vertices in the simple mesh.
    !!
    !! Params:
    !!     zz               -- input (all remaining fields are output)
    !!     changes          -- 1 if partition was changed, 0 otherwise 
    !!     numGidEntries    -- Number of integers used for a global ID 
    !!     numLidEntries    -- Number of integers used for a local ID 
    !!     numImport        -- Number of vertices to be sent to me 
    !!     importGlobalGids -- Global IDs of vertices to be sent to me 
    !!     importLocalGids  -- Local IDs of vertices to be sent to me 
    !!     importProcs      -- Process rank for source of each incoming vertex 
    !!     importToPart     -- New part for each incoming vertex 
    !!     numExport        -- Number of vertices I must send to other processes
    !!     exportGlobalGids -- Global IDs of the vertices I must send 
    !!     exportLocalGids  -- Local IDs of the vertices I must send 
    !!     exportProcs      -- Process to which I send each of the vertices 
    !!     exportToPart     -- Part to which each vertex will belong 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ierr = Zoltan_LB_Partition(zz, zolt%changes, zolt%numGidEntries, zolt%numLidEntries, &
                               zolt%numImport, zolt%importGlobalGids, zolt%importLocalGids, zolt%importProcs, zolt%importToPart, &
                               zolt%numExport, zolt%exportGlobalGids, zolt%exportLocalGids, zolt%exportProcs, zolt%exportToPart)


    call Zoltan_Destroy(zz)

  end subroutine partitionParticlesWithRCB

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Frees arrays allocated by Zoltan_LB_Partition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltanCleanup()

    integer :: error

    error = Zoltan_LB_Free_Part(zolt%importGlobalGids, zolt%importLocalGids, zolt%importProcs, zolt%importToPart)
    error = Zoltan_LB_Free_Part(zolt%exportGlobalGids, zolt%exportLocalGids, zolt%exportProcs, zolt%exportToPart)

  end subroutine zoltanCleanup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumObjs(data, ierr)
    ! Local declarations
    integer(Zoltan_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    zoltNumObjs = particles%nLParticles
    ierr = ZOLTAN_OK

  end function zoltNumObjs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGetObjs (data, nGIDs, nLIDs, GIDs, LIDs, wgt_dim, obj_wgts, ierr)
    integer(ZOLTAN_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(in)  :: nGIDs, nLIDs
    integer(ZOLTAN_INT), intent(out) :: GIDs(*), LIDs(*)
    integer(ZOLTAN_INT), intent(in)  :: wgt_dim 
    real(ZOLTAN_FLOAT), intent(out)  :: obj_wgts(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    integer :: i
    do i= 1, particles%nLParticles
      GIDs(i) = particles%GIDs(i)
      LIDs(i) = i
    end do
    
    ierr = ZOLTAN_OK

  end subroutine zoltGetObjs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumGeom(data, ierr)
    integer(Zoltan_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT)              :: ierr

    zoltNumGeom = 3
    ierr = ZOLTAN_OK
  end function zoltNumGeom

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGeom(data, nGIDs, nLIDs, GID, LID, geom_vec, ierr)
    integer(Zoltan_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(in)  :: nGIDs, nLIDs, GID, LID
    real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    geom_vec(1) =  particles%x(LID)
    geom_vec(2) =  particles%y(LID)
    geom_vec(3) =  particles%z(LID)

    ierr = ZOLTAN_OK

  end subroutine zoltGeom
  function zoltObjSize(data, nGIDS,  nLIDs, GIDs, LIDs, ierr)
    integer(Zoltan_INT)              :: zoltObjSize
    integer(Zoltan_INT), intent(in)  :: data(*)
    INTEGER(Zoltan_INT), INTENT(IN)  :: nGIDS,  nLIDs
    INTEGER(Zoltan_INT), INTENT(IN)  :: GIDs(*), LIDs(*)
    INTEGER(Zoltan_INT), INTENT(OUT) :: ierr

    zoltObjSize=3*nLIDs*storage_size(particles%x(1))/8
    ierr=ZOLTAN_OK
  end function zoltObjSize

  subroutine zoltObjPack(data, num_gid_entries, num_lid_entries, global_id, local_id, dest, size, buf, ierr)
    INTEGER(Zoltan_INT), INTENT(IN)  :: data  
    INTEGER(Zoltan_INT), INTENT(IN)  :: num_gid_entries, num_lid_entries  
    INTEGER(Zoltan_INT), INTENT(IN)  :: global_id(*)  
    INTEGER(Zoltan_INT), INTENT(IN)  :: local_id(*)  
    INTEGER(Zoltan_INT), INTENT(IN)  :: dest, size  
    INTEGER(Zoltan_INT), INTENT(OUT) :: buf(*)
    INTEGER(Zoltan_INT), INTENT(OUT) :: ierr 
    integer :: c=0
    real(Zoltan_DOUBLE) :: r(3)
    c=c+1
    print *,c,dest,size,local_id(1)
    
    ierr=ZOLTAN_OK
  end subroutine zoltObjPack
  subroutine zoltObjUnPack(data, num_gid_entries, global_id, size, buf, ierr)
    INTEGER(Zoltan_INT),INTENT(IN)   :: data(*)  
    INTEGER(Zoltan_INT), INTENT(IN)  :: num_gid_entries  
    INTEGER(Zoltan_INT), INTENT(IN)  :: global_id(*)  
    INTEGER(Zoltan_INT), INTENT(IN)  :: size  
    INTEGER(Zoltan_INT), INTENT(IN)  :: buf(*)
    INTEGER(Zoltan_INT), INTENT(OUT) :: ierr 
    integer :: c=0
    c=c+1
    !print *,c,buf(1:size),size
    ierr=ZOLTAN_OK
  end subroutine zoltObjUnPack
end module zoltanRCB
