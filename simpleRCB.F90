program simpleRCB
   use mpi
   use zoltan
   use zoltanRCB

   implicit none

   integer(Zoltan_INT) :: error
   real(Zoltan_FLOAT)  :: version

   interface

     subroutine readInputObjects(fname,ps)
        use zoltanRCB 
        implicit none
        character (len=*), intent(in)     :: fname
        type(particlesType),intent(inout) :: ps
     end subroutine readInputObjects

   end interface

   call MPI_Init(error)
   error = Zoltan_Initialize(version)

   call readInputObjects("particles.xyz",particles)
   print *, particles%gids
   call partitionParticlesWithRCB()

   call visualizePartition(particles)

   deallocate(particles%labels)
   deallocate(particles%z,particles%y,particles%x)
   deallocate(particles%GIDs)
  
   call zoltanCleanUp()
   call MPI_Finalize(error)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program simpleRCB

subroutine readInputObjects(fname, ps)
  use mpi
  use zoltanRCB 
  implicit none

  character (len=*),intent(in) :: fname
  type(particlesType),intent(inout) :: ps

  integer :: fnum = 101, i, currIndx
  integer :: myRank, numProcs, mpi_ierr
  character(len=6) :: label
  real :: tmpX, tmpY,tmpZ

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpi_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, mpi_ierr)

  open(unit = fnum, file = fname)
  read(fnum,*) ps%nGParticles

  ps%nLParticles=0

  do i = 1, ps%nGParticles, 1
    !! assumes i start at 1, gives round robin initial distribution
    if ( MOD(i-1,numProcs) == myRank) then
      ps%nLParticles = ps%nLParticles + 1
    end if
  end do

  allocate(ps%GIDs(ps%nLParticles))
  allocate(ps%x(ps%nLParticles))
  allocate(ps%y(ps%nLParticles))
  allocate(ps%z(ps%nLParticles))
  allocate(ps%labels(ps%nLParticles))

  read(fnum,*) label

  currIndx = 1
  do i = 1, ps%nGParticles, 1
    read(fnum,*) label, tmpX, tmpY,tmpZ

    if ( MOD(i-1,numProcs) == myRank) then
      ps%GIDs(currIndx) = i
      ps%labels(currIndx) = trim(label)
      ps%x(currIndx) = tmpX
      ps%y(currIndx) = tmpY
      ps%z(currIndx) = tmpZ
      currIndx = currIndx + 1
    end if
  end do

  close(fnum)
end subroutine readInputObjects

subroutine showSimpleParticlePartitions(myProc,parts,ps,label)
  use mpi
  use zoltanRCB
  implicit none
  integer,intent(in)              :: myProc
  integer,intent(in)              :: parts(*)
  type(particlesType), intent(in) :: ps
  character(len=*),intent(in)     :: label

  !! Local variables
  integer :: i, j, part, ierr
  do i=1,ps%nLParticles
     write(*,*)label,ps%GIDs(i),parts(i)
  enddo 
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
end subroutine showSimpleParticlePartitions

subroutine visualizePartition(ps)
  use mpi
  use zoltanRCB
  implicit none
  type(particlesType),intent(inout) :: ps

  integer :: parts(ps%nLParticles)
  integer :: myrank, i, error

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, error)


  do i=1, ps%nLParticles, 1
    parts(i) = myRank;
  end do
  if (myRank== 0) then
    write (*,*) 'Particles assignments before calling Zoltan'
  end if
  call showSimpleParticlePartitions(myRank, parts,ps,"B");
  do i=1, zolt%numExport, 1
    parts(zolt%exportLocalGids(i)) = zolt%exportToPart(i)
  end do


  if (myRank == 0) then
    write (*,*) 'Particles assignments after calling Zoltan'
  end if

  call showSimpleParticlePartitions(myRank,parts,ps,"A")

end subroutine visualizePartition
