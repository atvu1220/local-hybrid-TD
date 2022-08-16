program testmpi
      use mpi
      implicit none
      integer:: i,my_rank,procnum,ierr, num_loops, flag,status, remain
      integer, parameter:: nz=98
      integer:: sent(nz), transform(nz) 
      integer, allocatable:: procarray(:), recvbuf(:), tail(:) 
      
      do i=1,nz
            sent(i) = i
      enddo
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      num_loops = nz/procnum
      remain=mod(nz,procnum)
      
      allocate(procarray(num_loops))
      flag = 1
      do i = my_rank*num_loops + 1, num_loops*(my_rank + 1)
            procarray(flag) = sent(i)
            flag = flag + 1
      enddo
   
      if ((remain /= 0) .and. (my_rank .eq. procnum-1)) then
            allocate(tail(remain))
            flag = 1
            do i = num_loops*procnum + 1, nz
    
                  tail(flag) = sent(i)
                  flag = flag +1
            enddo
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (my_rank .eq. 0) allocate(recvbuf(num_loops*procnum))
      call MPI_GATHER(procarray,num_loops,MPI_INT,recvbuf,num_loops,MPI_INT,0,MPI_COMM_WORLD,ierr)
      deallocate(procarray)

      if (remain /= 0) then
            allocate(procarray(remain))

            if (my_rank .eq. procnum-1) then
                  call MPI_SEND(tail,remain,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
                  deallocate(tail)
            endif
            if (my_rank .eq. 0) then
                  call MPI_RECV(procarray,remain,MPI_INTEGER,procnum-1,1,MPI_COMM_WORLD,status,ierr)
            endif
!            call MPI_SENDRECV(tail,remain,MPI_INTEGER,0,1,procarray,remain,MPI_INTEGER,procnum-1,1,MPI_COMM_WORLD,status,ierr)
      endif

      if (my_rank .eq. 0) then
            do i=1,procnum*num_loops
                  transform(i) = recvbuf(i)
            enddo
            if (remain /= 0) then
                  flag =1
                  do i=procnum*num_loops+1,nz
                        transform(i) = procarray(flag)
                        flag = flag +1
                  enddo
                  deallocate(procarray)
            endif
            deallocate(recvbuf)
            write(*,*) transform
      endif
!      deallocate(procarray,recvbuf,tail)
      
    
      call MPI_FINALIZE(ierr)
      write(*,*) 'End program'
            
      
     

end program testmpi

      