program test_read

   implicit none
   integer, parameter :: VLI_K = selected_int_kind (18)
   integer, parameter :: DR_K = selected_real_kind (14)

   integer (VLI_K) :: i,iw
   real (DR_K) :: a, b, c, d
   real wf_bry(1550)
   real wa_bry(1550)
   real df

   !open (unit=15, file="FLASH_RIP/input_spectrum.txt", status='old',    &
   !          access='sequential', form='formatted', action='read' )
   open(15,file='./FLASH_RIP/input_spectrum.txt',form='formatted',status='old')

   do iw=1,1550
     read(15,*) wf_bry(iw), wa_bry(iw)
     if (iw.gt.1) then
       df=wf_bry(iw)-wf_bry(iw-1)
     endif
   enddo
   write(*,*) df
   
   !do iw=1,1550
   !  write(*,*) wf_bry(iw), wa_bry(iw)
   !enddo

end program test_read
