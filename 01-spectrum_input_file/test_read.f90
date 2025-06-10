program test_read

   implicit none
   integer, parameter :: VLI_K = selected_int_kind (18)
   integer, parameter :: DR_K = selected_real_kind (14)

   integer (VLI_K) :: i,iw
   real (DR_K) :: a, b, c, d
   real wf_bry(1550)
   real wa_bry(1550)

   open (unit=15, file="input_spectrum.txt", status='old',    &
             access='sequential', form='formatted', action='read' )

   do iw=1,1550
     read(15,*) wf_bry(iw), wa_bry(iw)
   enddo  
   
   do iw=1,1550
     write(*,*) wf_bry(iw), wa_bry(iw)
   enddo  

end program test_read