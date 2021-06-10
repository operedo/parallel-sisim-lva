! C functions declaration
interface
    subroutine dijkstra_c(n,a,b,c,d,e,f,g) bind(C, name="dijkstra")
        use iso_c_binding
        implicit none
        integer(C_INT), intent(in) :: n
        integer(C_INT), intent(in) :: a
        integer(C_INT), intent(in) :: b(*)
        integer(C_INT), intent(in) :: c(*)
        real(C_DOUBLE), intent(in) :: d(*) 
        integer(C_INT), intent(in) :: e
        integer(C_INT), intent(in) :: f(*)
        real(C_DOUBLE), intent(in) :: g(*) 
    end subroutine


end interface
