!--------------------------------------------------------------------------------------------------------------------------------------

subroutine remove_integer(list, element)
    implicit none

    !Declaramos los argumentos de entrada y salida.
    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    if (allocated(list)) then
            if (size(list) > 1) then
                list = pack(list, mask=list /= element)
            else
                deallocate(list)
            end if
    end if
    end subroutine remove_integer
    
    subroutine add_integer(list, element)
    implicit none

    !Declaramos los argumentos de entrada y de salida.
    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    !Declaramos la copia de la lista.
    integer, dimension(:), allocatable :: copy_list

    !Declaramos las variables que vamos a usar en la subrutina.
    integer :: isize

    if (allocated(list)) then
        isize = size(list)
        allocate(copy_list(isize+1))

        copy_list(1:isize) = list
        copy_list(1+isize) = element

        deallocate(list)
        call move_alloc(copy_list, list)
    else
        allocate(list(1))
        list(1) = element
    end if
    end subroutine add_integer

!--------------------------------------------------------------------------------------------------------------------------------------

INTEGER FUNCTION findloc(list, element)
    implicit none
    integer, dimension(:), allocatable, intent(in) :: list
    integer, intent(in) :: element
    integer i

    findloc = 0

    if (allocated(list)) then
        do i = 1, size(list)
            if (list(i) .eq. element) then
            findloc = i
            exit
            endif
        enddo
    endif
END FUNCTION findloc


!--------------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE ORDENA(list)
    implicit none
    interface
        subroutine remove_integer(list, element)
            integer, dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: element
        end subroutine remove_integer
        subroutine add_integer(list, element)
            integer, dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: element
        end subroutine add_integer
    end interface
    integer, dimension(:), allocatable, intent(inout) :: list
    integer, dimension(:), allocatable :: list1, list2
    integer i, j

    allocate(list1(size(list)))
    list1 = list

    do i = 1, size(list)
        j = minval(list1)
        call add_integer(list2, j)
        call remove_integer(list1, j)
    enddo

    list = list2
END SUBROUTINE ORDENA