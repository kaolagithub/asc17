subroutine gather_out(ncid, filename, action, args)
    implicit none
    character(len=*), intent(in) :: filename, action, args
    integer, intent(in) :: ncid

    integer :: ierr, stat(mpi_status_size)

    real, allocatable :: t_v3(:), tt_v3(:), tt_loc3(:, :, :)
    integer :: t_loc2(2), t_loc3(3)
    integer :: t_ixl, t_iyl, t_kb, t_is, t_ie, t_js, t_je
    integer :: t_im, t_ixoverlay, pex
    integer :: max_ixl, max_iyl, max_kb, max_loc3(3)
    integer, dimension(npe) :: ixoverlayss, imss, ixlxx, iylxx
    integer, dimension(npe) :: kbss, isss, jsss, iess, jess
    integer, dimension(npe) :: t_loc2ss(2), t_loc3ss(3)

    integer :: ncount, index, tt_is, tt_js, tt_kb
    if (myid /= 0) then
        ncount = ixl * iyl * kb
        allocate(t_v3(ncount), tt_v3(ncount))
        call mpi_send(ixoverlay, 1, mpi_integer, 0, 298, mpi_comm_ympi, ierr)
        call mpi_send(im, 1, mpi_integer, 0, 299, mpi_comm_ympi, ierr)
        call mpi_send(ixl, 1, mpi_integer, 0, 300, mpi_comm_ympi, ierr)
        call mpi_send(iyl, 1, mpi_integer, 0, 301, mpi_comm_ympi, ierr)
        call mpi_send(kb , 1, mpi_integer, 0, 302, mpi_comm_ympi, ierr)
        call mpi_send(is , 1, mpi_integer, 0, 303, mpi_comm_ympi, ierr)
        call mpi_send(js , 1, mpi_integer, 0, 304, mpi_comm_ympi, ierr)
        call mpi_send(ie , 1, mpi_integer, 0, 305, mpi_comm_ympi, ierr)
        call mpi_send(je , 1, mpi_integer, 0, 306, mpi_comm_ympi, ierr)
        call mpi_send(loc2(1), 2, mpi_integer, 0, 307, mpi_comm_ympi, ierr)
        call mpi_send(loc3(1), 3, mpi_integer, 0, 308, mpi_comm_ympi, ierr)
    else
        max_ixl = ixl
        max_iyl = iyl
        max_kb  = kb
        max_loc3(1) = ie - is + loc3(1)
        max_loc3(2) = je - js + loc3(2)
        max_loc3(3) = kb -  1 + loc3(3)
        if(ixoverlay /= 0 .and. loc2(1) == 1)then
            max_loc3(1) = max(max_loc3(1), im)
        endif
        do pex = 1, npe - 1
            call mpi_recv(t_ixoverlay, 1, mpi_integer, pex, 298, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_im, 1, mpi_integer, pex, 299, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_ixl, 1, mpi_integer, pex, 300, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_iyl, 1, mpi_integer, pex, 301, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_kb , 1, mpi_integer, pex, 302, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_is , 1, mpi_integer, pex, 303, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_js , 1, mpi_integer, pex, 304, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_ie , 1, mpi_integer, pex, 305, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_je , 1, mpi_integer, pex, 306, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_loc2(1), 2, mpi_integer, pex, 307, mpi_comm_ympi, stat, ierr)
            call mpi_recv(t_loc3(1), 3, mpi_integer, pex, 308, mpi_comm_ympi, stat, ierr)
            ixoverlayss(pex) = t_ixoverlay
            imss(pex)        = t_im
            ixlxx(pex)       = t_ixl
            iylxx(pex)       = t_iyl
            kbss(pex)        = t_kb
            isss(pex)        = t_is
            jsss(pex)        = t_js
            iess(pex)        = t_ie
            jess(pex)        = t_je
            t_loc2ss(pex, :) = t_loc2(:)
            t_loc3ss(pex, :) = t_loc3(:)
            max_ixl = max(max_ixl, t_ixl)
            max_iyl = max(max_iyl, t_iyl)
            max_kb  = max(max_kb , t_kb )
            max_loc3(1) = max(max_loc3(1), t_ie - t_is + t_loc3(1))
            max_loc3(2) = max(max_loc3(2), t_je - t_je + t_loc3(2))
            max_loc3(3) = max(max_loc3(3), t_kb - 1 + t_loc3(3))
            if(t_ixoverlay /= 0 .and. t_loc2(1) == 1)then
                max_loc3(1) = max(max_loc3(1), t_im)
            endif
        enddo
        ncount = max_ixl * max_iyl * max_kb
        allocate(t_v3(ncount), tt_v3(ncount), tt_loc3(max_loc3(1), max_loc3(2), max_loc3(3)))
    endif


    index = 1
    do 123 t_kb =  1, kb
    do 123 t_js = js, je
    do 123 t_is = is, ie
        t_v3(index) = v3(t_is, t_js, t_kb)
        index = index + 1
    123 continue
    if(ixoverlay /= 0 .and. loc2(1) == 1)then
        index = 1
        do 124 t_kb =  1, kb
        do 124 t_js = js, je
        do 124 t_is = is, is+ixoverlay
            tt_v3(index) = v3(t_is, t_js, t_kb)
            index = index + 1
        124 continue
    endif

    if (myid /= 0) then
        ! call mpi_send(ixoverlay, 1, mpi_integer, 0, 298, mpi_comm_ympi, ierr)
        ! call mpi_send(im, 1, mpi_integer, 0, 299, mpi_comm_ympi, ierr)
        ! call mpi_send(ixl, 1, mpi_integer, 0, 300, mpi_comm_ympi, ierr)
        ! call mpi_send(iyl, 1, mpi_integer, 0, 301, mpi_comm_ympi, ierr)
        ! call mpi_send(kb , 1, mpi_integer, 0, 302, mpi_comm_ympi, ierr)
        ! call mpi_send(is , 1, mpi_integer, 0, 303, mpi_comm_ympi, ierr)
        ! call mpi_send(js , 1, mpi_integer, 0, 304, mpi_comm_ympi, ierr)
        ! call mpi_send(ie , 1, mpi_integer, 0, 305, mpi_comm_ympi, ierr)
        ! call mpi_send(je , 1, mpi_integer, 0, 306, mpi_comm_ympi, ierr)
        ! call mpi_send(loc2(1), 2, mpi_integer, 0, 307, mpi_comm_ympi, ierr)
        ! call mpi_send(loc3(1), 3, mpi_integer, 0, 308, mpi_comm_ympi, ierr)
        ncount = (ie - is + 1) * (je - js + 1) * kb
        call mpi_send(t_v3(1), ncount, mpi_real, 0, 309, mpi_comm_ympi, ierr)
        if(ixoverlay /= 0 .and. loc2(1) == 1)then
            ncount = (ixoverlay + 1) * (je - js + 1) * kb
            call mpi_send(tt_v3(1), ncount, mpi_real, 0, 310, mpi_comm_ympi, ierr)
        endif
    else
        tt_is = 0
        tt_js = 0
        tt_kb = 0
        ncount = (ie - is + 1) * (je - js + 1) * kb
        do index = 1, ncount
            tt_loc3(loc3(1) + tt_is, loc3(2) + tt_js, loc3(3) + tt_kb) = t_v3(index)
            tt_is = tt_is + 1
            if (tt_is > ie - is) then
                tt_is = 0
                tt_js = tt_js + 1
            endif
            if (tt_js > je - js) then
                tt_js = 0
                tt_kb = tt_kb + 1
            endif
        enddo
        if(ixoverlay /= 0 .and. loc2(1) == 1)then
            tt_is = 0
            tt_js = 0
            tt_kb = 0
            ncount = (ixoverlay + 1) * (je - js + 1) * kb
            do index = 1, ncount
                tt_loc3(im - ixoverlay + tt_is, loc3(2) + tt_js, loc3(3) + tt_kb) = tt_v3(index)
                tt_is = tt_is + 1
                if (tt_is > ixoverlay) then
                    tt_is = 0
                    tt_js = tt_js + 1
                endif
                if (tt_js > je - js) then
                    tt_js = 0
                    tt_kb = tt_kb + 1
                endif
            enddo
        endif
        do pex = 1, npe - 1
            ! call mpi_recv(t_ixoverlay, 1, mpi_integer, pex, 298, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_im, 1, mpi_integer, pex, 299, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_ixl, 1, mpi_integer, pex, 300, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_iyl, 1, mpi_integer, pex, 301, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_kb , 1, mpi_integer, pex, 302, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_is , 1, mpi_integer, pex, 303, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_js , 1, mpi_integer, pex, 304, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_ie , 1, mpi_integer, pex, 305, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_je , 1, mpi_integer, pex, 306, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_loc2(1), 2, mpi_integer, pex, 307, mpi_comm_ympi, stat, ierr)
            ! call mpi_recv(t_loc3(1), 3, mpi_integer, pex, 308, mpi_comm_ympi, stat, ierr)
            t_ixoverlay = ixoverlayss(pex) 
            t_im        = imss(pex)        
            t_ixl       = ixlxx(pex)       
            t_iyl       = iylxx(pex)       
            t_kb        = kbss(pex)        
            t_is        = isss(pex)        
            t_js        = jsss(pex)        
            t_ie        = iess(pex)        
            t_je        = jess(pex)        
            t_loc2(:)   = t_loc2ss(pex, :) 
            t_loc3(:)   = t_loc3ss(pex, :) 
            ncount = (t_ie - t_is + 1) * (t_je - t_js + 1) * t_kb
            call mpi_recv(t_v3(1), ncount, mpi_real, pex, 309, mpi_comm_ympi, stat, ierr)
            if(t_ixoverlay /= 0 .and. t_loc2(1) == 1)then
                ncount = (t_ixoverlay + 1) * (t_je - t_js + 1) * t_kb
                call mpi_recv(tt_v3(1), ncount, mpi_real, pex, 310, mpi_comm_ympi, stat, ierr)
            endif
            tt_is = 0
            tt_js = 0
            tt_kb = 0
            ncount = (t_ie - t_is + 1) * (t_je - t_js + 1) * t_kb
            do index = 1, ncount
                tt_loc3(t_loc3(1) + tt_is, t_loc3(2) + tt_js, t_loc3(3) + tt_kb) = t_v3(index)
                tt_is = tt_is + 1
                if (tt_is > t_ie - t_is) then
                    tt_is = 0
                    tt_js = tt_js + 1
                endif
                if (tt_js > t_je - t_js) then
                    tt_js = 0
                    tt_kb = tt_kb + 1
                endif
            enddo
            if(t_ixoverlay /= 0 .and. t_loc2(1) == 1)then
                tt_is = 0
                tt_js = 0
                tt_kb = 0
                ncount = (t_ixoverlay + 1) * (t_je - t_js + 1) * t_kb
                do index = 1, ncount
                    tt_loc3(t_im - t_ixoverlay + tt_is, t_loc3(2) + tt_js, t_loc3(3) + tt_kb) = tt_v3(index)
                    tt_is = tt_is + 1
                    if (tt_is > t_ixoverlay) then
                        tt_is = 0
                        tt_js = tt_js + 1
                    endif
                    if (tt_js > t_je - t_js) then
                        tt_js = 0
                        tt_kb = tt_kb + 1
                    endif
                enddo
            endif
        enddo
        call open_nc(ncid, filename, action)
        call writenc(ncid, args    , tt_loc3(:, :, :), locs = [1, 1, 1])
        call close_nc(ncid)
    endif
end subroutine gather_out