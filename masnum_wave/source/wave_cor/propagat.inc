!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: propagat

subroutine propagat

!-------------------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------------------

integer :: j, k, i, iyyz, ixx, jyy, ixx1, jyy1, iwk, iwk1, i1, jth, jth1
real :: x0, y0, d0, dddx0, dddy0, duxdx0, duxdy0, duydy0, duydx0, th0
!  real :: sinth, costh, deltt, wk0, ws0, dk0, cg, cgx, cgy, delts, wecs
real :: sinth, costh, wk0, ws0, dk0, cg, cgx, cgy 
  !angx_d, angy_d, , delts, radx_d, rady_d, wecs
real :: xx, yy, x1, x2, y1, y2, dsidd
real :: ssr1, ssr2, ssrwk, ssrth, wks, fien, fien1, wk1, wk2, dtth0, ths
real :: e11, e12, e13, e14, e1, e21, e22, e23, e24, e2
real :: e32, e33, e34, e3, e41, e42, e43, e44, e4, e31
real :: th1, th2, exxyy, x_d, y_d

real :: dx, dy
  
integer :: iii

!-------------------------------------------------------------------------------
!--------------------------------update-----------------------------------------
! ======== 2017.2.15 15:33 ==========
integer :: iwkleft, iwkrigh, iwkmid
real :: r_tmp, t_tand_rs
real :: sin_th(jl), cos_th(jl), tand_ic(iys:iyl), tand_rs(jl, iys:iyl)
real :: dsidd_x0(kl), dsidd_y0(kl)
real :: yy_t(kl), xx_t(kl)
real :: deltt_rs, deltt_ic, deltt_rs_ic(iys:iyl)
real :: costh_x0, costh_y0, sinth_x0, sinth_y0
real :: sin_rs(jl), cos_rs(jl), sin_rs_ic(jl, iys:iyl), cos_rs_ic(jl, iys:iyl)
real :: cos_rs_k(kl, jl), sin_rs_ic_k(kl, jl, iys:iyl)

deltt_rs = deltt/rs*180./pi

do ic = iys, iyl
    deltt_rs_ic(ic) = deltt/rslat(ic)*180./pi
    tand_ic(ic) = tand(y(ic)) / rs
enddo

do  j=1,jl
    th0=thet(j)
    sin_th(j)=sin(th0)
    cos_th(j)=cos(th0)
    sin_rs(j) = sin_th(j) / rs
    cos_rs(j) = cos_th(j) / rs
    do ic = iys, iyl
        tand_rs(j, ic) = cos_th(j) * tand_ic(ic)
        sin_rs_ic(j, ic) = sin_th(j) / rslat(ic)
        cos_rs_ic(j, ic) = cos_th(j) / rslat(ic)
        do k = 1, kl
            cos_rs_k(k, j) = cos_rs(j) / wk(k)
            sin_rs_ic_k(k, j, ic) = sin_rs_ic(j, ic) / wk(k)
        enddo
    enddo
enddo
! ===================================
!-------------------------------------------------------------------------------

do 100 ic=iys,iyl
do 100 ia=ixs,ixl

    if(nsp(ia,ic).ne.1) cycle

    x0=x(ia)
    y0=y(ic)
    d0=d(ia,ic)

    dx = deltx(ia)
    dy = delty(ic)

!    dddx0=(d(ia+1,ic)-d(ia-1,ic))/(2*dx)
!    dddy0=(d(ia,ic+1)-d(ia,ic-1))/(2*dy)
!yinxq:    if(glbflag == 0)then ! --- for global mdel. 2010.04.09
!yinxq:      if(ia == ixs)then
!yinxq:        dddx0=(d(ia+1,ic)-d(ixl-1,ic))/(2.*dx)
!yinxq:      elseif(ia == ixl)then
!yinxq:        dddx0=(d(ixs+1,ic)-d(ia-1,ic))/(2.*dx)
!yinxq:      else
!yinxq:        dddx0=(d(ia+1,ic)-d(ia-1,ic))/(2.*dx)
!yinxq:      endif
!yinxq:    else
!yinxq:      if(ia == ixs)then
!yinxq:        dddx0=(d(ia+1,ic)-d(ia,ic))/(dx)
!yinxq:      elseif(ia == ixl)then
!yinxq:        dddx0=(d(ia,ic)-d(ia-1,ic))/(dx)
!yinxq:      else
!yinxq:        dddx0=(d(ia+1,ic)-d(ia-1,ic))/(2.*dx)
!yinxq:      endif
!yinxq:    endif
!yinxq:    dddy0=(d(ia,ic+1)-d(ia,ic-1))/(2.*dy)

    dddx0 = dddx(ia, ic) ! yinxq
    dddy0 = dddy(ia, ic) ! yinxq

! ========== 2017.2.14 15:54 ==========
    duxdx0=uxx(ia,ic)
    duxdy0=uxy(ia,ic)
    duydx0=uyx(ia,ic)
    duydy0=uyy(ia,ic)

    deltt_ic = deltt_rs_ic(ic)

    do k = 1, kl
        wk0=wk(k)
        ws0=zpi*wf(k,ia,ic)
        dk0=d0*wk0
        r_tmp = cosh(dk0)
        if (dk0.lt.40.) then
            dsidd=0.5 * g / r_tmp * wk0**2 / ws0 / r_tmp
        else
            dsidd=0.
        endif
        dsidd_x0(k) = dsidd * dddx0 * 180. / pi
        dsidd_y0(k) = dsidd * dddy0 * 180. / pi
    enddo
! ===============================

    do 200 j=1,jl
        th0=thet(j)
        sinth=sin_th(j)
        costh=cos_th(j)

! ========== 2017.2.14 15:54 ==========
        t_tand_rs = tand_rs(j, ic)
        costh_x0 = costh * duxdx0 * 180. / pi
        costh_y0 = costh * duxdy0 * 180. / pi
        sinth_x0 = sinth * duydx0 * 180. / pi
        sinth_y0 = sinth * duydy0 * 180. / pi
! =====================================

        do 500 k=1,kl
!yinxq:        deltt=delttm*60.
            wk0=wk(k)
            cg=ccg(k,ia,ic)
            cgx=cg*costh
            cgy=cg*sinth

        !******  1.  "the calculation of wave engery-current spreading"

        !      wecs=sqrt((cgx+ux)**2+(cgy+uy)**2)
!        wecs=sqrt((cgx+ux(ia,ic))**2+(cgy+uy(ia,ic))**2)
!        delts=wecs*deltt
!        x_d=delts*costh
!        y_d=delts*sinth
!        radx_d=x_d/rslat(ic)
!        rady_d=y_d/rs
!        angx_d=radx_d*180./pi
!        !angy_d=y_d/rs*180./pi
!        angy_d=rady_d*180./pi
!        delts = sqrt((cgx+ux(ia,ic))**2+(cgy+uy(ia,ic))**2)*deltt*180./pi
!        angx_d = delts*costh/rslat(ic)
!        angy_d = delts*sinth/rs
!        if(abs(angx_d).gt.dx.or.abs(angy_d).gt.dy) then
!          write(6, *)angx_d, dx, ia, ic, x(ia), y(ic), d(ia,ic), &
!          sqrt((cgx+ux(ia,ic))**2+(cgy+uy(ia,ic))**2), cg
!          write(6,*) 'be careful !!! program have to stop !!!'
!          stop
!        endif

        !----------------------------------
! ======================= 2017.2.13 23:15 =======================
            ! xx=x0-deltt*(cgx+ux(ia,ic))/rslat(ic)*180./pi
            ! yy=y0-deltt*(cgy+uy(ia,ic))/rs*180./pi
            xx=x0-deltt_ic*(cgx+ux(ia,ic))
            yy=y0-deltt_rs*(cgy+uy(ia,ic))
! =====================================================================

!        if(xx.lt.0.)xx=xx+360.
!        if(xx.ge.360.)xx=xx-360.
!        if(yy.lt.0.)yy=yy+180.
!        if(yy.ge.180.)yy=yy-180.
!        ixx=int(xx/dx)+1
!        jyy=int(yy/dy)+1
!        ixx=int(xx/dx)+1

!-------------------------------------------------------------------------------

            if(glbflag == 0)then
                if(xx < x(ixs))xx = xx + 360.
                if(xx > x(ixl))xx = xx - 360.
            else
                if(xx.lt.x(ixs))xx=x(ixs)
                if(xx.ge.x(ixl))xx=x(ixl)
            endif
            ixx = 0

!----------------------------------use binary search----------------------------
          !   tmp = ixx
          !   do iii = ixs, ixl-1
          !       if(xx >= x(iii) .and. xx <= x(iii+1))then
          !           ixx = iii; exit
          !       endif
          !   enddo
! 2017.2.10 22:00
            iwkleft = ixs
            iwkrigh = ixl - 1
            if (iwkleft <= iwkrigh) then
                do iii = 1, 999
                    if (iwkleft > iwkrigh) exit
                    iwkmid = shiftr(iwkrigh + iwkleft, 1)
                    if(xx <= x(iwkmid)) then
                        iwkrigh = iwkmid - 1
                    else
                        iwkleft = iwkmid + 1
                    endif
                enddo
                ixx = iwkleft - 1
                if (ixx < ixs) ixx = ixs
                if (ixx > ixl - 1) ixx = ixl - 1
            endif
!-------------------------------------------------------------------------------

            if(ixx < 1 .or. ixx > ixl) cycle
            
            if(yy < y(iys)) yy = y(iys)
            if(yy > y(iyl)) yy = y(iyl)
            jyy = 0

!-----------------------------------use binary search---------------------------
          !  tmp = jyy
          !  do iii = iys, iyl-1
          !       if(yy >= y(iii) .and. yy <= y(iii+1))then
          !           jyy = iii; 
          !           exit
          !       endif
          !   enddo
! 2017.2.10 22:07
            iwkleft = iys
            iwkrigh = iyl - 1
            if (iwkleft <= iwkrigh) then
                do iii = 1, 999
                    if (iwkleft > iwkrigh) exit
                    iwkmid = shiftr(iwkrigh + iwkleft, 1)
                    if(yy <= y(iwkmid)) then
                        iwkrigh = iwkmid - 1
                    else
                        iwkleft = iwkmid + 1
                    endif
                enddo
                jyy = iwkleft - 1
                if (jyy < iys) jyy = iys
                if (jyy > iyl - 1) jyy = iyl - 1
            endif
!------------------------------------------------------------------------------- 

            if(jyy < 0 .or. jyy > iyl)cycle

!-------------------------------------------------------------------------------

            ixx1=ixx+1
            jyy1=jyy+1

            x1=x(ixx)
            x2=x(ixx1)
            y1=y(jyy)
            y2=y(jyy1)

            !******  2.  "the effect of refraction caused by topography and current"
            !      ddds= dddx0*costh+dddy0*sinth
            !      dddn=-dddx0*sinth+dddy0*costh

            !      duxds= duxdx0*costh+duxdy0*sinth
            !      duyds= duydx0*costh+duydy0*sinth

            !      duxdn=-duxdx0*sinth+duxdy0*costh
            !      duydn=-duydx0*sinth+duydy0*costh

            !      dudsl=duxds*costh+duyds*sinth
            !      dudnl=duxdn*costh+duydn*sinth
            !      ssrwk=-(dsidd*ddds+wk0*dudsl)
            !      ssrth=-(dsidd*dddn+wk0*dudnl)/wk0

!yinxq 2012/9/26 15:30:26        ssr1=(dsidd*dddx0+wk0*costh*duxdx0+wk0*sinth*duydx0)
!yinxq 2012/9/26 15:30:31        ssr2=(dsidd*dddy0+wk0*costh*duxdy0+wk0*sinth*duydy0)
! ================================ 2017.2.15 15:35 ================================
            ! ssr1=(dsidd*dddx0+wk0*costh*duxdx0+wk0*sinth*duydx0)*180./pi
            ! ssr2=(dsidd*dddy0+wk0*costh*duxdy0+wk0*sinth*duydy0)*180./pi
            ! ssrwk=-(ssr1*costh/rslat(ic)+ssr2*sinth/rs)
            ! ssrth=(ssr1*sinth/rslat(ic)-ssr2*costh/rs)/wk0
            ! ssrth=ssrth-cg*costh*tand(y(ic))/rs !?yinxq

            ssr1 = dsidd_x0(k) + wk0 * costh_x0 + wk0 * sinth_x0
            ssr2 = dsidd_y0(k) + wk0 * costh_y0 + wk0 * sinth_y0
            ssrwk=-(ssr1*cos_rs_ic(j, ic)      + ssr2 * sin_rs(j))
            ssrth= (ssr1*sin_rs_ic_k(k, j, ic) - ssr2 * cos_rs_k(k, j)) - cg * t_tand_rs
! =================================================================================

            wks=wk0-deltt*ssrwk
            if (wks.lt.0.) wks=0.
            if (wks.le.wkmin) then
                iwk=1
                iwk1=1
                fien=0.
                fien1=1.
                wk1=0.
                wk2=wk(iwk1)
            else
                if (wks.lt.wk(kld)) then

                !======= THE FOLLOWING BELONGS TO THE LAGFD-WAM WAVE MODEL
                !            iwk=int(log(wks/wkmin)/log(pwk))+1
                !==========================================================
                !            if (wks.lt.wk(iwk)) iwk=iwk-1
                !            if (wks.gt.wk(iwk+1)) iwk=iwk+1

                !========   The following is modified by Yang Yongzeng =====c

!------------------------use binary search------------------------
             !   tmp = iwk
             !   do iyyz=1,kld
             !     if(wks.ge.wk(iyyz).and.wks.lt.wk(iyyz+1)) iwk=iyyz
             !   enddo
! 2017.2.10 22.11
                    iwkleft = 2
                    iwkrigh = kld + 1
                    if (iwkleft <= iwkrigh) then
                        do iii = 1, 999
                            if (iwkleft > iwkrigh) exit
                            iwkmid = shiftr(iwkrigh + iwkleft, 1)
                            if (wk(iwkmid) > wks) then
                                iwkrigh = iwkmid - 1
                            else 
                                iwkleft = iwkmid + 1
                            endif
                        enddo
                        iwk = iwkleft - 1
                    endif
!---------------------------------------------------------------------
                    iwk1=iwk+1
!=====================================================================c

                    wk1=wk(iwk)
                    wk2=wk(iwk+1)
                    if (iwk.lt.kl) then
                        iwk1=iwk+1
                        fien=1.
                        fien1=1.
                    else
                        i=iwk-kl+1
                        i1=i+1
                         fien=wkh(i)
                        fien1=wkh(i1)
                         iwk=kl
                        iwk1=kl
                    endif
                else
                    wks=wk(kld)
                    iwk=kl
                    iwk1=kl
                    i=kld-kl+1
                    i1=i+1
                    fien=wkh(i)
                    fien1=wkh(i1)
                endif
            endif

            dtth0=deltt*ssrth
            !      if (abs(dtth0).gt.pi2) then
            !         dtth0ab=abs(dtth0)
            !         sig=dtth0/dtth0ab
            !         dtth0l=amax1(dtth0ab,pi2)
            !         dtth0=sig*dtth0l
            !      endif
            ths=th0-dtth0
            if (ths.ge.(-1.0)*zpi.and.ths.lt.(2.0)*zpi) then
                goto 1200
            else
                ths=th0

!yinxq 2012/9/26 15:37:56          write(6,*) ia,ic,k,j
!yinxq 2012/9/26 15:38:01          write(6,*) ths,th0,zpi
!yinxq 2012/9/26 15:38:07          write(6,*) 'The wave is spinning ,I can not computer it !!!'
!yinxq 2012/9/26 15:38:12          stop

            endif
            1200      continue
            if (ths.ge.zpi) ths=ths-zpi
            if (ths.lt. 0.) ths=ths+zpi
            jth=int(ths/deltth)+1
            jth1=jth+1
            if(jth.eq.jlp1)then
                jth=jl
                jth1=1
            endif
            if (jth1.eq.jlp1) jth1=1

            th1=thet(jth)
            th2=thet(jth+1)

            !******  3.  "determing the wave energy at the physical space point (xx,yy)
            !******                  and the wave space point (wks,ths)"

            !*       3.1 "at physical space point (x1,y1)"

            e11=ee(iwk ,jth ,ixx,jyy)*fien
            e12=ee(iwk1,jth ,ixx,jyy)*fien1
            e13=ee(iwk ,jth1,ixx,jyy)*fien
            e14=ee(iwk1,jth1,ixx,jyy)*fien1

            call inter(wk1,wk2,th1,th2,e11,e12,e13,e14,wks,ths,e1)

            !*       3.2 "at physical space point (x2,y1)"

            e21=ee(iwk ,jth ,ixx1,jyy)*fien
            e22=ee(iwk1,jth ,ixx1,jyy)*fien1
            e23=ee(iwk ,jth1,ixx1,jyy)*fien
            e24=ee(iwk1,jth1,ixx1,jyy)*fien1

            call inter(wk1,wk2,th1,th2,e21,e22,e23,e24,wks,ths,e2)

            !*       3.3 "at physical space point (x1,y2)"

            e31=ee(iwk ,jth ,ixx,jyy1)*fien
            e32=ee(iwk1,jth ,ixx,jyy1)*fien1
            e33=ee(iwk ,jth1,ixx,jyy1)*fien
            e34=ee(iwk1,jth1,ixx,jyy1)*fien1

            call inter(wk1,wk2,th1,th2,e31,e32,e33,e34,wks,ths,e3)

            !*       3.4 "at physical space point (x2,y2)"

            e41=ee(iwk ,jth ,ixx1,jyy1)*fien
            e42=ee(iwk1,jth ,ixx1,jyy1)*fien1
            e43=ee(iwk ,jth1,ixx1,jyy1)*fien
            e44=ee(iwk1,jth1,ixx1,jyy1)*fien1

            call inter(wk1,wk2,th1,th2,e41,e42,e43,e44,wks,ths,e4)

            !*       3.5 "at physical space point (xx,yy)"

            call inter(x1,x2,y1,y2,e1,e2,e3,e4,xx,yy,exxyy)

    !        e(k,j,ia,ic)=exxyy
            e(k,j,ia,ic)=max(exxyy, small)
            !        write(6,*) e(k,j,ia,ic)

        500      continue
    200      continue
100      continue

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine propagat

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
