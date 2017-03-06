!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mean2
 ! ¼ÆËãÆ×¾Ø£ºÓÃÓÚ¼ÆËãÔ´º¯Êý

subroutine mean2

implicit none

!-------------------------------------------------------------------------------

integer :: k, k1, i, i1, j
real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1

!-------------------------------------------------------------------------------

! ============== 2017.2.13 18:34 ==============
real :: t_ae, t_asi, t_ark, t_awf, t_awk
real :: sqr_wkk(kld), sqr_wkk1(kld)
do k = 1, kld
    k1=k+1
    wkk=wk(k)
    wkk1=wk(k1)
    sqr_wkk(k) = sqrt(wkk)
    sqr_wkk1(k) = sqrt(wkk1)
enddo
! =============================================

do 100 ic=iys,iyl
do 100 ia=ixs,ixl

    ae(ia,ic)=0.
    asi(ia,ic)=0.
    awf(ia,ic)=0.
    awk(ia,ic)=0.
    ark(ia,ic)=0.

    if(nsp(ia,ic).ne.1) cycle

! ============== 2017.2.13 18:33 ==============
    do 202 k=1,kld
        k1=k+1
        i=k-kl+1
        i1=i+1
        dwkk=dwk(k)
        wfk=wf(k,ia,ic)
        wfk1=wf(k1,ia,ic)
        !      wfk=fr(k)
        !      wfk1=fr(k1)
        wsk=zpi*wfk
        wsk1=zpi*wfk1
        wkk=wk(k)
        wkk1=wk(k1)

! ================= 2017.2.13 20:44 =================
        t_ae = 0.
        t_asi = 0.
        t_ark = 0.
        t_awf = 0.
        t_awk = 0.
        ekj = 1.
        ekj1 = 1.
        
        ! if (k.ge.kl) then
        !     ekj=wkh(i)
        !     ekj1=wkh(i1)
        ! endif

        do 200 j=1,jl
            if (k.lt.kl) then
                ekj=e(k,j,ia,ic)
                ekj1=e(k1,j,ia,ic)
            else
                ekj=e(kl,j,ia,ic)*wkh(i)
                ekj1=e(kl,j,ia,ic)*wkh(i1)
            endif

            ! ekj=ekj * e(k,j,ia,ic)
            ! ekj1=ekj1 * e(k1,j,ia,ic)

            ! ae(ia,ic)=ae(ia,ic)+(ekj+ekj1)*dwkk
            ! awf(ia,ic)=awf(ia,ic)+(ekj*wfk+ekj1*wfk1)*dwkk
            ! asi(ia,ic)=asi(ia,ic)+(ekj/wsk+ekj1/wsk1)*dwkk
            ! awk(ia,ic)=awk(ia,ic)+(ekj*wkk+ekj1*wkk1)*dwkk
            ! ark(ia,ic)=ark(ia,ic)+(ekj/sqrt(wkk)+ekj1/sqrt(wkk1))*dwkk

            t_ae=t_ae+ekj+ekj1
            t_awf=t_awf+ekj*wfk+ekj1*wfk1
            t_asi=t_asi+ekj/wsk+ekj1/wsk1
            t_awk=t_awk+ekj*wkk+ekj1*wkk1
            t_ark=t_ark+ekj/sqr_wkk(k)+ekj1/sqr_wkk1(k)
        200 continue
        ae(ia,ic) = ae(ia,ic) + t_ae * dwkk
        awf(ia, ic) = awf(ia, ic) + t_awf * dwkk
        asi(ia, ic) = asi(ia, ic) + t_asi * dwkk
        awk(ia, ic) = awk(ia, ic) + t_awk * dwkk
        ark(ia, ic) = ark(ia, ic) + t_ark * dwkk
! =================================================================
    202 continue

    asi(ia,ic)=ae(ia,ic)/asi(ia,ic)
    awf(ia,ic)=awf(ia,ic)/ae(ia,ic)
    awk(ia,ic)=awk(ia,ic)/ae(ia,ic)
    ark(ia,ic)=(ae(ia,ic)/ark(ia,ic))**2

  100      continue

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine mean2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
