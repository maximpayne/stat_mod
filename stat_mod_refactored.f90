program stat_mod
    implicit none
    integer, parameter :: MAX_M = 50, MAX_K = 50, MAX_G=20
    real :: lambda(MAX_M), nu(MAX_M), mg(MAX_G)
    real :: kg, ns, ns1, lambd, nui
    integer :: mdi(MAX_M)
    logical :: md(MAX_K,MAX_M), c(MAX_M), ch, cp, cq
    real :: t(MAX_M), tau(MAX_M), pb(MAX_M), nel(MAX_M)
    real :: a(10)
    integer :: i,j,n,nn1,nn2,m,kcm,ki,jn,neli
    real :: dd, delta, ta, zdelta, wd
    real :: amb,u1
    real :: ts,th,tm,ti,tsl
    real :: sqt, l
    integer :: k
    real :: r,d,w
    integer :: ns_int
    logical :: otv
    real :: z

    write(*,*) 'Исходные данные'
    write(*,*) '---------------'
    write(*,*) 'Введите допустимую ошибку'
    read(*,*) dd
    write(*,100) dd
100 format(1x, 'Допустимая ошибка ', f6.4)
    write(*,*) 'Введите доверительную вероятность'
    read(*,*) delta
    write(*,101) delta
101 format(1x, 'Доверительная вероятность ', f6.4)
    write(*,*) 'Введите проолжительность режима'
    read(*,*) ta
    write(*,102) ta
102 format(1x, 'Продолжительность режима', f6.0)
    zdelta=2569.449*delta**3-7196.125*delta**2+6722.0792*delta-2092.4842
    write(*,*) 'Введите количество элементов системы'
    read(*,*) m
    write(*,*) 'Количество элементов системы', m
    write(*,*) 'Bведите количество условий отказа системы'
    read(*,*) kcm
    write(*,*) 'Количество условий отказа системы ', kcm
    write(*,*) 'Введите количество установленных источников'
    read(*,*) ki
    write(*,*) 'Количество установленных источников', ki
    write(*,*) 'Введите потребляемую в режиме мощность'
    read(*,*) wd
    write(*,103) wd
103 format(1x, 'Потребляемая в режиме мощность ', f7.0)
    write(*,*) 'Введите массив источников'
    read(*,*) (mg(i),i=1,ki)
    write(*,*) 'Массив источников:'
    do i=1,ki
        write(*,106) i, mg(i)
    end do
106 format(1x, 'Мощность ',i2, '-го источника ',f7.0)

    jn=0
    do
        write(*,*) 'Введите номер элемента'
        read(*,*) neli
        if(neli==0) cycle
        jn=jn+1
        nel(jn)=neli
        write(*,*) 'Введите: лямбда, ну, Рв'
        read(*,*) lambd, nui, pb(jn)
        lambda(neli)=lambd
        nu(neli)=nui
        if(jn>=m) exit
    end do

    write(*,*) 'Диагностическая матрица'
    write(*,*) 'Вввод Диагностической матрицы:'
    write(*,*) '---------------'
    do i=1,kcm
        write(*,110) i
110     format(1x, 'Введите номера столбцов' ,i2,'-ой строки')
        read(*,*) nn1,nn2
        mdi(:)=0
        mdi(nn1)=1
        mdi(nn2)=1
        write(*,7) (mdi(j),j=1,m)
7       format(1x,79i1)
        do j=1,m
            md(i,j)=mdi(j)==1
        end do
    end do

    write(*,*) 'При определении источников'
    write(*,*) 'проверять состояние их щитов?'
    write(*,*) 'Да - введите 1'
    write(*,*) 'Нет - введите 0'
    read(*,*) otv

    write(*,*) '-----------------------------------'
    write(*,*) '|  №  |  лямбда |    ню   |  Рв   |'
    write(*,*) '-----------------------------------'
    do j=1,m
        amb=lambda(j)*1e6
        u1=nu(j)*1e3
        write(*,401) j, amb,u1,pb(j)
    end do
401 format(1x, '| ',i3, ' | ', f7.2, ' | ',f7.2,' | ',f5.3, ' |')
    write(*,*) '-----------------------------------'
    write(*,*) '     РЕЗУЛЬТАТЫ МОДЕЛИРОВАНИЯ'
    write(*,*) '---------------------------------------'
    write(*,*) '| Д |  дельта | n исп. |  кг  |    В  |'
    write(*,*) '---------------------------------------'

    ts=0.0
    ns_int=100
    sqt=0.0
    l=0.0

    do
        do n=1,100
        do i=1,m
            c(i)=.true.
            t(i)=0.0
        end do
        th=0.0
        tsl=0.0
        cp=.true.
        cq=.true.
        do i=1,m
            do
                call random_number(z)
                if(z>0.0) exit
            end do
            tau(i)=-log(z)/lambda(i)
            t(i)=tau(i)
        end do

        do
            tm=t(1)
            k=1
            do i=2,m
                if(t(i)<tm) then
                    tm=t(i)
                    k=i
                end if
            end do

            if(tm>ta) then
                if(.not.cp) then
                    ti=ta-th
                    ts=ts+ti
                    tsl=tsl+ti
                end if
                exit
            end if

            if(.not.c(k)) then
                c(k)=.true.
                call random_number(z)
                tau(k)=-log(z)/lambda(k)
            else
                c(k)=.false.
                call random_number(z)
                if(z<pb(k)) then
                    call random_number(z)
                    tau(k)=-log(z)/nu(k)
                else
                    tau(k)=ta
                end if
            end if
            t(k)=t(k)+tau(k)

            if(c(k).and.cp) then
                ch=.true.
            else
                ch=.false.
                do j=1,kcm
                    ch=.true.
                    do i=1,m
                        if(.not.c(i) .or. .not.md(j,i)) then
                            ch=.false.
                            exit
                        end if
                    end do
                    if(ch) exit
                end do
                w=0.0
                if(otv) then
                    do i=1,ki
                        if(c(i).and.c(ki+i)) w=w+mg(i)
                    end do
                else
                    do i=1,ki
                        if(c(i)) w=w+mg(i)
                    end do
                end if
                if(w<wd) ch=.false.
            end if

            if(cp.and..not.ch) then
                ti=tm-th
                ts=ts+ti
                tsl=tsl+ti
            else if(.not.cp.and.ch) then
                th=tm
                cq=.false.
            end if
            cp=ch
        end do
        sqt=sqt+tsl*tsl
        if(.not.cq) l=l+1
    end do

        kg=l-ts/ta/ns_int
        r=1.-l/ns_int
        d=zdelta*sqrt(r*(l-r)/ns_int)
        write(*,33) d, delta, ns_int, kg, r
        if(l>0 .and. d<=dd) exit
        ns_int=ns_int+100
    end do
33  format(1x,'| ',f7.5,' | ',f5.3, ' | ',f6.0,' | ', f7.5,' | ', f7.5,' |')
    write(*,*) '---------------------------------------'
end program stat_mod
