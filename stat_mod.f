!gfortran, gcc version 5.4.0 20160609

!Программа статистического моделирования
program stat_mod

        !Ввод истодини данных      
        real lambda(50) , nu(50), mg(20), kg, ns, ns1, lambd, nui
        integer*2 mdi(50)
        integer*2 iy
        logical md(50,50),c(50),ch,cp,cq
        dimension t(50),tau(50),pb(50),nel(50),a(10)
        write(*,*) 'Исходные данные'
        write(*,*) '---------------'
        write(*,*) 'Введите допустимую ошибку'
        read(*,*) dd
        write(*,100) dd
100     format(1x, 'Допустимая ошибка ', f6.4)
        write(*,*) 'Введите доверительную вероятность'
        read(*,*) delta
        write(*,101) delta
101     format(1x, 'Доверительная вероятность ', f6.4)
        write(*,*) 'Введите проолжительность режима'
        read(*,*) ta
        write(*,102) ta
102     format(1x, 'Продолжительность режима', f6.0)    
        !Вычисление функции Z по данным б   
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
103     format(1x, 'Потребляемая в режиме мощность ', f7.0)
        write(*,*) 'Введите массив источников'      
        read(*,*) (mg(i),i=1,ki)
        write(*,*) 'Массив источников:'
        do 105 iii=1, ki
        write(*,106) iii, mg(iii)    
106     format(1x, 'Мощность ',i2, '-го источника ',f7.0)
105     continue
        write(*,*) 'Ввод лямбда, ну, Рв'
        write(*,*) '---------------'
1105    jn=0
1102    continue
        write(*,*) 'Введите номер элемента'
        read(*,*) neli
        if(neli.eq.0) goto 1102
        jn=jn+1
        nel(jn)=neli
        goto 1101
1101    continue
        write(*,*) 'Введите: лямбда, ну, Рв'
        read(*,*) lambd, nui,pbi
        do 1103 jj=1,jn
        kkk=nel(jj)
        lambda(kkk)=lambd
        nu(kkk)=nui
1103    pb(kkk)=pbi
        sum=sum+jn
        if(sum.lt.m) goto 1105
        write(*,*) 'Диагностическая матрица'
        write(*,*) 'Вввод Диагностической матрицы:'
        write(*,*) '---------------'    
        do 10 i=1, kcm
        write(*,110) i
110     format(1x, 'Введите номера столбцов' ,i2,'-ой строки')
        read(*,*) nn1,nn2
        do 2001 ne=1,m
2001    mdi(ne)=0
        mdi(nn1)=1
        mdi(nn2)=1
        write(*,7) (mdi(j),j=1,m)
7       format(1x,79i1)
        do 10 j=1,m
        if(mdi(j)) 8,8,9
8       md(i,j)=.false.
        goto 10
9       md(i,j)=.true.
10      continue
        write(*,*) 'При определении источников'
        write(*,*) 'проверять состояние их щитов?'
        write(*,*) 'Да - введите 1'
        write(*,*) 'Нет - введите 0'
        read(*,*) otv
        write(*,*) '-----------------------------------'
        write(*,*) '|  №  |  лямбда |    ню   |  Рв   |'
        write(*,*) '-----------------------------------'      
        do 400 j=1, m
        amb=lambda(j)*1e6
        u1=nu(j)*1e3
        write(*, 401) j, amb,u1,pb(j)
401     format(1x, '| ',i3, ' | ', f7.2, ' | ',f7.2,' | ',f5.3, ' |')
400     continue
        write(*,*) '-----------------------------------' 
        write(*,*) '     РЕЗУЛЬТАТЫ МОДЕЛИРОВАНИЯ'       
        write(*,*) '---------------------------------------'
        write(*,*) '| Д |  дельта | n исп. |  кг  |    В  |'
        write(*,*) '---------------------------------------'
  
        iy=0
        ts=0    ! Установка начальных значений
        ns=100  !времени, количества спттавки
        sqt=0.  ! и кол-ва отказов
        l=0.
12      continue
        do 32 n=1, 100 !цикл испытаний (походов судна)
            do 13 i=1, m !цикл перебора элементов системы
				c(i)=.true.
13          continue
            t(i)=0.
            th=0
            tsl=0
            cp=.true.
            cq=.true.
            do 14 i=1,m ! Определение времени
3001			continue
				call random_number(z) ! до наступления
				if(z.le.0) goto 3001 !отказов каждого
				tau(i)=-alog(z)/lambda(i) !из элементов системы.
14          continue
			t(i)=tau(i)     !определение минимального интервала времени
15          continue
			tm=t(1) ! (ts тут или нет хз) отказов элементов и номера 
            k=1 !элемента отказавшего первым
            do 16 i=2, m
				if(t(i).ge.tm) goto 16
				tm=t(i)
				k=i
16          continue
            if(tm.le.ta) goto 37
				if(cp) goto 31
				ti=ta-th  !счет времени неработоспособного состяния системы
				ts=ts+ti
				tsl=tsl+ti
				goto 31
37          continue
			if(c(k)) goto 38
				c(k)=.true.    !если у элемента к произошло восстановление, моделируется интервал времени до очередного отказа        
				call random_number(z)
				tau(k)=-alog(z)/lambda(k)
				go to 21
38          continue
			c(k)=.false. ! если у элемента произошел отказ, определяется возможность его(отказа) устранения
			call random_number(z)
            
            !write(*,*) 'd1:',z, pb(k)
            if(z.ge.pb(k)) goto 20
				call random_number(z)
                    !write(*,*) 'd2:',z, pb(k)
					tau(k)=-alog(z)/nu(k)
				goto 21
20          continue
			tau(k)=ta
21			continue          
			t(k)=t(k)+tau(k)
            if(c(k).and.cp) goto 22
				go to 23
22			continue          
			ch=.true.
            go to 28
23			continue          
			do 25 j=1, kcm
				do 24 i=1,m
					if(.not.c(i).or..not.md(j,i)) go to 24 !определение состояния системы по результату сравненения ее состояния со строками диагностической матрицы
						ch=.true.
						go to 25
24          	continue
				ch=.false.
				go to 28
25          continue
26          continue
			w=0
            do 27 i=1, ki
                if(otv) 2000,2000,3000
3000		        continue        
                    if(c(i).and.c(ki+i)) w=w+mg(i) !Суммирование мощности исправных источников
                    goto 27
2000		        continue        
                    if(c(i)) w=w+mg(i)
27          continue
            ch=.true.
            if (w.lt.wd) ch=.false.
28          continue
			if (cp.or..not.ch) goto 29 !Опеределение интервала времени неработоспособного состояния
				ti=tm-th
				ts=ts+ti
				tsl=tsl+ti
29          continue
			if(.not.cp.or.ch) goto 30
				th=tm       ! Фиксируется время неработоспособного состояния
				cq=.false.
30			continue          
			cp=ch 
            go to 15
31			continue          
			sqt=sqt+tsl*tsl 
            if(.not.cq) l=l+1   
32      continue
        kg=l-ts/ta/ns  !Вычисление коэф. готовности системы
        r=1.-l*1./ns  ! Вычисление вероятности безотказной работы 
        d=zdelta*sqrt(r*(l-r)/ns)  !значение доверительной ошибки
        if(l.gt.0.and.d.le.dd) go to 34
            write(*,33) d, delta, ns, kg, r
33          format(1x,'| ',f7.5,' | ',f5.3, ' | ' ,f6.0,' | ', f7.5,' | ', f7.5, ' |')
            ns=ns+100
            go to 12
34      continue
        write(*,33) d, delta, ns, kg, r
        write(*,*) '---------------------------------------'   
        stop
end program stat_mod
