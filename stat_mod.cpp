#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include "parse.h"

inline std::string fortranF(double v, int width, int prec) {
    std::ostringstream ss;
    if(prec==0) ss << std::showpoint;
    ss << std::fixed << std::setprecision(prec) << v;
    std::string s = ss.str();
    if((int)s.size() > width) return std::string(width, '*');
    if((int)s.size() < width) s = std::string(width - s.size(),' ') + s;
    return s;
}

int main() {
    const int MAX_M = 50;
    const int MAX_K = 50;
    const int MAX_G = 20;

    std::vector<double> lambda(MAX_M+1, 0.0);
    std::vector<double> nu(MAX_M+1, 0.0);
    std::vector<double> mg(MAX_G+1, 0.0);
    std::vector<int> mdi(MAX_M+1, 0);
    std::vector<std::vector<bool>> md(MAX_K+1, std::vector<bool>(MAX_M+1,false));
    std::vector<bool> c(MAX_M+1, false);
    std::vector<double> t(MAX_M+1,0.0), tau(MAX_M+1,0.0), pb(MAX_M+1,0.0), nel(MAX_M+1,0.0);

    double dd, delta, ta, zdelta, wd;
    int m, kcm, ki;

    std::cout << "Исходные данные" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Введите допустимую ошибку" << std::endl;
    std::cin >> dd;
    std::cout << " Допустимая ошибка " << std::fixed << std::setprecision(4) << dd << std::endl;
    std::cout << "Введите доверительную вероятность" << std::endl;
    std::cin >> delta;
    std::cout << " Доверительная вероятность " << std::fixed << std::setprecision(4) << delta << std::endl;
    std::cout << "Введите проолжительность режима" << std::endl;
    std::cin >> ta;
    std::cout << " Продолжительность режима " << std::showpoint << std::fixed << std::setprecision(0) << ta << std::endl;
    std::cout << std::noshowpoint;

    zdelta = 2569.449*std::pow(delta,3) - 7196.125*delta*delta + 6722.0792*delta - 2092.4842;

    std::cout << "Введите количество элементов системы" << std::endl;
    std::cin >> m;
    std::cout << " Количество элементов системы " << m << std::endl;
    std::cout << "Bведите количество условий отказа системы" << std::endl;
    std::cin >> kcm;
    std::cout << " Количество условий отказа системы " << kcm << std::endl;
    std::cout << "Введите количество установленных источников" << std::endl;
    std::cin >> ki;
    std::cout << " Количество установленных источников " << ki << std::endl;
    std::cout << "Введите потребляемую в режиме мощность" << std::endl;
    std::cin >> wd;
    std::cout << " Потребляемая в режиме мощность " << std::showpoint << std::fixed << std::setprecision(0) << wd << std::endl;
    std::cout << std::noshowpoint;
    std::cout << "Введите массив источников" << std::endl;
    std::string mgLine;
    std::getline(std::cin, mgLine);
    if(mgLine.empty()) std::getline(std::cin, mgLine);
    auto vals=parseNumbers(mgLine);
    for(size_t idx=0; idx<vals.size() && idx<static_cast<size_t>(ki); ++idx) mg[idx+1]=vals[idx];
    std::cout << " Массив источников:" << std::endl;
    for(int i=1;i<=ki;i++) {
        std::cout << " Мощность " << std::setw(2) << i << "-го источника " << std::showpoint << std::fixed << std::setprecision(0) << mg[i] << std::endl;
    }
    std::cout << std::noshowpoint;

    int jn=0, neli; double lambd, nui, pbi;
    while(true) {
        std::cout << "Введите номер элемента" << std::endl;
        std::cin >> neli;
        if(neli==0) break;
        jn++;
        nel[jn]=neli;
        std::cout << "Введите: лямбда, ну, Рв" << std::endl;
        { std::string line; std::getline(std::cin, line); if(line.empty()) std::getline(std::cin, line); auto arr=parseNumbers(line); if(arr.size()>0) lambd=arr[0]; if(arr.size()>1) nui=arr[1]; if(arr.size()>2) pbi=arr[2]; }
        lambda[neli]=lambd;
        nu[neli]=nui;
        pb[neli]=pbi;
        if(jn>=m) break;
    }

    std::cout << "Диагностическая матрица" << std::endl;
    std::cout << "Вввод Диагностической матрицы:" << std::endl;
    std::cout << "---------------" << std::endl;
    int nn1, nn2;
    for(int i=1;i<=kcm;i++) {
        std::cout << "Введите номера столбцов" << std::setw(2) << i << "-ой строки" << std::endl;
        { std::string line; std::getline(std::cin,line); if(line.empty()) std::getline(std::cin,line); auto arr=parseNumbers(line); nn1=arr.size()>0?static_cast<int>(arr[0]):0; nn2=arr.size()>1?static_cast<int>(arr[1]):0; }
        for(int j=1;j<=m;j++) mdi[j]=0;
        mdi[nn1]=1; mdi[nn2]=1;
        for(int j=1;j<=m;j++) {
            md[i][j]=mdi[j]==1;
        }
    }

    std::cout << "При определении источников" << std::endl;
    std::cout << "проверять состояние их щитов?" << std::endl;
    std::cout << "Да - введите 1" << std::endl;
    std::cout << "Нет - введите 0" << std::endl;
    int otv_input; std::cin >> otv_input;
    bool otv = otv_input!=0;

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "|  №  |  лямбда |    ню   |  Рв   |" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    for(int j=1;j<=m;j++) {
        double amb=lambda[j]*1e6;
        double u1=nu[j]*1e3;
        std::cout << "| " << std::setw(3) << j << " | "
                  << std::setw(7) << std::fixed << std::setprecision(2) << amb << " | "
                  << std::setw(7) << std::fixed << std::setprecision(2) << u1 << " | "
                  << std::setw(5) << std::fixed << std::setprecision(3) << pb[j] << " |" << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "     РЕЗУЛЬТАТЫ МОДЕЛИРОВАНИЯ" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "| Д |  дельта | n исп. |  кг  |    В  |" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(0.0,1.0);

    double ts=0.0, th=0.0, tm=0.0, ti=0.0, tsl=0.0;
    double sqt=0.0, l=0.0;
    int ns_int=100;

    while(true) {
        for(int n=1;n<=100;n++) {
            for(int i=1;i<=m;i++) {
                c[i]=true; t[i]=0.0;
            }
            th=0.0; tsl=0.0; bool cp=true; bool cq=true;
            for(int i=1;i<=m;i++) {
                double z;
                do { z = dist(rng); } while(z<=0.0);
                tau[i] = -std::log(z)/lambda[i];
                t[i] = tau[i];
            }
            while(true) {
                tm=t[1]; int k=1;
                for(int i=2;i<=m;i++) {
                    if(t[i]<tm) { tm=t[i]; k=i; }
                }
                if(tm>ta) {
                    if(!cp) {
                        ti=ta-th; ts+=ti; tsl+=ti;
                    }
                    break;
                }
                double z;
                if(!c[k]) {
                    c[k]=true;
                    do { z=dist(rng); } while(z<=0.0);
                    tau[k] = -std::log(z)/lambda[k];
                } else {
                    c[k]=false;
                    z = dist(rng);
                    if(z<pb[k]) {
                        do { z=dist(rng); } while(z<=0.0);
                        tau[k] = -std::log(z)/nu[k];
                    } else {
                        tau[k] = ta;
                    }
                }
                t[k]+=tau[k];
                bool ch;
                if(c[k] && cp) {
                    ch=true;
                } else {
                    double w=0.0;
                    if(otv) {
                        for(int i=1;i<=ki;i++) {
                            if(c[i] && c[ki+i]) w+=mg[i];
                        }
                    } else {
                        for(int i=1;i<=ki;i++) {
                            if(c[i]) w+=mg[i];
                        }
                    }
                    ch=true;
                    if(w<wd) ch=false;
                }
                if(!cp && ch) {
                    ti=tm-th; ts+=ti; tsl+=ti;
                } else if(cp && !ch) {
                    th=tm; cq=false;
                }
                cp=ch;
            }
            sqt += tsl*tsl;
            if(!cq) l += 1.0;
        }
        double kg = l - ts/ta/ns_int;
        double r = 1.0 - l/ns_int;
        double d = zdelta*std::sqrt(r*(l-r)/ns_int);
    std::cout << "| " << fortranF(d,7,5)
              << " | " << fortranF(delta,5,3)
              << " | " << fortranF(static_cast<double>(ns_int),6,0)
              << " | " << fortranF(kg,7,5)
              << " | " << fortranF(r,7,5) << " |" << std::endl;
        if(l>0.0 && d<=dd) break;
        ns_int += 100;
        if(ns_int>1000) break;
    }
    std::cout << "---------------------------------------" << std::endl;
    return 0;
}

