#pragma once



#define OMEGA_FRAC 0.9   // ?_*/?_crit в [0, 1)
#define R_POLE     1.0   // полярный радиус (масштаб Rb=15)
#define Omega 0.0
#define N 100// 620 // 1240 // 620 // 402 // 350 // 7167 //1792 //1792                 // Количество ячеек по x
#define M 50 // // 256 // 448 // 512 // 256 // 256 // 256 // //1280 //1280                 // Количество ячеек по y  - должно делиться на 4 для квадруполя и на 2 для диполя
// Квадруполярная сетка работает только для 256
#define K (N*M)                // Количество ячеек в сетке
#define Rb (15.0)  // (6.0)             // Внешний радиус сетки
#define qb (1.01) // (1.005) // (1.01) // (1.01999785) // (1.02)            // Сгущение сетки каждая следующая ширина на (qb - 1)% больше предыдущей
#define dphi (pi/M)           
#define qphi (1.025) // (1.0125) // (1.025) // (1.02)              // можно настроить; при большом M брать близким к 1
#define M_HALF (M / 2)           // предполагаем, что M чётное
#define NQ (M / 4)          // M должно делиться на 4
#define BX 16   // Размеры блока потоков  
#define BY 16   // Размеры блока потоков 
#define HALO 2
#define SHARED_I (BX + 2*HALO)   // = 20
#define SHARED_J (BY + 2*HALO)   // = 20


// Следующие определения для ядра  compute_cell_ez_and_slopes
#define SLOPE_BX BX   // TODO(user): согласуйте с реальным блоком запуска
#define SLOPE_BY BY
#define SLOPE_SHARED_I (SLOPE_BX + 2)   // +1 halo слева, +1 halo справа
#define SLOPE_SHARED_J (SLOPE_BY + 2)   // +1 halo снизу, +1 halo сверху


// Предполагается, что индексы ячеек i (по радиусу) и j (по углу) отсчитываются от 0
// 
// ----------------- Радиальное разбиение -----------------

// Ширина первой радиальной ячейки (аналитически из суммы геометрической прогрессии)
#define DR1 ((Rb - 1.0) * (qb - 1.0) / (pow(qb, N) - 1.0))

// Радиус i-й границы (i = 0..N)
#define R_EDGE(i) (1.0 + DR1 * (pow(qb, (i)) - 1.0) / (qb - 1.0))

// Ширина i-й радиальной ячейки (i = 0..N-1)
#define DR(i) (DR1 * pow(qb, (i)))

// Геометрический центр
#define R_CENTER(i, j) (0.5 * (R_EDGE(i) + R_EDGE(i + 1)))


// Радиус i-й ?-границы при угле phi
#define R_EDGE_BF(i, phi) ( Roche_R(phi) + \
                            (Rb - Roche_R(phi)) * (pow(qb,(i)) - 1.0) / (pow(qb,N) - 1.0) )

// «Радиальная» ширина i-й ячейки при угле phi (для ?-граней!)
#define DR_BF(i, phi) ( (Rb - Roche_R(phi)) * (qb - 1.0) * pow(qb,(i)) / (pow(qb,N) - 1.0) )

// Приблизительный центр ячейки (r в точке ?_center)
#define R_CENTER_BF(i, j) ( 0.5 * (R_EDGE_BF(i, PHI_CENTER(j)) + R_EDGE_BF(i+1, PHI_CENTER(j))) )


// ----------------- Угловое разбиение неравномерное для диполя ----------------------------------------
// Минимальный угловой шаг (у phi = 0)
#define DPHI_MIN ((pi / 2.0) * (qphi - 1.0) / (pow(qphi, M_HALF) - 1.0))

// Ширина j-й угловой ячейки (j = 0..M-1)
#define DPHI(j) ((j) < M_HALF ? \
                 DPHI_MIN * pow(qphi, M_HALF - 1 - (j)) : \
                 DPHI_MIN * pow(qphi, (j) - M_HALF))

// Угловая граница с индексом k (k = 0..M)
#define PHI_EDGE(k) ((k) <= M_HALF ? \
                     (-(pi) / 2.0 + DPHI_MIN * (pow(qphi, M_HALF) - pow(qphi, M_HALF - (k))) / (qphi - 1.0)) : \
                     (DPHI_MIN * (pow(qphi, (k) - M_HALF) - 1.0) / (qphi - 1.0)))

// Нижняя и верхняя границы j-й ячейки
#define PHI_LEFT(j)  (PHI_EDGE(j))
#define PHI_RIGHT(j) (PHI_EDGE((j) + 1))

// Центр j-й ячейки по углу (средний угол)
#define PHI_CENTER(j) (0.5 * (PHI_LEFT(j) + PHI_RIGHT(j)))

// Внутренняя дуга (при r = R_EDGE(i))
#define L_INNER(i,j) (R_EDGE(i) * DPHI(j))

// Внешняя дуга (при r = R_EDGE(i+1))
#define L_OUTER(i,j) (R_EDGE(i + 1) * DPHI(j))


// Левый и правый радиальные отрезки (равны DR(i))
#define L_LEFT(i) (DR(i))
#define L_RIGHT(i) (DR(i))

// ----------------- Площадь (объём) ячейки -----------------

//#define CELL_AREA(i) (0.5 * (R_EDGE(i + 1) * R_EDGE(i + 1) - R_EDGE(i) * R_EDGE(i)) * dphi)   // равномерный угол
#define CELL_AREA(i,j) (0.5 * (R_EDGE(i + 1) * R_EDGE(i + 1) - R_EDGE(i) * R_EDGE(i)) * DPHI(j))   // неравномерный угол



#define krit 0.25 // 0.45 // 0.1 // 0.25
#define ggg (5.0/3.0)          // Показатель адиабаты
#define ga (5.0/3.0)          // Показатель адиабаты
#define g1 (ga - 1.0)
#define gg1 (ga - 1.0)
#define g2 (ga + 1.0)
#define gg2 (ga + 1.0)
#define gp ((g2/ga)/2.0)
#define gm ((g1/ga)/2.0)
#define gga ga
#define kv(x) ((x)*(x))
#define kv3(x) ((x)*(x)*(x))
#define skk(u,v,w,bx,by,bz) ( (u)*(bx) + (v)*(by) + (w)*(bz) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))
#define U8(ro, p, u, v, w, bx, by, bz)  ( (p) / (ggg - 1.0) + 0.5 * (ro) * kvv(u,v,w) + kvv(bx,by,bz) / 25.13274122871834590768)

#define pi 3.14159265358979323846
#define PI 3.14159265358979323846
#define cpi4 12.56637061435917295384
#define cpi8 25.13274122871834590768
#define spi4 ( 3.544907701811032 )

#define M_inf  0.0 // 0.7 // 0.4 // 0.8
#define M_infty  M_inf
#define phi_0  4.0 // 2.0 // 17.0 // 4.878 //1.627
#define alpha  45.0
#define M_alf  8.0 // 12.0
#define epsilon_ (1.0/M_alf)
#define step  70000
#define omega 0.0 //6 //1600
#define M_0  10.0

#define ddist 0.3 // 0.45 // 0.65 // 0.85

#define kk_ 196.0
//#define chi 17.0
#define rr_0  0.07    // Внутренняя сфера (на которой ставятся граничные условия - реально она дальше)


extern __device__ int sign_(const double& x);
extern __device__ double minmod_(double x, double y);
extern __device__ double linear_(double x1, double t1, double x2, double t2, double x3, double t3, double y);
extern __device__ void linear2_(double x1, double t1, double x2, double t2, double x3, double t3, double y1, double y2,//
    double& A, double& B);
extern __device__ void TVD(const double2& s_1, const double2& s_2, const double2& s_3, const double2& s_4, const double2& s_5,//
    const double2& s_6, const double2& s_7, const double2& s_8, const double2& s_9, double2& s12,//
    double2& s13, double2& s14, double2& s15, double2& s21, double2& s31, double2& s41, double2& s51, double dx, double dy, bool zero);
extern __device__ double HLLC_Korolkov_2D(const double2& Ls, const double2& Lu, const double2& Rs, const double2& Ru,//
    const double n1, const double n2, double2& Ps, double2& Pu, const double rad);
extern __device__ double HLLCQ_Korolkov_2D(const double2& Ls, const double2& Lu, const double2& Rs, const double2& Ru,//
    const double& LQ, const double& RQ, double n1, double n2, double2& Ps, double2& Pu, double& PQ, double rad);
extern __device__ double HLLCQ_Aleksashov(const double2& Ls, const double2& Lu, const double2& Rs, const double2& Ru,//
    const double& LQ, const double& RQ, double n1, double n2, double2& Ps, double2& Pu, double& PQ, double rad);
extern __device__ double HLLC_Aleksashov_2D(double2& Ls, double2& Lu, double2& Rs, double2& Ru,//
    double n1, double n2, double2& Ps, double2& Pu, double rad);
__device__ double HLLDQ_Korolkov(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
    const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
    const double& Bx_R, const double& By_R, const double& Bz_R, double* P, double& PQ, const double& n1, const double& n2, const double& n3, double& SL_, double& SR_, const double& rad, int metod, double x = 0.0, double y = 0.0);