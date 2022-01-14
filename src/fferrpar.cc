

#include "smpar.h"
#include "myfun.h"
#include "fferrpar.h"

//Bs->Kstr,ll
BdtoKstrll_fferrpar::BdtoKstrll_fferrpar(){
    tpd=pow(mBd()+mKst(),2.0);
    tmd=pow(mBd()-mKst(),2.0);
    tzd=tpd*(1.0-sqrt(1.0-tmd/tpd));

    m_PSbs = 5.366, m_Vbs = 5.415, m_Abs = 5.829;

    A0_a0 = 0.369196, A0_a1 = -1.36584, A0_a2 = 0.128191,
    A1_a0 = 0.29725, A1_a1 = 0.392378, A1_a2 = 1.18916,
    A12_a0 = 0.265375, A12_a1 = 0.533638, A12_a2 = 0.483166,
    V_a0 = 0.376313, V_a1 = -1.16597, V_a2 = 2.42443,
    T1_a0 = 0.312055, T1_a1 = -1.00893, T1_a2 = 1.5272,
    T2_a0 = 0.312055, T2_a1 = 0.496846, T2_a2 = 1.61431,
    T23_a0 = 0.667412, T23_a1 = 1.31812, T23_a2 = 3.82334;

    eA0_a0 = 0.0289419, eA0_a1 = 0.256849, eA0_a2 = 1.63438,
    eA1_a0 = 0.026356, eA1_a1 = 0.187894 , eA1_a2 = 1.02531,
    eA12_a0 = 0.0208033, eA12_a1 = 0.128772, eA12_a2 = 0.656273,
    eV_a0 = 0.0332944, eV_a1 = 0.261268, eV_a2 = 1.53102,
    eT1_a0 = 0.027496, eT1_a1 = 0.189575, eT1_a2 = 1.63965,
    eT2_a0 = 0.027496, eT2_a1 = 0.166219, eT2_a2 = 0.803783,
    eT23_a0 = 0.0633333, eT23_a1 = 0.222129, eT23_a2 = 2.20338;}

double BdtoKstrll_fferrpar::cov_A0[][mxdm] = {{0.000837635, 0.00471065, 0.00272131}, {0.00471065, 0.0659712, 0.204148}, {0.00272131, 0.204148, 2.67119}};
double BdtoKstrll_fferrpar::cov_A1[][mxdm] = {{0.000694641, 0.00347885, 0.00415084}, {0.00347885, 0.0353042, 0.144041}, {0.00415084, 0.144041, 1.05126}};
double BdtoKstrll_fferrpar::cov_A12[][mxdm] = {{0.000432775, 0.00165482, 0.000605487}, {0.00165482, 0.0165823, 0.059165}, {0.000605487, 0.059165, 0.430695}};
double BdtoKstrll_fferrpar::cov_V[][mxdm] = {{0.00110852, 0.00658825, -0.0202371}, {0.00658825, 0.0682608, 0.0138459}, {-0.0202371, 0.0138459, 2.34402}};
double BdtoKstrll_fferrpar::cov_T1[][mxdm] = {{0.000756032, 0.00372113, -0.0213639}, {0.00372113, 0.0359388, 0.00028928}, {-0.0213639, 0.00028928, 2.68844}};
double BdtoKstrll_fferrpar::cov_T2[][mxdm] = {{0.000756032, 0.00316639, -0.00503577}, {0.00316639, 0.0276287, 0.0511444}, {-0.00503577, 0.0511444, 0.646067}};
double BdtoKstrll_fferrpar::cov_T23[][mxdm] = {{0.0040111, 0.00460369, -0.0993374}, {0.00460369, 0.0493411, 0.140689}, {-0.0993374, 0.140689, 4.85487}};

////######refer to "MCError.nb"
double BdtoKstrll_fferrpar::cen_FF[] = {0.369196, -1.36584, 0.128191, 0.29725, 0.392378, 1.18916, /*0.265375,*/ 0.533638, 0.483166, 0.376313, -1.16597, 2.42443,
                                        0.312055, -1.00893, 1.5272, /*0.312055,*/ 0.496846, 1.61431, 0.667412, 1.31812, 3.82334};

double BdtoKstrll_fferrpar::unc_FF[] = {0.0289419, 0.256849, 1.63438, 0.026356, 0.187894, 1.02531, /*0.0208033,*/ 0.128772, 0.656273, 0.0332944, 0.261268, 1.53102,
                                        0.027496, 0.189575, 1.63965, /*0.027496,*/ 0.166219, 0.803783, 0.0633333, 0.222129, 2.20338};

double BdtoKstrll_fferrpar::chd_cov[][MXdm] =
    {{0.0289419, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.162762, 0.198695, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0940266, 0.950422, 1.32629, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.00162188, 0.00891175, 0.0102406, 0.0225327, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0317625, 0.0886649, 0.10637, 0.0686952, 0.101982, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.246737, 0.224089, 0.657956, -0.221199, 0.603478, 0.306829, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0784861, 0.0616464, 0.00207627, -0.0293886, 0.00707759, -0.0116239, 0.074624, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0232996, 0.378901, 0.264747, -0.114241, 0.0456552, 0.13923, 0.420329, 0.0727913, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.00130952, 0.0109986, 0.0144095, 0.0249514, -0.000190672, 0.0000253352, -0.0018253, -0.00129291, 0.0122715, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0465494, 0.148918, 0.125378, 0.104186, 0.070245, -0.0265486, 0.0128345, -0.00552823, 0.0418479, 0.0987799, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.299193, -0.0961989, 0.0639645, -0.654258, 0.185784, 0.151818, 0.131343, -0.0957442, -0.32762, 0.778137, 1.00811, 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.000453363, 0.00820388, 0.010197, 0.0210307, -0.000111574, -0.000397777, -0.000461424, -0.000455808, 0.00500313, 0.00177127, -0.000787501, 0.0106242,
     0., 0., 0., 0., 0., 0., 0.},
    {0.00858101, 0.102937, 0.0831889, 0.0718722, 0.0381135, -0.0298683, 0.000354509, -0.00864035, 0.0312983, 0.0648098, -0.0171633, 0.0203925, 0.0698091,
     0., 0., 0., 0., 0., 0.},
    {-0.259814, -0.174126, -0.110607, -0.760751, -0.266877, -0.149465, -0.0498694, -0.21357, 0.0440043, 0.42681, 0.159031, -0.353074, 0.989719, 0.73782,
     0., 0., 0., 0., 0.},
    {0.0249766, 0.0859914, 0.0802007, 0.0624074, 0.0720445, -0.019893, -0.00292433, -0.00669054, 0.0161842, 0.0303456, -0.0061512, 0.0137381, 0.0184727,
     0.00155864, 0.043486, 0., 0., 0., 0.},
    {0.0250422, 0.167756, 0.149335, -0.326021, 0.245417, -0.0999542, 0.00833677, -0.110622, 0.00710868, 0.233226, 0.0379235, -0.151696, 0.347655, 0.208349,
     0.286429, 0.284555, 0., 0., 0.},
    {0.0285681, 0.0235534, 0.0328473, 0.0155197, 0.00710868, 0.00376057, -0.00214834, 0.00401288, 0.00451454, 0.00521994, -0.000461975, 0.0000946369,
     0.00122964, -0.00202422, -0.00182327, -0.00416626, 0.0340714, 0., 0.},
    {0.0233036, 0.134117, 0.0190439, 0.0161651, 0.0325043, -0.02255, 0.0698608, 0.0136259, 0.00191628, -0.0102979, 0.0133567, -0.0141345, 0.00767429,
     0.00890779, 0.0129547, 0.000645929, -0.00178231, 0.150857, 0.},
    {-0.765568, 0.0744983, -0.862889, -0.51616, -0.366398, -0.22148, 0.421026, -0.142154, -0.033986, -0.0336954, 0.192015, 0.0914068,
     0.439198, 0.520461, 0.213893, 0.657043, -0.995098, 0.916255, 0.230145}};


double BdtoKstrll_fferrpar::zd(double qsq){
    return (sqrt(tpd-qsq)-sqrt(tpd-tzd))/(sqrt(tpd-qsq) +sqrt(tpd-tzd));}

double BdtoKstrll_fferrpar::mPSbs(){return m_PSbs;}
double BdtoKstrll_fferrpar::mVbs(){return m_Vbs;}
double BdtoKstrll_fferrpar::mAbs(){return m_Abs;}

double BdtoKstrll_fferrpar::A0a0(){return mnd_cov(cen_FF,chd_cov,0);}
double BdtoKstrll_fferrpar::A0a1(){return mnd_cov(cen_FF,chd_cov,1);}
double BdtoKstrll_fferrpar::A0a2(){return mnd_cov(cen_FF,chd_cov,2);}
double BdtoKstrll_fferrpar::A1a0(){return mnd_cov(cen_FF,chd_cov,3);}
double BdtoKstrll_fferrpar::A1a1(){return mnd_cov(cen_FF,chd_cov,4);}
double BdtoKstrll_fferrpar::A1a2(){return mnd_cov(cen_FF,chd_cov,5);}
double BdtoKstrll_fferrpar::A12a0(){return A12_a0 + ((A0a0()-A0_a0)/eA0_a0)*eA12_a0;}
double BdtoKstrll_fferrpar::A12a1(){return mnd_cov(cen_FF,chd_cov,6);}
double BdtoKstrll_fferrpar::A12a2(){return mnd_cov(cen_FF,chd_cov,7);}
double BdtoKstrll_fferrpar::Va0(){return mnd_cov(cen_FF,chd_cov,8);}
double BdtoKstrll_fferrpar::Va1(){return mnd_cov(cen_FF,chd_cov,9);}
double BdtoKstrll_fferrpar::Va2(){return mnd_cov(cen_FF,chd_cov,10);}
double BdtoKstrll_fferrpar::T1a0(){return mnd_cov(cen_FF,chd_cov,11);}
double BdtoKstrll_fferrpar::T1a1(){return mnd_cov(cen_FF,chd_cov,12);}
double BdtoKstrll_fferrpar::T1a2(){return mnd_cov(cen_FF,chd_cov,13);}
double BdtoKstrll_fferrpar::T2a0(){return T1a0();}
double BdtoKstrll_fferrpar::T2a1(){return mnd_cov(cen_FF,chd_cov,14);}
double BdtoKstrll_fferrpar::T2a2(){return mnd_cov(cen_FF,chd_cov,15);}
double BdtoKstrll_fferrpar::T23a0(){return mnd_cov(cen_FF,chd_cov,16);}
double BdtoKstrll_fferrpar::T23a1(){return mnd_cov(cen_FF,chd_cov,17);}
double BdtoKstrll_fferrpar::T23a2(){return mnd_cov(cen_FF,chd_cov,18);}


//Bs->phi,ll
Bstophill_fferrpar::Bstophill_fferrpar(){
    tps=pow(mBs()+mKst(),2.0);
    tms=pow(mBs()-mKst(),2.0);
    tzs=tps*(1.0-sqrt(1.0-tms/tps));

    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;

    V_a0=0.376313,V_a1=-1.16597,V_a2=2.42443,
    A0_a0=0.369196,A0_a1=-1.36584,A0_a2=0.128191,
    A1_a0=0.29725,A1_a1=0.392378,A1_a2=1.18916,
    A12_a0=0.265375,A12_a1=0.533638,A12_a2=0.483166,
    T1_a0=0.312055,T1_a1=-1.00893,T1_a2=1.5272,
    T2_a0=0.312055,T2_a1=0.496846,T2_a2=1.61431,
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;

    eV_a0=0.0376313,eV_a1=0.116597,eV_a2=0.242443,
    eA0_a0=0.0369196,eA0_a1=0.136584,eA0_a2=0.0128191,
    eA1_a0=0.029725,eA1_a1=0.0392378,eA1_a2=0.118916,
    eA12_a0=0.0265375,eA12_a1=0.0533638,eA12_a2=0.0483166,
    eT1_a0=0.0312055,eT1_a1=0.100893,eT1_a2=0.15272,
    eT2_a0=0.0312055,eT2_a1=0.0496846,eT2_a2=0.161431,
    eT23_a0=0.0667412,eT23_a1=0.131812,eT23_a2=0.382334;}

double Bstophill_fferrpar::zs(double qsq){
    return (sqrt(tps-qsq)-sqrt(tps-tzs))/(sqrt(tps-qsq) +sqrt(tps-tzs));}

double Bstophill_fferrpar::mPSbs(){return m_PSbs;}
double Bstophill_fferrpar::mVbs(){return m_Vbs;}
double Bstophill_fferrpar::mAbs(){return m_Abs;}

double Bstophill_fferrpar::Va0(){return mnd(V_a0,eV_a0);}
double Bstophill_fferrpar::Va1(){return mnd(V_a1,eV_a1);}
double Bstophill_fferrpar::Va2(){return mnd(V_a2,eV_a2);}
double Bstophill_fferrpar::A0a0(){return mnd(A0_a0,eA0_a0);}
double Bstophill_fferrpar::A0a1(){return mnd(A0_a1,eA0_a1);}
double Bstophill_fferrpar::A0a2(){return mnd(A0_a2,eA0_a2);}
double Bstophill_fferrpar::A1a0(){return mnd(A1_a0,eA1_a0);}
double Bstophill_fferrpar::A1a1(){return mnd(A1_a1,eA1_a1);}
double Bstophill_fferrpar::A1a2(){return mnd(A1_a2,eA1_a2);}
double Bstophill_fferrpar::A12a0(){return mnd(A12_a0,eA12_a0);}
double Bstophill_fferrpar::A12a1(){return mnd(A12_a1,eA12_a1);}
double Bstophill_fferrpar::A12a2(){return mnd(A12_a2,eA12_a2);}
double Bstophill_fferrpar::T1a0(){return mnd(T1_a0,eT1_a0);}
double Bstophill_fferrpar::T1a1(){return mnd(T1_a1,eT1_a1);}
double Bstophill_fferrpar::T1a2(){return mnd(T1_a2,eT1_a2);}
double Bstophill_fferrpar::T2a0(){return mnd(T2_a0,eT2_a0);}
double Bstophill_fferrpar::T2a1(){return mnd(T2_a1,eT2_a1);}
double Bstophill_fferrpar::T2a2(){return mnd(T2_a2,eT2_a2);}
double Bstophill_fferrpar::T23a0(){return mnd(T23_a0,eT23_a0);}
double Bstophill_fferrpar::T23a1(){return mnd(T23_a1,eT23_a1);}
double Bstophill_fferrpar::T23a2(){return mnd(T23_a2,eT23_a2);}


//Bd->K,ll
BdtoKll_fferrpar::BdtoKll_fferrpar(){
    tp=pow(mBd()+mK(),2.0);
    tm=pow(mBd()-mK(),2.0);
    tz=tp*(1.0-sqrt(1.0-tm/tp));

    lcsr_m_Vbs=5.412,lcsr_m_Sbs=5.630;

    c_fp0=0.32909,c_fp1=-0.866947,c_fp2=0.00609567,c_fz0=0.0,c_fz1=0.195117,
    c_fz2=-0.446126,c_ft0=0.299383,c_ft1=-0.773546,c_ft2=0.00955438;

    ec_fp0=0.032909,ec_fp1=0.0866947,ec_fp2=0.000609567,ec_fz0=0.00,ec_fz1=0.0195117,
    ec_fz2=0.0446126,ec_ft0=0.0299383,ec_ft1=0.0773546,ec_ft2=0.000955438;

    lat_m_Vbs=5.4154,lat_m_Sbs=5.711;

    b0_p=0.466,b1_p=-0.885,b2_p=-0.213,b0_z=0.292,b1_z=0.281,b2_z=0.150,
    b0_t=0.460,b1_t=-1.089,b2_t=-1.114;

    eb0_p=0.0466,eb1_p=0.0885,eb2_p=0.0213,eb0_z=0.0292,eb1_z=0.0281,eb2_z=0.0150,
    eb0_t=0.0460,eb1_t=0.1089,eb2_t=0.1114;}

double BdtoKll_fferrpar::z(double qsq){return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));}
double BdtoKll_fferrpar::lcsr_mVbs(){return lcsr_m_Vbs;}
double BdtoKll_fferrpar::lcsr_mSbs(){return lcsr_m_Sbs;}

double BdtoKll_fferrpar::cfp0(){return mnd(c_fp0,ec_fp0);}
double BdtoKll_fferrpar::cfp1(){return mnd(c_fp1,ec_fp1);}
double BdtoKll_fferrpar::cfp2(){return mnd(c_fp2,ec_fp2);}
double BdtoKll_fferrpar::cfz0(){return mnd(c_fz0,ec_fz0);}
double BdtoKll_fferrpar::cfz1(){return mnd(c_fz1,ec_fz1);}
double BdtoKll_fferrpar::cfz2(){return mnd(c_fz2,ec_fz2);}
double BdtoKll_fferrpar::cft0(){return mnd(c_ft0,ec_ft0);}
double BdtoKll_fferrpar::cft1(){return mnd(c_ft1,ec_ft1);}
double BdtoKll_fferrpar::cft2(){return mnd(c_ft2,ec_ft2);}

double BdtoKll_fferrpar::lat_mVbs(){return lat_m_Vbs;}
double BdtoKll_fferrpar::lat_mSbs(){return lat_m_Sbs;}

double BdtoKll_fferrpar::b0p(){return mnd(b0_p,eb0_p);}
double BdtoKll_fferrpar::b1p(){return mnd(b1_p,eb1_p);}
double BdtoKll_fferrpar::b2p(){return mnd(b2_p,eb2_p);}
double BdtoKll_fferrpar::b0z(){return mnd(b0_z,eb0_z);}
double BdtoKll_fferrpar::b1z(){return mnd(b1_z,eb1_z);}
double BdtoKll_fferrpar::b2z(){return mnd(b2_z,eb2_z);}
double BdtoKll_fferrpar::b0t(){return mnd(b0_t,eb0_t);}
double BdtoKll_fferrpar::b1t(){return mnd(b1_t,eb1_t);}
double BdtoKll_fferrpar::b2t(){return mnd(b2_t,eb2_t);}












