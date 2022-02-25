#include "fferrpar.h"

//Bs->Kstr,ll
BdtoKstrll_fferrpar::BdtoKstrll_fferrpar(){
    tpd=pow(mBd()+mKst(),2.0);
    tmd=pow(mBd()-mKst(),2.0);
    tzd=tpd*(1.0-sqrt(1.0-tmd/tpd));

    m_PSbs = 5.336, m_Vbs = 5.412, m_Abs = 5.829;
}

////######refer to "CholeskyDecomp.nb"
//central vector:\mu
double BdtoKstrll_fferrpar::cen_FF[] = {0.369196, -1.36584, 0.128191, 0.29725, 0.392378, 1.18916, /*0.265375,*/ 0.533638, 0.483166, 0.376313, -1.16597, 2.42443,
                    0.312055, -1.00893, 1.5272, /*0.312055,*/ 0.496846, 1.61431, 0.667412, 1.31812, 3.82334};
//uncertainty vector:\sigma
double BdtoKstrll_fferrpar::unc_FF[] = {0.0289419, 0.256849, 1.63438, 0.026356, 0.187894, 1.02531, /*0.0208033,*/ 0.128772, 0.656273, 0.0332944, 0.261268, 1.53102,
                    0.027496, 0.189575, 1.63965, /*0.027496,*/ 0.166219, 0.803783, 0.0633333, 0.222129, 2.20338};
//Cholesky Decomposition of the covariance matrix
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

double BdtoKstrll_fferrpar::A0a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,0);}
double BdtoKstrll_fferrpar::A0a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,1);}
double BdtoKstrll_fferrpar::A0a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,2);}
double BdtoKstrll_fferrpar::A1a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,3);}
double BdtoKstrll_fferrpar::A1a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,4);}
double BdtoKstrll_fferrpar::A1a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,5);}
double BdtoKstrll_fferrpar::A12a0(double unv[]){return A0a0(unv)*(pow(mBd(),2.0)-pow(mKst(),2.0))/(8.0*mBd()*mKst());}//Eq. 17 of arXiv:1503.05534
double BdtoKstrll_fferrpar::A12a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,6);}
double BdtoKstrll_fferrpar::A12a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,7);}
double BdtoKstrll_fferrpar::Va0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,8);}
double BdtoKstrll_fferrpar::Va1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,9);}
double BdtoKstrll_fferrpar::Va2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,10);}
double BdtoKstrll_fferrpar::T1a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,11);}
double BdtoKstrll_fferrpar::T1a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,12);}
double BdtoKstrll_fferrpar::T1a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,13);}
double BdtoKstrll_fferrpar::T2a0(double unv[]){return T1a0(unv);}//Eq. 17 of arXiv:1503.05534
double BdtoKstrll_fferrpar::T2a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,14);}
double BdtoKstrll_fferrpar::T2a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,15);}
double BdtoKstrll_fferrpar::T23a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,16);}
double BdtoKstrll_fferrpar::T23a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,17);}
double BdtoKstrll_fferrpar::T23a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,18);}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
Bstophill_fferrpar::Bstophill_fferrpar(){
    tps=pow(mBs()+mphi(),2.0);
    tms=pow(mBs()-mphi(),2.0);
    tzs=tps*(1.0-sqrt(1.0-tms/tps));

    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;
}

////######refer to "CholeskyDecomp.nb"
//central vector:\mu
double Bstophill_fferrpar::cen_FF[] = {0.421328, -0.976454, 3.2714, 0.288007, 0.350826, 1.69688, /*0.267053,*/ 0.954402, 2.15263, 0.364478, -1.22389, 3.74061,
                    0.299475, -1.1013, 0.58459, /*0.299475,*/ 0.403564, 1.03987, 0.65233, 2.09622, 6.73572};
//uncertainty vector:\sigma
double Bstophill_fferrpar::unc_FF[] = {0.0240514, 0.23757, 1.35909, 0.0105759, 0.103813, 0.790297, /*0.0152447,*/ 0.125879, 0.47881, 0.0141353, 0.164176, 1.72682,
                    0.0120688, 0.0835241, 1.00324, /*0.0120688,*/ 0.0803563, 0.608872, 0.0357469, 0.330082, 1.79526};
//Cholesky Decomposition of the covariance matrix
double Bstophill_fferrpar::chd_cov[][MXdm] =
    {{0.0240514, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.163293, 0.172554, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {-0.475215, 0.676253, 1.07888, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0.},
    {0.00218642, -0.00197483, 0.00342496, 0.00956237, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {-0.0576101, -0.0326219, 0.000174754, 0.0351823, 0.0718067, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {-0.472299, -0.28095, 0.23947, -0.127972, 0.451391, 0.212353, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.0972527, 0.0451694, -0.0316008, -0.0205043, -0.0118269, -0.0110278, 0.0516394, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.203317, 0.233571, 0.0438621, -0.0914031, -0.093997, 0.0486694, 0.329933, 0.0550212, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.00265482, -0.00217211, 0.00622126, 0.00306886, -0.000528257, 0.00110301, -0.000239527, 0.00151878, 0.0116644, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {-0.0693557, -0.0524, -0.00273018, 0.0388572, 0.0539049, -0.0223373, -0.00321119, 0.0329628, 0.0467491, 0.105799, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {-0.788385, -0.828968, 0.00831865, 0.215627, 0.263561, -0.0923733, 0.00957827, 0.109689, -0.287164, 0.91525, 0.785072, 0., 0., 0., 0., 0., 0., 0., 0.},
    {0.00104471, -0.0000830111, 0.00328639, 0.0021387, 0.00222748, 0.00208524, -0.00243335, 0.000919938, 0.00486207, 0.000834824, -0.00208809, 0.00918745,
        0., 0., 0., 0., 0., 0., 0.},
    {-0.0248073, -0.00493503, -0.00851714, 0.0171202, 0.0384485, -0.00269388, -0.0218469, 0.0127849, 0.000661001, 0.0237323, -0.0174053, 0.0288518, 0.0463202,
        0., 0., 0., 0., 0., 0.},
    {0.0673318, -0.0815351, -0.19776, -0.143991, -0.232604, 0.0413882, -0.110348, -0.119071, -0.39298, -0.306509, 0.177287, -0.338206, 0.353364, 0.578129,
        0., 0., 0., 0., 0.},
    {-0.0367423, -0.0234062, -0.0107079, 0.0243058, 0.0402967, -0.00759992, -0.0167143, 0.00452568, 0.00111411, 0.00961089, -0.00369442, 0.00686867,
        0.00901169, 0.00420305, 0.0402382, 0., 0., 0., 0.},
    {-0.119938, -0.163033, -0.150601, 0.00982777, -0.00233557, -0.0186944, -0.0751897, -0.0586033, -0.18926, -0.138564, 0.0923495, -0.224566, 0.107444, 0.240453,
        0.258727, 0.217396, 0., 0., 0.},
    {0.0251552, 0.00769765, -0.00525245, 0.0014018, -0.00808339, 0.00323918, -0.00310673, -0.00219015, 0.00462386, -0.0112424, -0.00198258, 0.00251315, 0.00149322,
        0.00422493, 0.00154043, 0.00445653, 0.0162984, 0., 0.},
    {0.240301, 0.125205, -0.0548407, -0.0358827, -0.0542308, 0.0033148, 0.0587943, 0.00334095, 0.0236515, -0.0661997, 0.0141647, -0.0504951, -0.00243042,
        0.0336393, -0.00886644, -0.00027137, 0.0521349, 0.114857, 0.},
    {1.08095, 0.572936, -0.0735966, -0.411486, -0.439335, 0.124, 0.389559, 0.0129582, 0.0233923, -0.364223, 0.232819, -0.237481, -0.0124845, 0.475836, -0.098462,
        0.2244, -0.203678, 0.776448, 0.12808}};

double Bstophill_fferrpar::zs(double qsq){
    return (sqrt(tps-qsq)-sqrt(tps-tzs))/(sqrt(tps-qsq) +sqrt(tps-tzs));}

double Bstophill_fferrpar::mPSbs(){return m_PSbs;}
double Bstophill_fferrpar::mVbs(){return m_Vbs;}
double Bstophill_fferrpar::mAbs(){return m_Abs;}

double Bstophill_fferrpar::A0a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,0);}
double Bstophill_fferrpar::A0a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,1);}
double Bstophill_fferrpar::A0a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,2);}
double Bstophill_fferrpar::A1a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,3);}
double Bstophill_fferrpar::A1a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,4);}
double Bstophill_fferrpar::A1a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,5);}
double Bstophill_fferrpar::A12a0(double unv[]){return A0a0(unv)*(pow(mBd(),2.0)-pow(mKst(),2.0))/(8.0*mBd()*mKst());}//Eq. 17 of arXiv:1503.05534
double Bstophill_fferrpar::A12a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,6);}
double Bstophill_fferrpar::A12a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,7);}
double Bstophill_fferrpar::Va0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,8);}
double Bstophill_fferrpar::Va1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,9);}
double Bstophill_fferrpar::Va2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,10);}
double Bstophill_fferrpar::T1a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,11);}
double Bstophill_fferrpar::T1a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,12);}
double Bstophill_fferrpar::T1a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,13);}
double Bstophill_fferrpar::T2a0(double unv[]){return T1a0(unv);}//Eq. 17 of arXiv:1503.05534
double Bstophill_fferrpar::T2a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,14);}
double Bstophill_fferrpar::T2a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,15);}
double Bstophill_fferrpar::T23a0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,16);}
double Bstophill_fferrpar::T23a1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,17);}
double Bstophill_fferrpar::T23a2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,18);}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
BdtoKll_fferrpar::BdtoKll_fferrpar(){
    tp=pow(mBd()+mK(),2.0);
    tm=pow(mBd()-mK(),2.0);
    tz=tp*(1.0-sqrt(1.0-tm/tp));

    m_Vbs=5.412,m_Sbs=5.630;
}

////######refer to "CholeskyDecomp.nb"
//central vector:\mu
double BdtoKll_fferrpar::cen_FF[] = {0.195117, -0.446126, 0.299383, -0.773546, 0.00955438, 0.32909, -0.866947, 0.00609567};
//uncertainty vector:\sigma
double BdtoKll_fferrpar::unc_FF[] = {0.168097, 0.408946, 0.025872, 0.150128, 0.871576, 0.0277282, 0.138144, 0.750826};
//Cholesky Decomposition of the covariance matrix
double BdtoKll_fferrpar::chd_cov[][MXdm] =
    {{0.168097, 0., 0., 0., 0., 0., 0., 0.},
    {0.344397, 0.220517, 0., 0., 0., 0., 0., 0.},
    {0.0148271, -0.0122488, 0.0173056, 0., 0., 0., 0., 0.},
    {-0.00791751, 0.00696742, -0.0572808, 0.13837, 0., 0., 0., 0.},
    {-0.215666, 0.19584, -0.438077, 0.361481, 0.593464, 0., 0., 0.},
    {0.016741, -0.0171203, 0.0082008, 0.00286598, 0.00132161, 0.0108753, 0., 0.},
    {0.0102572, 0.0286941, -0.00663218, 0.0147008, 0.00223952, -0.0279217, 0.130807, 0.},
    {-0.170731, 0.411061, -0.0572078, 0.0390342, 0.0421714, -0.156681, 0.490221, 0.306888}};

double BdtoKll_fferrpar::z(double qsq){return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));}
double BdtoKll_fferrpar::mVbs(){return m_Vbs;}
double BdtoKll_fferrpar::mSbs(){return m_Sbs;}

double BdtoKll_fferrpar::fza1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,0);}
double BdtoKll_fferrpar::fza2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,1);}
double BdtoKll_fferrpar::fTa0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,2);}
double BdtoKll_fferrpar::fTa1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,3);}
double BdtoKll_fferrpar::fTa2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,4);}
double BdtoKll_fferrpar::fpa0(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,5);}
double BdtoKll_fferrpar::fpa1(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,6);}
double BdtoKll_fferrpar::fpa2(double unv[]){return mnd_cov(cen_FF,chd_cov,unv,7);}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////













