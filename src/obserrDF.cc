#include "myfun.h"
#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"
#include "obserrDF.h"




//Bd->K,ll
double BdtoKll_obserrDF::Er_diffWidth(double qsq, double ml){
    return sqrt (4*pow (Erfz (qsq), 2)*
    pow (2*ml*(pow (mBd (), 2) - pow (mK (), 2) +
         qsq)*((C10Im () +
            C10RHIm ())*(((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                 pow (mK (), 2)))/(2.*(mb () - ms ())*
               fp (qsq)) + (C10Im ()*(pow (mBd (), 2) -
                 pow (mK (), 2))*ml)/(qsq*fp (qsq))) + (C10Re () +
            C10RHRe ())*(((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                 pow (mK (), 2)))/(2.*(mb () - ms ())*
               fp (qsq)) + (C10Re ()*(pow (mBd (), 2) -
                 pow (mK (), 2))*ml)/(qsq*fp (qsq)))) +
      qsq*(pow (betal (qsq, ml),
           2)*((pow (CSIm () + CSRHIm (), 2)*
               pow (pow (mBd (), 2) - pow (mK (), 2), 2)*
               fz (qsq))/(2.*pow (mb () - ms (), 2)*
               pow (fp (qsq), 2)) + (pow (CSRe () + CSRHRe (), 2)*
               pow (pow (mBd (), 2) - pow (mK (), 2), 2)*
               fz (qsq))/(2.*pow (mb () - ms (), 2)*
               pow (fp (qsq), 2))) +
         2*(((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                 pow (mK (), 2)))/(2.*(mb () - ms ())*
               fp (qsq)) + (C10Im ()*(pow (mBd (), 2) -
                 pow (mK (), 2))*ml)/(qsq*
               fp (qsq)))*(((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                 pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
               fp (qsq)) +
            C10Im ()*
             ml*(-1 + ((pow (mBd (), 2) - pow (mK (), 2))*(-1 +
                    fz (qsq)/fp (qsq)))/qsq)) +
         2*(((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                 pow (mK (), 2)))/(2.*(mb () - ms ())*
               fp (qsq)) + (C10Re ()*(pow (mBd (), 2) -
                 pow (mK (), 2))*ml)/(qsq*
               fp (qsq)))*(((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                 pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
               fp (qsq)) +
            C10Re ()*
             ml*(-1 + ((pow (mBd (), 2) - pow (mK (), 2))*(-1 +
                    fz (qsq)/fp (qsq)))/qsq))), 2)*
    pow (nf (qsq, ml), 2)*pow (xiP (qsq), 4) +
   4*pow (Erfp (qsq), 2)*
    pow (nf (qsq, ml)*
       pow (xiP (qsq),
        2)*(2*ml*(pow (mBd (), 2) - pow (mK (), 2) +
            qsq)*((C10Im () +
               C10RHIm ())*(-((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                     pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  pow (fp (qsq),
                   2)) - (C10Im ()*(pow (mBd (), 2) - pow (mK (), 2))*
                  ml*fz (qsq))/(qsq*pow (fp (qsq), 2))) + (C10Re () +
               C10RHRe ())*(-((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                     pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  pow (fp (qsq),
                   2)) - (C10Re ()*(pow (mBd (), 2) - pow (mK (), 2))*
                  ml*fz (qsq))/(qsq*pow (fp (qsq), 2)))) +
         qsq*(pow (betal (qsq, ml),
              2)*(-(pow (CSIm () + CSRHIm (), 2)*
                   pow (pow (mBd (), 2) - pow (mK (), 2), 2)*
                   pow (fz (qsq), 2))/(2.*pow (mb () - ms (), 2)*
                  pow (fp (qsq), 3)) - (pow (CSRe () + CSRHRe (), 2)*
                  pow (pow (mBd (), 2) - pow (mK (), 2), 2)*
                  pow (fz (qsq), 2))/(2.*pow (mb () - ms (), 2)*
                  pow (fp (qsq), 3))) +
            2*(-((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                    pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  pow (fp (qsq),
                   2)) - (C10Im ()*(pow (mBd (), 2) - pow (mK (), 2))*
                  ml*fz (qsq))/(qsq*
                  pow (fp (qsq),
                   2)))*(((CPIm () + CPRHIm ())*(pow (mBd (), 2) -
                    pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  fp (qsq)) +
               C10Im ()*
                ml*(-1 + ((pow (mBd (), 2) - pow (mK (), 2))*(-1 +
                    fz (qsq)/fp (qsq)))/qsq)) +
            2*(-((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                    pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  pow (fp (qsq),
                   2)) - (C10Re ()*(pow (mBd (), 2) - pow (mK (), 2))*
                  ml*fz (qsq))/(qsq*
                  pow (fp (qsq),
                   2)))*(((CPRe () + CPRHRe ())*(pow (mBd (), 2) -
                    pow (mK (), 2))*fz (qsq))/(2.*(mb () - ms ())*
                  fp (qsq)) +
               C10Re ()*
                ml*(-1 + ((pow (mBd (), 2) - pow (mK (), 2))*(-1 +
                    fz (qsq)/fp (qsq)))/qsq))) + (lambda (qsq)*((-16*
                 CTIm ()*ml*
                 fT (qsq)*(C9Im () +
                   C9RHIm () + (8*CTIm ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPIm (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*pow (fp (qsq), 2)) - (16*CTRe ()*ml*
                 fT (qsq)*(C9Re () +
                   C9RHRe () + (8*CTRe ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPRe (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*pow (fp (qsq), 2))))/4.) + (nf (qsq, ml)*
         pow (xiP (qsq),
          2)*(qsq*((-8*pow (CT5Im (), 2)*pow (betal (qsq, ml), 2)*
                 pow (fT (qsq), 2)*
                 lambda (qsq))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 3)) - (8*pow (CT5Re (), 2)*
                 pow (betal (qsq, ml), 2)*pow (fT (qsq), 2)*
                 lambda (qsq))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 3)) +
              pow (betal (qsq, ml),
                2)*((-8*pow (CTIm (), 2)*pow (betal (qsq, ml), 2)*
                    pow (fT (qsq), 2)*
                    lambda (qsq))/(pow (mBd () + mK (), 2)*
                    pow (fp (qsq), 3)) - (8*pow (CTRe (), 2)*
                    pow (betal (qsq, ml), 2)*pow (fT (qsq), 2)*
                    lambda (qsq))/(pow (mBd () + mK (), 2)*
                    pow (fp (qsq), 3)))) - (pow (betal (qsq, ml), 2)*
              lambda (qsq)*((-16*CTIm ()*ml*
                   fT (qsq)*(C9Im () +
                    C9RHIm () + (8*CTIm ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPIm (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                    mK ())*pow (fp (qsq), 2)) - (16*CTRe ()*ml*
                   fT (qsq)*(C9Re () +
                    C9RHRe () + (8*CTRe ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPRe (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                    mK ())*pow (fp (qsq), 2))))/4. +
           2*ml*betal (qsq, ml)*
            sqrt (lambda (qsq))*((-16*pow (CTIm (), 2)*ml*
                 betal (qsq, ml)*pow (fT (qsq), 2)*
                 sqrt (lambda (qsq)))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 3)) - (16*pow (CTRe (), 2)*ml*
                 betal (qsq, ml)*pow (fT (qsq), 2)*
                 sqrt (lambda (qsq)))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 3)) - (2*CTIm ()*betal (qsq, ml)*
                 fT (qsq)*
                 sqrt (lambda (qsq))*(C9Im () +
                   C9RHIm () + (8*CTIm ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPIm (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*pow (fp (qsq), 2)) - (2*CTRe ()*
                 betal (qsq, ml)*fT (qsq)*
                 sqrt (lambda (qsq))*(C9Re () +
                   C9RHRe () + (8*CTRe ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPRe (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*pow (fp (qsq), 2)))))/3., 2) +
   4*pow (ErfT (qsq), 2)*
    pow ((nf (qsq, ml)*
         pow (xiP (qsq),
          2)*(qsq*((8*pow (CT5Im (), 2)*pow (betal (qsq, ml), 2)*
                 fT (qsq)*lambda (qsq))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 2)) + (8*pow (CT5Re (), 2)*
                 pow (betal (qsq, ml), 2)*fT (qsq)*
                 lambda (qsq))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 2)) +
              pow (betal (qsq, ml),
                2)*((8*pow (CTIm (), 2)*pow (betal (qsq, ml), 2)*
                    fT (qsq)*lambda (qsq))/(pow (mBd () + mK (), 2)*
                    pow (fp (qsq), 2)) + (8*pow (CTRe (), 2)*
                    pow (betal (qsq, ml), 2)*fT (qsq)*
                    lambda (qsq))/(pow (mBd () + mK (), 2)*
                    pow (fp (qsq), 2)))) - (pow (betal (qsq, ml), 2)*
              lambda (qsq)*((16*CTIm ()*
                   ml*(C9Im () +
                    C9RHIm () + (8*CTIm ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPIm (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                    mK ())*fp (qsq)) + (16*CTRe ()*
                   ml*(C9Re () +
                    C9RHRe () + (8*CTRe ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPRe (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                    mK ())*fp (qsq))))/4. +
           2*ml*betal (qsq, ml)*
            sqrt (lambda (qsq))*((16*pow (CTIm (), 2)*ml*
                 betal (qsq, ml)*fT (qsq)*
                 sqrt (lambda (qsq)))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 2)) + (16*pow (CTRe (), 2)*ml*
                 betal (qsq, ml)*fT (qsq)*
                 sqrt (lambda (qsq)))/(pow (mBd () + mK (), 2)*
                 pow (fp (qsq), 2)) + (2*CTIm ()*betal (qsq, ml)*
                 sqrt (lambda (qsq))*(C9Im () +
                   C9RHIm () + (8*CTIm ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPIm (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*fp (qsq)) + (2*CTRe ()*betal (qsq, ml)*
                 sqrt (lambda (qsq))*(C9Re () +
                   C9RHRe () + (8*CTRe ()*ml*
                    fT (qsq))/((mBd () + mK ())*fp (qsq)) + (2*mb ()*
                    TauPRe (qsq))/(mBd ()*xiP (qsq))))/((mBd () +
                   mK ())*fp (qsq)))))/
       3. + (lambda (qsq)*nf (qsq, ml)*
         pow (xiP (qsq),
          2)*((16*CTIm ()*
              ml*(C9Im () +
                C9RHIm () + (8*CTIm ()*ml*fT (qsq))/((mBd () + mK ())*
                   fp (qsq)) + (2*mb ()*TauPIm (qsq))/(mBd ()*
                   xiP (qsq))))/((mBd () + mK ())*fp (qsq)) + (16*
              CTRe ()*
              ml*(C9Re () +
                C9RHRe () + (8*CTRe ()*ml*fT (qsq))/((mBd () + mK ())*
                   fp (qsq)) + (2*mb ()*TauPRe (qsq))/(mBd ()*
                   xiP (qsq))))/((mBd () + mK ())*fp (qsq))))/4., 2)); }

/*double BdtoKll_obserr::diffBrnch(double qsq, double ml){
    return tauBd()/2.0*diffWidth(qsq,ml); }

double BdtoKll_obserr::diffAFB(double qsq, double ml){
    return blNP(qsq,ml)/diffWidth(qsq,ml); }

double BdtoKll_obserr::diffFH(double qsq, double ml){
    return 2.0*(alNP(qsq,ml)+clNP(qsq,ml))/diffWidth(qsq,ml); }*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

