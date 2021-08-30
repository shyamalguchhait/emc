#include <complex>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
typedef std::complex<float> fcmplx;
#define NMOMLEN 20201 /* =2*MAXTRM+1, PMON is of size [4][NMOMLEN] */
using namespace std;
/*
     You need special compiler options to promote real to double precision etc !!!

      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE ) 
      LOGICAL  ANYANG, PERFCT, PRNT(*)
      INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
      REAL     GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,
     &         XMU(*), XX
      COMPLEX  CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
*/
extern "C"
{
    void miev0_(float *XX, fcmplx *CREFIN, long int *PERFCT, float *MIMCUT, long int *ANYANG, int *NUMANG, float *XMU, int *NMOM, int *IPOLZN, long int *MOMDIM, long int *PRNT, float *QEXT, float *QSCA, float *GQSC, float *PMOM, fcmplx *SFORW, fcmplx *SBACK, fcmplx *S1, fcmplx *S2, fcmplx *TFORW, fcmplx *TBACK, float *SPIKE);
};
int main(void)
{
    long int MAXANG, MOMDIM;
    MAXANG = 7;
    MOMDIM = 200;
    long int ANYANG, PERFCT, PRNT[2];
    int NUMANG, NMOM, IPOLZN;
    float GQSC, MIMCUT, QEXT, QSCA, SPIKE, PMOM[MOMDIM][4], XMU[MAXANG], XX;
    fcmplx CREFIN, SFORW, SBACK, S1[MAXANG], S2[MAXANG], TFORW[2], TBACK[2];

    // LOCAL VARIABLES
    long int NOPMOM;
    int i, j, k, NCAS, NPQUAN, time0, time1, cntrat, maxcnt;
    float DEGPOL, FNORM, I1, I2, INTEN, PI, QABS, TESTIN;

    //local array
    float ANGLE[MAXANG];

    int NCASES;
    NCASES = 19;

    // Arrays in common
    long int TESTAN[NCASES];
    int TESTIP[NCASES];
    float ESTXX[NCASES], TESTQE[NCASES], TESTQS[NCASES], TESTGQ[NCASES], TESTPM[MOMDIM][4][NCASES];

    PI = M_PI;
    NOPMOM = 0;

    MIMCUT = 1e-6;
    NUMANG = MAXANG;
    for (i = 0; i <= NUMANG; i++)
    {
        ANGLE[i] = (i)*180 / (NUMANG - 1);
        XMU[i] = cos((PI / 180) * ANGLE[i]);
        //printf("%f\t", XMU[i]);
        printf("%d\t", i);
    }

    XX = .001;
    CREFIN = fcmplx(1., 0.);
    PERFCT = 0;
    ANYANG = 1;
    IPOLZN = +1234;
    NMOM = 1;

    if (NOPMOM == 1)
    {
        NMOM = 0;
    }

    PRNT[1] = 1;
    PRNT[2] = 1;
    //printf("%lf + i %lf", creal(CREFIN), cimag(CREFIN));
    //printf("%ld\n", MOMDIM);
    miev0_(&XX, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, PMOM[4], &SFORW, &SBACK, S1, S2, TFORW, TBACK, &SPIKE);
    //printf("%ld", sizeof(S1));
    for (i = 0; i <= NUMANG - 1; i++)
    {
        cout << ANGLE[i] << S1[i].real() << S1[i].imag() << endl;
    }
    //PRNT[1] = 0;
    //PRNT[2] = 0;

    ofstream f("data.dat");
    int N = 6000;
    float n = 4, r = .1, l[N], xx[N], q[N];
    CREFIN = fcmplx(n, 0);
    float x;
    for (i = 0; i <= N; i++)
    {
        l[i] = .4 + (.0001) * i;
        xx[i] = 2 * PI * r / l[i];
        //printf("%f\t", l[i]);
        x = xx[i];
        miev0_(&x, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, PMOM[4], &SFORW, &SBACK, S1, S2, TFORW, TBACK, &SPIKE);
        q[i] = QEXT;
        f << l[i] << "\t" << xx[i] << "\t" << q[i] << endl;
    }
    f.close();

    return 0;
}