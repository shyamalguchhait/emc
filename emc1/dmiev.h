#include <complex>
#include <assert.h>
typedef std::complex<float> fcmplx;
#define NMOMLEN 20201 /* =2*MAXTRM+1, PMON is of size [4][NMOMLEN] */

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

/* calls miev0 with one given mu=cos(theta) */
void amiev(float xx, fcmplx crefin, float mu, fcmplx *s1, fcmplx *s2)
{
    long int PERFCT = 0;
    float MIMCUT = 1e-8;
    long int ANYANG = 1;
    int NUMANG = 1;
    int NMOM = 0;
    int IPOLZN = 0;
    long int MOMDIM = 1;
    long int PRNT[2] = {0, 0};
    float QEXT, QSCA, GQSC, PMOM, SPIKE;
    fcmplx SFORW, SBACK, TFORW[2], TBACK[2];
    miev0_(&xx, &crefin, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, &mu, &NMOM, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, &PMOM, &SFORW, &SBACK, s1, s2, TFORW, TBACK, &SPIKE);
};

/* calls miev0 with one given mu=cos(theta) and obtain Legendre coefficients of S21 */
/* pmon is a pointer to array pmom[4][nmom]. numang must be odd */
void mievp(float xx, fcmplx crefin, int numang, fcmplx *s1, fcmplx *s2, float pmom[4][NMOMLEN])
{
    long int PERFCT = 0;
    float MIMCUT = 1e-8;
    long int ANYANG = 0;
    int IPOLZN = 3; /* calculates S21 */
    long int MOMDIM = NMOMLEN - 1;
    long int PRNT[2] = {0, 0};
    float QEXT, QSCA, GQSC, SPIKE;
    fcmplx SFORW, SBACK, TFORW[2], TBACK[2];
    int nmom = (int)floor(2. * (xx + 4. * pow(xx, 1. / 3.) + 2.)); /* counting from 0 to nmom inclusively in Fortran */
    float step = 2. / (numang - 1);
    float *mu = new float[numang];
    for (int i = 0; i < numang; i++)
        mu[i] = 1. - i * step;
    assert(nmom > 2 * (xx + 4 * pow(xx, 1 / 3.) + 1));
    miev0_(&xx, &crefin, &PERFCT, &MIMCUT, &ANYANG, &numang, mu, &nmom, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, &pmom[0][0], &SFORW, &SBACK, s1, s2, TFORW, TBACK, &SPIKE);
    delete[] mu;
};

/* calls miev0 with n given mu=cos(theta) */
void miev(float xx, fcmplx crefin, int numang, float *mu, fcmplx *s1, fcmplx *s2)
{
    long int PERFCT = 0;
    float MIMCUT = 1e-8;
    long int ANYANG = 1;
    int NMOM = 0;
    int IPOLZN = 0;
    long int MOMDIM = 1;
    long int PRNT[2] = {0, 0};
    float QEXT, QSCA, GQSC, PMOM, SPIKE;
    fcmplx SFORW, SBACK, TFORW[2], TBACK[2];
    miev0_(&xx, &crefin, &PERFCT, &MIMCUT, &ANYANG, &numang, mu, &NMOM, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, &PMOM, &SFORW, &SBACK, s1, s2, TFORW, TBACK, &SPIKE);
};

/* obtain the basic info */
void mievinfo(float xx, fcmplx crefin, float *qext, float *qsca, float *g)
{
    long int PERFCT = 0;
    float MIMCUT = 1e-8;
    long int ANYANG = 1;
    int NUMANG = 1;
    float XMU = 0;
    int NMOM = 0;
    int IPOLZN = 0;
    long int MOMDIM = 1;
    long int PRNT[2] = {0, 0};
    float gqsc, PMOM, SPIKE;
    fcmplx SFORW, SBACK, S1, S2, TFORW[2], TBACK[2];
    miev0_(&xx, &crefin, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, &XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, qext, qsca, &gqsc, &PMOM, &SFORW, &SBACK, &S1, &S2, TFORW, TBACK, &SPIKE);
    *g = gqsc / (*qsca);
};
