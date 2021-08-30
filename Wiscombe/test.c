#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

void miev0_(float *XX, float complex *CREFIN, long int *PERFCT, float *MIMCUT, long int *ANYANG, int *NUMANG, float *XMU, int *NMOM, int *IPOLZN, long int *MOMDIM, long int *PRNT, float *QEXT, float *QSCA, float *GQSC, float *PMOM, float complex *SFORW, float complex *SBACK, float complex *S1, float complex *S2, float complex *TFORW, float complex *TBACK, float *SPIKE);

int main()
{
    long int MAXANG, MOMDIM;
    MAXANG = 7;
    MOMDIM = 200;
    long int ANYANG, PERFCT, PRNT[2];
    int NUMANG, NMOM, IPOLZN;
    float GQSC, MIMCUT, QEXT, QSCA, SPIKE, PMOM[MOMDIM][4], XMU[MAXANG], XX;
    float complex CREFIN, SFORW, SBACK, S1[MAXANG], S2[MAXANG], TFORW[2], TBACK[2];

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
    CREFIN = CMPLX(1., 0.);
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
        printf("%f \t %e + %e i\n", ANGLE[i], creal(S1[i]), cimag(S1[i]));
    }
    //PRNT[1] = 0;
    //PRNT[2] = 0;

    FILE *f = fopen("data.dat", "wb");
    int N = 6000;
    float n = 4, r = .1, l[N], xx[N], q[N];
    CREFIN = CMPLX(n, 0);
    float x;
    for (i = 0; i <= N; i++)
    {
        l[i] = .4 + (.0001) * i;
        xx[i] = 2 * PI * r / l[i];
        //printf("%f\t", l[i]);
        x = xx[i];
        miev0_(&x, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, &QEXT, &QSCA, &GQSC, PMOM[4], &SFORW, &SBACK, S1, S2, TFORW, TBACK, &SPIKE);
        q[i] = QEXT;
        fprintf(f, "%e\t%e\t%e\n", l[i], xx[i], q[i]);
    }
    fclose(f);

    return 0;
}
