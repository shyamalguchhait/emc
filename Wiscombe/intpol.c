#include <math.h>
#include <stdio.h>

double intp(double arrx[], double arry[], double x, int len)
{
    int i;
    double min, max;
    //int len = sizeof(arrx) / sizeof(double);
    double near[len], nearlow, nearhigh;
    double xa, xb;
    int ia, ib;
    double y;
    for (i = 0; i < len; i++)
    {
        //printf("%lf\t", arrx[i]);
    }

    if (len < 0)
    {
        return -1;
    }
    else
    {
        min = arrx[0];
        max = arrx[0];
        for (i = 1; i < len; i++)
        {
            if (arrx[i] < min)
            {
                min = arrx[i];
            }
            if (arrx[i] > max)
            {
                max = arrx[i];
            }
        }
        //printf("%lf\t%lf", min, max);

        if (x > min && x < max)
        {
            for (i = 0; i < len; i++)
            {
                near[i] = arrx[i] - x;
                //printf("%lf\t", near[i]);
            }
            for (i = 0; i <= len; i++)
            {
                if (near[i] == 0)
                {
                    //printf("%lf\t", near[i]);
                    return arry[i];
                }
            }
            nearlow = near[0];
            nearhigh = near[len - 1];
            //printf("nearlow = %lf\tnearhigh = %lf", nearlow, nearhigh);
            for (i = 1; i < len; i++)
            {

                if (near[i] < 0 && near[i] > nearlow)
                {
                    nearlow = near[i];
                    ia = i;
                }
            }
            for (i = len - 2; i >= 0; i--)
            {
                if (near[i] > 0 && near[i] < nearhigh)
                {
                    nearhigh = near[i];
                    ib = i;
                }
            }
        }
        //printf("nearlow = %lf\tnearhigh = %lf", nearlow, nearhigh); /*
        y = arry[ib] + (arry[ib] - arry[ia]) * (x - arrx[ia]) / (arrx[ib] - arrx[ia]);
        return y;
    }
}

int main()
{
    int i;
    double arrx[] = {1, 2, 3, 4}, arry[] = {1, 0, .6, 1};
    double x = 2.5, y;
    int len = sizeof(arrx) / (sizeof(double));
    for (i = 0; i < 4; i++)
    {
        //   printf("%lf\t", arrx[i]);
    }
    y = intp(arrx, arry, x, len);
    printf("%lf\n", y);
    FILE *fo, fp;
    fo = fopen("Johnson.dat", "rb");

    i = 0;
    double w[100][3];
    while (true)
    {
        fgets(fo, "%d", w[i][0]);
        //printf("%lf\n", w[i]);
        i++;
    }

    return 0;
}