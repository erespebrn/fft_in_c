void rootsofunity(double *realtwiddle, double *imtwiddle, unsigned int logN)
{
    double twopi = 6.2831853071;
    unsigned int n, N;

    N = 1<<logN;

    for(n=0; n<(N>>1); n++)
    {
       realtwiddle[n] = cos(twopi*n/N);
       imtwiddle[n] = -sin(twopi*n/N);

    }
}

void fft(double *real, double *im, double *realtwiddle, double *imtwiddle, unsigned int N)
{
    unsigned int even, odd = 0, span, log, rootindex;
    double temp;
    log=0;

    for(span=N>>1; span; span>>=1, log++)
    {
        qDebug() << "";
        for(odd=span; odd<N; odd++)
        {
            odd = odd | span;
            even = odd ^ span;


            temp = real[even] + real[odd];
            real[odd] = real[even] - real[odd];
            real[even] = temp;

            temp = im[even] + im[odd];
            im[odd] = im[even] - im[odd];
            im[even] = temp;

            rootindex = (even<<log) & (N-1);
            temp = realtwiddle[rootindex] * real[odd] - imtwiddle[rootindex] * im[odd];
            im[odd] = realtwiddle[rootindex] * im[odd] + imtwiddle[rootindex] * real[odd];
            real[odd] = temp;
        }
     }
}

static inline void swap(unsigned int forward, unsigned int rev, double *real, double *im)
{
    double temp;

    temp = real[forward];
    real[forward] = real[rev];
    real[rev] = temp;

    temp = im[forward];
    im[forward] = im[rev];
    im[rev] = temp;
}

void complexbitrev(double *real, double *im, unsigned int logN)
{
    unsigned int i, forward, rev, zeros;
    unsigned int nodd, noddrev, N;
    unsigned int halfn, quartn, nmin1;

    N = 1<<logN;
    halfn = N>>1;
    quartn = N>>2;
    nmin1 = N-1;

    forward = halfn;
    rev = 1;

    for(i=quartn; i; i--)
    {

     nodd = ~i;                                  
     for(zeros=0; nodd&1; zeros++) nodd >>= 1;   
     forward ^= 2 << zeros;
     rev ^= quartn >> zeros;

        if(forward<rev)
        {
            swap(forward, rev, real, im);
            nodd = nmin1 ^ forward;
            noddrev = nmin1 ^ rev;
            swap(nodd, noddrev, real, im);
        }

        nodd = forward ^ 1;
        noddrev = rev ^ halfn;
        swap(nodd, noddrev, real, im);
    }
}
