#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3.h>
#include <mark5access.h>


static void usage(const char *prog_name)
{
    printf("Usage: %s file_name format nchan aver_time(ms) total_time(s)\n", prog_name);
}

static int spec(const char *filename, const char *format, int nchan, 
                double aver_time, double total_time, int if_num)
{
    struct mark5_stream *ms;
    long long offset = 0;
    int nint;
    double real_step;
    fftwf_plan *plan;
    fftwf_complex **zdata;
    float **data;
    double **spec;
    int i, j, k;
    int c;  /* Iteration over spectral channels */
    int step_num = 0;  /* Number of time steps */

    ms = new_mark5_stream_absorb(new_mark5_stream_file(filename, offset),
                                 new_mark5_format_generic_from_string(format));
    
    if(!ms){
        fprintf(stderr, "Error: problem opening %s\n", filename);

        return EXIT_FAILURE;
    }

    /*mark5_stream_print(ms);*/

    /* aver_time in ms */
    nint = (aver_time * 1e-3) * ms->samprate / (2 * nchan);
    fprintf(stderr, "nint = %d\n", nint);
    real_step = (double)(nint * 2 * nchan) / (double)ms->samprate;
    fprintf(stderr, "Real time step = %lf ms\n", real_step * 1e3);
    step_num = (int)(total_time / real_step);
    fprintf(stderr, "Number of time steps = %d\n", step_num);

    /* Prepare data arrays */
    spec = (double **)malloc(ms->nchan * sizeof(double *));
    data = (float **)malloc(ms->nchan * sizeof(float *));
    zdata = (fftwf_complex **)malloc(ms->nchan * sizeof(fftwf_complex *));
    plan = (fftwf_plan *)malloc(ms->nchan * sizeof(fftwf_plan));

    for(i = 0; i < ms->nchan; ++i){
        spec[i] = (double *)malloc(nchan * sizeof(double));
        zdata[i] = fftwf_alloc_complex(nchan+2);
        /* data[i] = (float *)calloc(2*nchan+2, sizeof(float)); */
        /* Use in-place FFT */
        data[i] = (float *)zdata[i];
        plan[i] = fftwf_plan_dft_r2c_1d(nchan*2, data[i], zdata[i], FFTW_MEASURE);
    }

    for(k = 0; k < step_num; ++k){
        /* Zero spec */
        /*for(i = 0; i < ms->nchan; ++i)*/
        i = if_num;
            for(c = 0; c < nchan; ++c)
                spec[i][c] = 0.0;

        for(j = 0; j < nint; ++j){
            int status;
            double re, im;
            
            status = mark5_stream_decode(ms, 2*nchan, data);
            if(status < 0){
                fprintf(stderr, "Error: mark5_stream_decode failed\n");
                break;
            }
            if(ms->consecutivefails > 5){
                fprintf(stderr, "Error: problem with data decoding\n");
                break;
            }
            /*for(i = 0; i < ms->nchan; ++i){*/
            i = if_num;
                fftwf_execute(plan[i]);
                for(c = 0; c < nchan; ++c){

                    re = creal(zdata[i][c]);
                    im = cimag(zdata[i][c]);
                    spec[i][c] += (re*re + im*im) / (double)(2*nchan);
                }
            /*}*/
        }

        for(c = 0; c < nchan; ++c){
            printf("%lf ", spec[if_num][c] / (double)nint);
        }
        putchar('\n');

    }

    for(i = 0; i < ms->nchan; ++i){
        fftwf_destroy_plan(plan[i]);
        free(spec[i]);
        fftwf_free(zdata[i]);
    }

    free(plan);
    free(spec);
    free(data);
    free(zdata);

    delete_mark5_stream(ms);

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
    int ret = 0;
    int nchan;
    double aver_time, total_time;

    if(argc != 6){
        usage(argv[0]);
        exit(0);
    }

    /* Check input parameters */
    nchan = atoi(argv[3]);
    if(nchan <=0 || nchan > 2<<15){
        fprintf(stderr, "Ivalid number of spectral channels: %d\n", nchan);

        exit(EXIT_FAILURE);
    }
    aver_time = atof(argv[4]);
    if(aver_time <= 0){
        fprintf(stderr, "Invalid average time: %lf\n", aver_time);

        exit(EXIT_FAILURE);
    }
    total_time = atof(argv[5]);
    if(total_time <= 0){
        fprintf(stderr, "Invalid total time: %lf\n", total_time);

        exit(EXIT_FAILURE);
    }

    ret = spec(argv[1], argv[2], nchan, aver_time, total_time, 1);

    return ret;
}
