/***************************************************************************
 *   Copyright (C) 2015 by Petr Voytsik                                    *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3.h>
#include <mark5access.h>


struct fft_data_t {
    fftwf_plan *plan;
    fftwf_complex **zdata;
    float **data;
    double **spec;
};

static void usage(const char *prog_name)
{
    printf("Usage: %s file_name format nchan aver_time(ms) total_time(s)\n\n", 
            prog_name);
    printf("file_name  - the name of the input file\n");
    printf("format     - mark5access data format in form <FORMAT>-<Mbps>-<nchan>-<nbit>\n");
    printf("nchan      - number of spectral channels\n");
    printf("aver_time  - approximate integration time per spectrum in milliseconds\n");
    printf("total_time - total time in seconds\n");
}

static void fft_data_init(struct fft_data_t *fft_data, unsigned nchan)
{
    unsigned i;

    fft_data->spec = (double **)malloc(nchan * sizeof(double *));
    fft_data->data = (float **)malloc(nchan * sizeof(float *));
    fft_data->zdata = (fftwf_complex **)malloc(nchan * sizeof(fftwf_complex *));
    fft_data->plan = (fftwf_plan *)malloc(nchan * sizeof(fftwf_plan));

    for(i = 0; i < nchan; ++i){
        fft_data->spec[i] = (double *)malloc(nchan * sizeof(double));
        fft_data->zdata[i] = fftwf_alloc_complex(nchan+2);
        /* data[i] = (float *)calloc(2*nchan+2, sizeof(float)); */
        /* Use in-place FFT */
        fft_data->data[i] = (float *)(fft_data->zdata[i]);
        fft_data->plan[i] = fftwf_plan_dft_r2c_1d(nchan*2, fft_data->data[i], 
                                                  fft_data->zdata[i], FFTW_MEASURE);
    }

}

static void fft_data_free(struct fft_data_t *fft_data, unsigned nchan)
{
    unsigned i;

    for(i = 0; i < nchan; ++i){
        fftwf_destroy_plan(fft_data->plan[i]);
        free(fft_data->spec[i]);
        fftwf_free(fft_data->zdata[i]);
    }

    free(fft_data->plan);
    free(fft_data->spec);
    free(fft_data->data);
    free(fft_data->zdata);

}

static int spec(const char *filename, const char *format, int nchan, 
                double aver_time, double total_time, int if_num)
{
    struct mark5_stream *ms;
    struct fft_data_t fft_data;
    long long offset = 0;
    int nint;
    double real_step;
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
    fft_data_init(&fft_data, ms->nchan);

    for(k = 0; k < step_num; ++k){
        /* Zero spec */
        /*for(i = 0; i < ms->nchan; ++i)*/
        i = if_num;
            for(c = 0; c < nchan; ++c)
                fft_data.spec[i][c] = 0.0;

        for(j = 0; j < nint; ++j){
            int status;
            double re, im;
            
            status = mark5_stream_decode(ms, 2*nchan, fft_data.data);
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
                fftwf_execute(fft_data.plan[i]);
                for(c = 0; c < nchan; ++c){

                    re = creal(fft_data.zdata[i][c]);
                    im = cimag(fft_data.zdata[i][c]);
                    fft_data.spec[i][c] += (re*re + im*im) / (double)(2*nchan);
                }
            /*}*/
        }

        for(c = 0; c < nchan; ++c){
            printf("%lf ", fft_data.spec[if_num][c] / (double)nint);
        }
        putchar('\n');

    }

    fft_data_free(&fft_data, ms->nchan);
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
