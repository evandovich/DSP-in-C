/////////////////////////////////////////////////////////////////////////////////////////////////
// File: MINILITE2.0V_Cprogram.c
//
// Author           : Evans Baidoo
// Email            : ebaidoo2@hhu.edu.cn
// Create Date      : 2021.08.20
// Module Name      : Signal_statistics
// Project Name     : miniDSP Xhorse
// Target Device    :
// Description      : This program computes the statistics of a continuous-Time/Analog signal. Signal can be real or
//                  : complex communication signal obtained from time domain devices such as Oscilloscope, Digitizer
//                  : or audio capture devices.
//                  : Notable computation includes: MEAN, STANDARD DEVIATION, VARIANCE & SNR
//
// Revision         : V2.0
// Additional Comments:
//////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX_FILE_NAME 25


double calc_signal_mean(double *signal_source_arr, int sig_length); // define mean function
double calc_signal_variance(double * signal_source_arr, double sig_mean, int sig_length);   // define variance function
double calc_signal_std(double sig_variance);                        // define signal standard deviation
double calc_signal_snr(double sig_mean, double sig_std);             // define signal to noise ratio
double calc_signal_skewness(double *signal_source_arr, double sig_mean, double std_deviation,  int sig_length);
double calc_signal_kurtosis(double *signal_source_arr, double sig_mean, int sig_length);
double INPUT_SIG_ARR[5000];
double MEAN;
double VARIANCE;
double STD;
double SNR;
double SKEWNESS;
double KURTOSIS;



int main()
{
    int SIG_LENGTH = 0;
    char filename[MAX_FILE_NAME];
    double num;

    FILE *fptr;
    /* Get file name from user. The file should be either in current folder.
       In case file is in a different destination, complete path should be provided */
    printf("Enter file name: ");
    fgets(filename, 25, stdin);
    filename[strlen(filename)-1] = '\0';

    /* Open the file and Check if file exists  */
    if ((fptr = fopen(filename, "r")) == NULL)          //run a statement until a condition is met or true
      {
        printf("\nCould not locate file %s.", filename);
        printf("\nMake sure the file name and destination address is correct");

        return 0;
      }

    printf("File %s loaded successfully\n ",filename);

    while (fscanf(fptr, "%lf",&num) != EOF)
    {

        INPUT_SIG_ARR[SIG_LENGTH] = num;
        SIG_LENGTH +=1;
    }


    MEAN = calc_signal_mean(&INPUT_SIG_ARR[0], SIG_LENGTH);
    VARIANCE = calc_signal_variance(&INPUT_SIG_ARR[0],MEAN,SIG_LENGTH);
    STD = calc_signal_std(VARIANCE);
    SNR = calc_signal_snr(MEAN, STD);
    SKEWNESS = calc_signal_skewness(&INPUT_SIG_ARR[0], MEAN, STD, SIG_LENGTH);
    KURTOSIS = calc_signal_kurtosis(&INPUT_SIG_ARR[0], MEAN, SIG_LENGTH);

    printf("\nMEAN : %f",MEAN);
    printf("\nVARIANCE : %f",VARIANCE);
    printf("\nSTANDARD DEVIATION : %f",STD);
    printf("\nSNR (dB) : %f",SNR);
    printf("\nSKEWNESS : %f",SKEWNESS);
    printf("\nKURTOSIS : %f",KURTOSIS);

    fclose(fptr);

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_mean
// Description:
//              This function obtains a data signal and performs the mean of the signal
//
//////////////////////////////////////////////////////////////////////////////////////////////////

double calc_signal_mean(double *signal_source_arr, int sig_length)
{
    double _mean = 0.0;
    for (int i =0;i<sig_length;i++)
    {
        _mean = _mean + signal_source_arr[i];

    }
    _mean = _mean/(double)sig_length;
    return _mean;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_variance
// Description:
//              This function obtains a data signal and performs the variance of the signal
//
//////////////////////////////////////////////////////////////////////////////////////////////////
double calc_signal_variance(double *signal_source_arr, double sig_mean, int sig_length)
{
    double _variance = 0.0;
    double  stan_dev = 0.0;
    for (int i =0; i<sig_length;i++)
    {
        stan_dev = signal_source_arr[i]-sig_mean;
        _variance = _variance +pow(stan_dev,2);
    }
    _variance =_variance /(sig_length-1);
    return _variance;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_std
// Description:
//              This function obtains a data signal and performs the standard deviation of the signal
//
//////////////////////////////////////////////////////////////////////////////////////////////////////
double calc_signal_std(double sig_variance)
{
    double _std = sqrt(sig_variance);
    return _std;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_snr
// Description:
//              This function obtains a data signal and performs the signal-to-noise ratio of the signal
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
double calc_signal_snr(double sig_mean, double sig_std)
{
     double _snr = 10*log10(sig_mean/sig_std);                      // log of negative product factor is invalid.
     return _snr;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_skewness
// Description:
//              This function obtains a data signal and performs the skewness of the signal
//              SKEWNESS  = measure the symmetry of the signal distribution about its mean
//              The signal skewness is approximately asymmetrical distribution if it > 0
//              The signal skewness is approximately symmetrical distribution if it < 0
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_signal_skewness(double *signal_source_arr, double sig_mean, double std_deviation,  int sig_length)
{
    double _skewness =  0.0;
    double  stan_ = 0.0;
    double  stan_num = 0.0;

    for (int i = 0; i<sig_length ; i++)
    {
        stan_ = signal_source_arr[i]-sig_mean;
        stan_num = stan_num + pow(stan_,3);
    }
    _skewness = stan_num /((sig_length-1)*pow(std_deviation,3));
    return _skewness;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function: calc_signal_kurtosis
// Description:
//              This function obtains a data signal and performs the skewness of the signal
//              KURTOSIS  = measure whether the tail of distribution extends beyond the  ±3
//              standard deviation of the mean or not.
//              The signal kurtosis distribution is heavy-tailed  if it > 3
//              The signal kurtosis is short-tailed distribution if it < 3
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_signal_kurtosis(double *signal_source_arr, double sig_mean, int sig_length)
{
    double _kurtosis =  0.0;
    double  second_moment = 0.0;
    double  fourth_moment = 0.0;
    double second_ = 0.0;   double fourth_ = 0.0;

    for (int i = 0; i<sig_length ; i++)
    {
        fourth_= signal_source_arr[i]-sig_mean;
        fourth_moment = fourth_moment + pow(fourth_,4);

        second_ = signal_source_arr[i]- pow(sig_mean,2);
        second_moment = second_moment + pow(second_,2);
    }
    _kurtosis = fourth_moment / second_moment ;
    return _kurtosis;

}
