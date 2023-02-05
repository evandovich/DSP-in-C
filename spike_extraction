
/////////////////////////////////////////////////////////////////////////////////////////////////
// File: SPIKR.0V_C++program.cpp
//
// Author           : Evans Baidoo
// Email            : ebaidoo2@hhu.edu.cn
// Create Date      : 2022.08.20
// Module Name      : Spike_Extration
// Project Name     : SpecTRUM
// Target Device    :
// Description      : This program performs spike extraction from a data sample collected from a sensor device
//                  : The signal is averaged and downsampled, and a window of the signal is analyse in 
//                  : in the Frequency spectrum. The higher spikes are then collected together with their freqs
//                  :
//
// Revision         : V2.0
// Additional Comments:
//////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#define SIZE 100
using namespace std;

// function to compute FFT
void dft_signal_calc(double *signal_source_arr, double *sig_output_real_arr, double *sig_output_imx_arr, int sig_Length )
{
    int i,j,k;
    double PI = 3.14159265359;
    for(j = 0; j<sig_Length/2; j++)
    {
        sig_output_real_arr[j] =0;
        sig_output_imx_arr[j] =0;

    }
    // Main DFT algorithm
    for(k = 0;k<sig_Length/2;k++)
    {
        for(i=0;i<sig_Length;i++)
        {
            sig_output_real_arr[k] = sig_output_real_arr[k] + signal_source_arr[i]*cos(2*PI*k*i/sig_Length);
            sig_output_imx_arr[k] = sig_output_imx_arr[k] - signal_source_arr[i]*sin(2*PI*k*i/sig_Length);
        }
    }
}

 //* Divide one integer by another, rounding towards minus infinity.
int div_floor(int x, int y) {
    int q = x/y;
    int r = x%y;
    if ((r!=0) && ((r<0) != (y<0))) --q;
    return q;
}

// m-point averaged and downsampled
void downsample(double *signal_source_arr, double *wav_avg, int samp_len, int bin)
{
    int j;
    for (int i = 1; i<=samp_len; i++)
     {
        double sum = 0;
         for (j = bin*(i-1)+1;j<bin*i;j++){
               sum += signal_source_arr[j];
         }
          wav_avg[i] = sum/j;
         }
}

//calculate median
        //void median(double *arr, int size){
        double median(double arr[], int size){
        sort(arr, arr+size);
        if (size % 2 != 0)
        {
         return (double)arr[size/2];
         //res = arr[size/2];
        }
        else{
            //res = (arr[(size-1)/2] + arr[size/2])/2.0;
            return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
        }

        }

// spike finding
    void spike_func(double *y_threshold, double *spike_hold_arr, int ii, int pp )
    {
        double thresh, abs_val[0];
        while (ii <= pp)
    {
        int maxx = max(1,ii-450);
        for (int j = maxx; j<=ii ; j++)
        {
            abs_val[j] = fabs(y_threshold[j]);
          //  cout <<"yval hold :" <<abs_val[j]<<'\n';
        }

       double midl = median(abs_val, ii);
        thresh = 1*midl;
        if(fabs(y_threshold[ii])> thresh)
        {
           for(int j = 0; j<= 20;j++)
        {
            if (ii+j <= pp)
            {
                spike_hold_arr[ii] = y_threshold[ii];
            //  cout <<"spike :" <<spike_hold_arr[ii] <<'\n';
                ii += 1;
            }
        }
        ii = ii + 80;
    }
    ii += 1;
    }
    }



    // Create a vector of evenly spaced numbers.
        vector<double> linespace(double min, double max, size_t N)
        {
            vector<double> linespace;
            double delta = (max-min)/double(N-1);
            for(int i=0; i<N; i++) {
            linespace.push_back(min + i*delta);
                        }
            return linespace;
                        }

    // return Largest Element in Array
    int largest(double arr[], int n)
{
    int i;
    // Initialize maximum element
    int max = arr[0];

    // Traverse array elements from second and compare
    // every element with current max
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];

    return max;
}


int main()
{
 
std::ifstream waveFile;                      //creates stream myFile
waveFile.open("sincos_wave_3Hz_16Hz.txt");  //opens .txt file. Input your file here

std::vector<double>wavedata;  //vector to store the numerical values in
double number = 0 ;
double sum = 0;
while(waveFile >> number)
    {
    wavedata.push_back(number);
    }
    double wavform[wavedata.size()];
    std::copy(wavedata.begin(), wavedata.end(), wavform);// Convert vector to array

    for (int i = 0; i< wavedata.size(); i++){
   // std::cout << wavedata[i] << '\n';
   // cout << wavform[i] << '\n';
    sum += wavedata[i];
    }                        //calculates sum

std::cout << "Average number: " << sum/wavedata.size() << std::endl;  //prints average


// analyze spike of the waveform in a window
 int T1 = 19, T2 = 49; //T = 2900000
 double wavedata_win[0];
 //vector<double>wavedata_win;                //  vector declaration
 int cnt = 0;
 for (int i = T1; i< T2; i++){
      wavedata_win[cnt] = wavform[i];

 cnt += 1;
 }
 cout <<wavedata_win[1]<< '\n';

int sig_len = cnt;
 double real_output [sig_len/2];
 double imx_output [sig_len/2];

    dft_signal_calc(&wavedata_win[0], &real_output[0], &imx_output[0], sig_len);
    for(int i=0;i<sig_len/2;i++)
    {
        cout <<"fft real :" <<real_output[i] <<'\n';


    }
    //*downsample by m
    int Fs = 5* 625;    //176400;          // sampling rate of waveform
    int m = 5;   // 20             // moving average bin;
    int Fsm =(Fs/m);
    int p = div_floor(sig_len,m);       //Window length
    double y_avg[p];
    downsample(&wavedata_win[0],&y_avg[0], p, m);

    for(int i=1;i<=p;i++)
    {
        cout <<"downsample wave :" <<y_avg[i] <<'\n';

    }
    //threshold to find spikes
    int i = 1;
    double spike_hold[0];
     for( int j = 1; j<=p; j++)
            {
                spike_hold[j] = 0;
            }

            spike_func(&y_avg[0],&spike_hold[0], i, p);

    for( int j = 1; j<=p; j++)
            {
                cout <<"spikes :" <<spike_hold[j] <<'\n';
            }

    // finding the frequencies of each spike
    int k = 1, ii_spike = 0;
    double spike[0];
     while (k < p)
     {
         int jj = 1, end_spike = 0;
         double spike_zero[0];
         while (((spike_hold[k] != 0)||(end_spike > 0))&& (k < p))
         {
             spike[jj] = spike_hold[k];
             k += 1;
             jj+= 1;
             if (spike_hold[k] != 0)
             {
                 end_spike = 3;
             }
             else
             {
                 end_spike -= 1;
             }
             if (end_spike == 0)
             {
                 if (jj > 5)
                 {
                 for( int j = 1; j<=jj; j++)
                {
                    spike_zero[j] = 0;
                }
                //Merge two arrays
                int spike_holdsize = sizeof(spike)/sizeof(spike[0]);
                double spike_merge [spike_holdsize+jj];
                for (int j = 1; j<=spike_holdsize; j++)
                {
                    spike_merge[j] = spike[j];
                }
                for (int j = spike_holdsize; j<=spike_holdsize+jj ;j++)
                {
                    spike_merge[j] = spike_zero[j];
                }

                int sig_len = spike_holdsize+jj;
                double spike_real_output [sig_len/2];
                double spike_imx_output [sig_len/2];

                dft_signal_calc(&spike_merge[0], &spike_real_output[0], &spike_imx_output[0], sig_len);

                int spike_size = sizeof(spike)/sizeof(spike[0]);
                int m1 = 0;
                double m2 = 1+spike_size/2;
                vector<double>  lin_vec = linespace(m1, m2, spike_size);
                double lin_vect[lin_vec.size()];
                std::copy(lin_vec.begin(), lin_vec.end(), lin_vect);// Convert vector to array


                double spike_fun[0];
        for (int j = 1; j<=lin_vec.size() ; j++)
        {
            spike_real_output[j] = fabs(spike_real_output[j]);
          //  cout <<"yval hold :" <<abs_val[j]<<'\n';
        }
            int lent = sizeof(lin_vect)/sizeof(lin_vect[0]);
            int index = largest(spike_real_output,lent);
            if (lin_vect[index]>20)
            {
               ii_spike += 1;
               spike_fun[ii_spike] = lin_vect[index]*2;

            }
                k += 50;
                 }
                 delete spike;
                 break;
             }
         }
         k += 1;
     }



    return 0;
}
