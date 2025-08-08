#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <complex>
#include <tuple>
#include <cmath>
#include <ctime>
#include "fftw3.h"
// #include "kiss_fft.h"
// #include "fftw3.h"
#include <windows.h>
#define M_PI 3.14159265358979323846

using namespace std;



// 递归 Cooley-Tukey FFT
// std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& input, int N) {
//     if ((N & (N - 1)) != 0) {
//         throw std::invalid_argument("Input size must be a power of 2");
//     }

//     std::vector<std::complex<double>> data(input);

//     // Bit-reversal permutation
//     int j = 0;
//     for (int i = 1; i < N; ++i) {
//         int bit = N >> 1;
//         while (j & bit) {
//             j ^= bit;
//             bit >>= 1;
//         }
//         j ^= bit;

//         if (i < j) {
//             std::swap(data[i], data[j]);
//         }
//     }

//     // Cooley-Tukey FFT
//     for (int len = 2; len <= N; len <<= 1) {
//         double angle = -2 * M_PI / len;
//         std::complex<double> wlen(cos(angle), sin(angle));
//         for (int i = 0; i < N; i += len) {
//             std::complex<double> w(1);
//             for (int j = 0; j < len / 2; ++j) {
//                 std::complex<double> u = data[i + j];
//                 std::complex<double> t = w * data[i + j + len / 2];
//                 data[i + j] = u + t;
//                 data[i + j + len / 2] = u - t;
//                 w *= wlen;
//             }
//         }
//     }

//     return data;
// }

// void fft(vector<complex<double>> &x, int N) {
//     if (N <= 1) return; // 递归终止条件

//     if(N > x.size())
//     {
//         for(int i = x.size(); i < N; i++) x.push_back(0.0 + 0.0 * 1i); // padding zero
//     }
//     // 拆分为偶数和奇数部分
//     vector<complex<double>> even(N / 2), odd(N / 2);
//     for (int i = 0; i < N / 2; i++) {
//         even[i] = x[i * 2];
//         odd[i] = x[i * 2 + 1];
//     }
 
//     // 递归计算 FFT
//     fft(even, N/2);
//     fft(odd, N/2);
 
//     // 组合结果
//     for (int k = 0; k < N / 2; k++) {
//         complex<double> t = polar(1.0, -2 * M_PI * k / N) * odd[k]; // 旋转因子 e^(-j2πk/N)
//         x[k] = even[k] + t;
//         x[k + N / 2] = even[k] - t;
//     }
// }



tuple<double, double> topn(vector<double> &data)
{
    double mm = 0;
    tuple<double, double> result;
    for(int i = 0; i < data.size(); i++)
    {
        if(data[i] > mm) 
        {
            mm = data[i];
            result = {data[i], i};
        }
    }
    return result;
}

vector<complex<double>> chirp(bool is_up,int sf,double bw,double fs,double h,double cfo,double tdelta,double tscale)
{
    // % chirp  Generate a LoRa chirp symbol
    // %
    // % input:
    // %     is_up: `true` if constructing an up-chirp
    // %            `false` if constructing a down-chirp
    // %     sf: Spreading Factor
    // %     bw: Bandwidth
    // %     fs: Sampling Frequency
    // %     h: Start frequency offset (0 to pow(2, sf)-1)
    // %     cfo: Carrier Frequency Offset
    // %     tdelta: Time offset (0 to 1/fs)
    // %     tscale: Scaling the sampling frequency
    // % output:
    // %     y: Generated LoRa symbol
    double N = pow(2, sf);
    double T = N/bw;
    int samp_per_sym = round(fs/bw*N);
    double h_orig = h;
    h = round(h);
    cfo = cfo + (h_orig - h) / N * bw;
    double k;
    double f0;
    double phi;
    if(is_up)
    {
        k = bw/T;
        f0 = -bw/2+cfo;
    }
    else
    {
        k = -bw/T;
        f0 = bw/2+cfo;
    }

    // % retain last element to calculate phase
    vector<double> t(samp_per_sym);
    for(int i = 0; i < samp_per_sym; i++) t[i] = (i*(N-h)/N/fs*tscale + tdelta);
    int snum = t.size();
    vector<complex<double>> c1;
    vector<complex<double>> c2;
    vector<complex<double>> y;
    for(int i = 0; i < round(samp_per_sym*(N-h)/N); i++) 
    {
        double tt = i/fs*tscale + tdelta;
        double temp = double(2*M_PI*(tt*(f0+k*T*h/N+0.5*k*tt)));
        c1.push_back(exp(1i*temp));
    }

    if(snum == 0) phi = 0;
    else phi = atan(c1[snum].imag() / c1[snum].real());

    for(int i = 0; i < round(samp_per_sym*h/N-1); i++) 
    {
        double tt = i/fs + tdelta;
        double temp = double((phi + 2*M_PI*(tt*(f0+0.5*k*tt))));
        c2.push_back(exp(1i*temp));
    }
    for(int i = 0; i < snum; i++)
    {
        y.push_back(c1[i]);
    }
    for(int i = 0; i < c2.size(); i++)
    {
        y.push_back(c2[i]);
    }
    return y;
};

std::vector<complex<double>> resample(const std::vector<complex<double>>& x, int p, int q) {
    //  上采样：插0
    std::vector<complex<double>> upsampled;
    upsampled.reserve(x.size() * p);
    for (auto val : x) {
        upsampled.push_back(val);
        for (int i = 1; i < p; ++i) upsampled.push_back(0.0 + 1i * 0.0);
    }

    //  下采样
    std::vector<complex<double>> y;
    y.reserve(round(upsampled.size() / q));
    for (int i = 0; i < upsampled.size(); i += q) {
        y.push_back(upsampled[i]);
    }

    return y;
}

class Param
{   
    public:
        double rf_freq;                   //% carrier frequency
        int sf;                        //% spreading factor (7,8,9,10,11,12)
        double bw;                        //% bandwidth (125kHz 250kHz 500kHz)
        double fs;                        //% sampling frequency
        int cr;                        //% code rate: (1:4/5 2:4/6 3:4/7 4:4/8)
        int payload_len;               //% payload length
        int has_header;                //% explicit header: 1, implicit header: 0
        int crc;                       //% crc = 1 if CRC Check is enabled else 0
        bool ldr;                       //% ldr = 1 if Low Data Rate Optimization is enabled else 0
        vector<uint8_t> whitening_seq;             //% whitening sequence
        // crc_generator;             //% CRC generator with polynomial x^16+x^12+x^5+1
        vector<int> header_checksum_matrix;    //% we use a 12 x 5 matrix to calculate header checksum
        int preamble_len;              //% preamble length
        vector<complex<double>> downchirp;                 //% ideal chirp with decreasing frequency from B/2 to -B/2
        vector<complex<double>> upchirp;                   //% ideal chirp with increasing frequency from -B/2 to B/2
        int sample_num;                //% number of sample points per symbol
        int bin_num;                   //% number of bins after FFT (with zero padding)
        int zero_padding_ratio;        //% FFT zero padding ratio
        int fft_len;                   //% FFT size
        int preamble_bin;              //% reference bin in current decoding window, used to eliminate CFO
        double cfo;                       //% carrier frequency offset
        bool hamming_decoding_en;       // enable hamming decoding

        Param(double rf_freq_i,int sf_i,double bw_i,double fs_i, bool ldr_i)
        {
            rf_freq = rf_freq_i;
            sf = sf_i;
            bw = bw_i;
            fs = fs_i;
            has_header = 1;
            crc = 1;
            hamming_decoding_en = true;
            zero_padding_ratio = 4;
            cfo = 0;

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // % The whitening sequence is generated by an LFSR
            // % x^8+x^6+x^5+x^4+1
            // % Use the code below to generate such sequence
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // % reg = 0xFF;
            // % for i = 1:255
            // %     fprintf("0x%x, ", reg);
            // %     reg = bitxor(bitshift(reg,1), bitxor(bitget(reg,8), bitxor(bitget(reg,6), bitxor(bitget(reg,5), bitget(reg,4)))));
            // % end
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            whitening_seq = {0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0xb, 0x17, 0x2f, 0x5e, 0xbc, 0x78, 0xf1, 0xe3, 0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40, 0x80, 0x1, 0x2, 0x4, 0x8, 0x11, 0x23, 0x47, 0x8e, 0x1c, 0x38, 0x71, 0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0, 0xc0, 0x81, 0x3, 0x6, 0xc, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4d, 0x9b, 0x37, 0x6e, 0xdc, 0xb9, 0x72, 0xe4, 0xc8, 0x90, 0x20, 0x41, 0x82, 0x5, 0xa, 0x15, 0x2b, 0x56, 0xad, 0x5b, 0xb6, 0x6d, 0xda, 0xb5, 0x6b, 0xd6, 0xac, 0x59, 0xb2, 0x65, 0xcb, 0x96, 0x2c, 0x58, 0xb0, 0x61, 0xc3, 0x87, 0xf, 0x1f, 0x3e, 0x7d, 0xfb, 0xf6, 0xed, 0xdb, 0xb7, 0x6f, 0xde, 0xbd, 0x7a, 0xf5, 0xeb, 0xd7, 0xae, 0x5d, 0xba, 0x74, 0xe8, 0xd1, 0xa2, 0x44, 0x88, 0x10, 0x21, 0x43, 0x86, 0xd, 0x1b, 0x36, 0x6c, 0xd8, 0xb1, 0x63, 0xc7, 0x8f, 0x1e, 0x3c, 0x79, 0xf3, 0xe7, 0xce, 0x9c, 0x39, 0x73, 0xe6, 0xcc, 0x98, 0x31, 0x62, 0xc5, 0x8b, 0x16, 0x2d, 0x5a, 0xb4, 0x69, 0xd2, 0xa4, 0x48, 0x91, 0x22, 0x45, 0x8a, 0x14, 0x29, 0x52, 0xa5, 0x4a, 0x95, 0x2a, 0x54, 0xa9, 0x53, 0xa7, 0x4e, 0x9d, 0x3b, 0x77, 0xee, 0xdd, 0xbb, 0x76, 0xec, 0xd9, 0xb3, 0x67, 0xcf, 0x9e, 0x3d, 0x7b, 0xf7, 0xef, 0xdf, 0xbf, 0x7e, 0xfd, 0xfa, 0xf4, 0xe9, 0xd3, 0xa6, 0x4c, 0x99, 0x33, 0x66, 0xcd, 0x9a, 0x35, 0x6a, 0xd4, 0xa8, 0x51, 0xa3, 0x46, 0x8c, 0x18, 0x30, 0x60, 0xc1, 0x83, 0x7, 0xe, 0x1d, 0x3a, 0x75, 0xea, 0xd5, 0xaa, 0x55, 0xab, 0x57, 0xaf, 0x5f, 0xbe, 0x7c, 0xf9, 0xf2, 0xe5, 0xca, 0x94, 0x28, 0x50, 0xa1, 0x42, 0x84, 0x9, 0x13, 0x27, 0x4f, 0x9f, 0x3f, 0x7f};

            // header_checksum_matrix = 
            //     {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            //     {1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
            //     {0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0},
            //     {0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1},
            //     {0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1}};
            // 00001110001
            header_checksum_matrix = {3840, 2273, 1178, 599, 303};
            

            // crc_generator = comm.CRCGenerator('Polynomial','X^16 + X^12 + X^5 + 1');

            preamble_len = 6;
            bin_num = pow(2, sf) * zero_padding_ratio;
            sample_num = 2 * pow(2, sf);
            fft_len =  sample_num * zero_padding_ratio;

            downchirp = chirp(false, sf, bw, 2 * bw, 0, cfo, 0, 1);
            upchirp = chirp(true, sf, bw, 2 * bw, 0, cfo, 0, 1);
            
            ldr = ldr_i;
        };
        ~Param()
        {

        };
        
};

int sync(int x, vector<complex<double>> &signal, Param &p);
int detect(int start_idx, vector<complex<double>> &signal, Param &p);
tuple<double, double> dechirp(int x, bool is_up, Param &p, vector<complex<double>> &signal);
bool parse_header(const vector<int>& symbols);
double calc_sym_num(Param &p);
vector<double> dynamic_compensation(const vector<double>& data, Param &p);
vector<int> demodulate(vector<complex<double>> &signal, Param &p);
vector<uint16_t> gray_coding(vector<int> &din, Param &p);
vector<uint16_t> diag_deinterleave(vector<uint16_t> symbols_g, int ppm);
uint16_t parity_fix(uint16_t p);
vector<uint16_t> hamming_decode(vector<uint16_t> code, int rdd);
vector<uint16_t> dewhiten(vector<uint16_t> &data, Param &p);
vector<uint16_t> decode(vector<int> &symbols_m, Param &p);
vector<complex<double>> read_file(const string& filename);


tuple<double, double> dechirp(int x, bool is_up, Param &p, vector<complex<double>> &signal)
{
    // double start_time = clock();
    fftw_plan p_fft;
    fftw_complex *din, *dout;
    din = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * p.fft_len);
	dout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * p.fft_len);
    if(!is_up) 
    {
        for(int i = 0; i < p.upchirp.size(); i++)
        {
            din[i][0] = (signal[x + i] * p.upchirp[i]).real();
            din[i][1] = (signal[x + i] * p.upchirp[i]).imag();
        }
    }
    else
    {
        for(int i = 0; i < p.upchirp.size(); i++)
        {
            din[i][0] = (signal[x + i] * p.downchirp[i]).real();
            din[i][1] = (signal[x + i] * p.downchirp[i]).imag();
        }
    }

    for(int i = p.upchirp.size(); i < p.fft_len; i++)
    {
        din[i][0] = 0;
        din[i][1] = 0;
    }
    

    // fft(c, p.fft_len);
    p_fft = fftw_plan_dft_1d(p.fft_len, din, dout, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_fft); /* repeat as needed */
	fftw_destroy_plan(p_fft);
	fftw_cleanup();

    vector<double> ft_1(p.bin_num);
    for (int i = 0; i < p.bin_num; i++) {
        // ft_1.push_back(sqrt(out[i].i * out[i].i + out[i].r * out[i].r) + sqrt(out[i+ p.bin_num].i * out[i+ p.bin_num].i + out[i+ p.bin_num].r * out[i+ p.bin_num].r));
        ft_1[i] = (sqrt(dout[i][0]*dout[i][0] + dout[i][1]*dout[i][1]) + sqrt(dout[i+p.bin_num][0]*dout[i+p.bin_num][0] + dout[i+p.bin_num][1]*dout[i+p.bin_num][1]));
    }
    tuple<double, double>pk = topn(ft_1);
    // double end_time = clock();
    // cout<<"time:"<<(end_time - start_time) /  CLOCKS_PER_SEC<<endl;
    return pk;
    
};


int detect(int start_idx, vector<complex<double>> &signal, Param &p)
{
    int ii = start_idx;
    vector<int> pk_bin_list;
    int x = 0;
    while(ii <= signal.size() - p.sample_num * p.preamble_len)
    {
        if(pk_bin_list.size() == p.preamble_len - 1)
        {
            x = ii - round((pk_bin_list.back()) / p.zero_padding_ratio * 2);
            return x;
        }
        tuple<double, double> pk0 = dechirp(ii, 1, p, signal);
        // cout << get<0>(pk0) << '\t' << get<1>(pk0)<<endl;
        if(!(pk_bin_list.empty()))
        {
            int bin_diff = fmod(abs(pk_bin_list.back() - get<1>(pk0)), p.bin_num);
            
            if(bin_diff > p.bin_num / 2) 
            {
                bin_diff = p.bin_num - bin_diff;
            }
            // cout<<bin_diff<<endl;
            if(bin_diff <= p.zero_padding_ratio) 
            {
                pk_bin_list.push_back(get<1>(pk0));
            }
            else
            {
                pk_bin_list.clear();
                pk_bin_list.push_back(get<1>(pk0));
            }
        }
        else 
        {
            pk_bin_list.push_back(get<1>(pk0));
        }

        ii += p.sample_num;
        // cout<<pk_bin_list.back()<<endl;
    }
    x = -1;
    return x;
}




int sync(int x, vector<complex<double>> &signal, Param &p)
{
    bool found = false;
    while(x < signal.size() - p.sample_num)
    {
        tuple<double, double> up_peak, down_peak;
        up_peak = dechirp(x, true, p, signal);
        down_peak = dechirp(x, false, p, signal);
        if(abs(get<0>(down_peak)) > abs(get<0>(up_peak)))
        {
            found = true;
        }
        x += p.sample_num;
        if(found) break;
        
    }
    if(!found) return -1;
    tuple<double, double> pkd = dechirp(x, false, p, signal);
    int to = 0;
    if(get<1>(pkd) > p.bin_num / 2)
    {
        to = round((get<1>(pkd) - p.bin_num) / p.zero_padding_ratio);
    }
    else
    {
        to = round((get<1>(pkd)) / p.zero_padding_ratio);
    }
    x += to;
    tuple<double, double> pku = dechirp(x - 4*p.sample_num, true, p, signal);
    p.preamble_bin = get<1>(pku);
    if(p.preamble_bin > p.bin_num / 2)
    {
        p.cfo = double(p.preamble_bin-p.bin_num)*p.bw/double(p.bin_num);
    }
    else
    {
        p.cfo = double(p.preamble_bin)*p.bw/double(p.bin_num);
    }
    pku = dechirp(x-p.sample_num, true, p, signal);
    pkd = dechirp(x-p.sample_num, false, p, signal);
    if(abs(get<0>(pku) > abs(get<0>(pkd))))
    {
        x += round(2.25*p.sample_num);
    }
    else
    {
        x += round(1.25*p.sample_num);
    }

    return x;
};




bool parse_header(vector<int>& symbols, Param &p)
{
    // vector<int> data;
    vector<uint16_t> symbols_g = gray_coding(symbols, p);
    // for(int i = 0 ;i < symbols_g.size(); i++) cout<<symbols_g[i]<<endl;
    // vector<uint16_t> head(symbols_g.begin(), symbols_g.begin() + 8);
    vector<uint16_t> codewords = diag_deinterleave(symbols_g, p.sf - 2);
    vector<uint16_t> nibbles;
    vector<int> header_checksum;
    vector<int> header_checksum_calc(5);
    nibbles = hamming_decode(codewords, 8);
    // for(int opop = 0; opop < nibbles.size(); opop++) cout<<nibbles[opop]<<endl;
    p.payload_len = int(nibbles[0] * 16 + nibbles[1]);
    // cout<<p.payload_len<<endl;
    p.crc = int(nibbles[2] & 1);
    // cout<<p.crc<<endl;

    p.cr = int(nibbles[2] >> 1);
    // cout<<p.cr<<endl;

    header_checksum.push_back(nibbles[3] & 1);
    for(int i = 3; i >= 0; i--) header_checksum.push_back((nibbles[4] >> i) & 1);
    // cout<<(nibbles[1]<<4)<<endl;
    int header_checksum_seq = (nibbles[0] << 8) + (nibbles[1] << 4) + nibbles[2];
    for(int j = 0; j < 5; j++)
    {
        header_checksum_calc[j] = p.header_checksum_matrix[j] & header_checksum_seq;
        int cnt = 0;
        for(int i = 0; i < 12; i++) if((header_checksum_calc[j] >> i) & 1 == 1) cnt += 1;
        header_checksum_calc[j] = cnt % 2;
    }
    // for(int opop = 0; opop < header_checksum_calc.size(); opop++) cout<<header_checksum_calc[opop]<<endl;
    if (header_checksum != header_checksum_calc)
    {
        
        throw("Invalid header checksum!");
        return false;
    }
    return true;

};


double calc_sym_num(Param &p)
{
    return double(8 + max((4+p.cr)*ceil(double((2*p.payload_len - p.sf+7+4*p.crc-5*(1-p.has_header))) / double(p.sf-2*p.ldr)), 0.0));
};


vector<double> dynamic_compensation(const vector<double>& data, Param &p)
{
    vector<double> symbols(data.size());
    for(int i = 0; i < data.size(); i++)
    {
        symbols[i] = (fmod((data[i] - (1 + i) * pow(2, p.sf) * p.cfo / p.rf_freq), (pow(2, p.sf))));
    }

            if(p.ldr)
            {
                double bin_offset = 0;
                double v_last = 1;
                double v;

                for(int i = 0; i < (symbols.size()); i++)
                {
                    v = symbols[i];
                    double bin_delta = fmod((fmod((v-v_last), 4) + 4), 4);
                    if(bin_delta < 2) bin_offset = bin_offset - bin_delta;
                    else bin_offset = bin_offset - bin_delta + 4;
                    v_last = v;
                    symbols[i] = fmod((v+bin_offset), (pow(2, p.sf)));
                }
            }
    return symbols;
}



vector<int> demodulate(vector<complex<double>> &signal, Param &p)
{
    p.cfo = 0;
    vector<int> symbols_m;
    vector<double> cfo_m;
    vector<double> netid_m;
    int x = 0;
    signal = resample(signal, 1, 2);
    while(x < signal.size())
    {
        
        x = detect(x, signal, p);
        if(x < 0) break;
        // cout<<x<<endl;
        x = sync(x, signal, p);
        // cout<<x<<endl;

        // double pk_netid1 = dechirp(round(x - 4.25*sample_num));
        // double pk_netid2 = dechirp(round(x - 4.25*sample_num));
        // netid_m.push_back([mod((pk_netid1(2)+ bin_num- preamble_bin)/ zero_padding_ratio, 2^ sf), mod((pk_netid2(2)+ bin_num- preamble_bin)/ zero_padding_ratio, 2^ sf)])
        vector<double> symbols;
        vector<int> head;
        vector<tuple<double, double>> pk_list;
        
        if (x > (signal.size() - 8*p.sample_num + 1)) return {-1};
        for(int ii = 0; ii < 8; ii++)
        {
            tuple<double, double> pk = dechirp(x+ii*p.sample_num, 1, p, signal);
            pk_list.push_back(pk);
            // cout<<'p' << '\t' << get<1>(pk)<<endl;
            symbols.push_back(((fmod(((get<1>(pk)+p.bin_num-p.preamble_bin)/p.zero_padding_ratio) , (pow(2,p.sf))))));
            head.push_back(int(round(fmod(((get<1>(pk)+p.bin_num-p.preamble_bin)/p.zero_padding_ratio) , (pow(2,p.sf))))));
            // cout<< (get<1>(pk)+p.bin_num-p.preamble_bin)/p.zero_padding_ratio << endl;
        };
        // for(int cnt = 0; cnt < symbols.size(); cnt++) cout<<symbols[cnt]<<endl;

        if(p.has_header)
        {
            bool is_valid = parse_header(head, p);
            if(!is_valid)
            {
                x = x + 7*p.sample_num;
                continue;
            };
        };
        double sym_num = calc_sym_num(p);
        // int sym_num = 28;
        // cout<<"sym_num  "<<sym_num<<endl;

        if (x > (signal.size() - sym_num * p.sample_num + 1)) return {-1};
        for(int ii = 8; ii<sym_num; ii++)
        {
            tuple<double, double> pk = dechirp(x+ii*p.sample_num, 1, p, signal);
            pk_list.push_back(pk);
            symbols.push_back(((fmod(((get<1>(pk)+p.bin_num-p.preamble_bin)/p.zero_padding_ratio) , (pow(2,p.sf))))));
        };

        x += sym_num * p.sample_num;
        // cout<<"cfo "<<p.cfo<<endl;
        symbols = dynamic_compensation(symbols, p);
        // for(int iii = 0; iii < symbols.size(); iii++) cout<<"sym  "<<symbols[iii]<<endl;

        for(int jj = 0; jj < symbols.size(); jj++)
        {
            symbols_m.push_back(round(fmod(symbols[jj] , (pow(2, p.sf)))));
            cfo_m.push_back(p.cfo);
        };
        if(symbols_m.empty()) printf("No preamble detected!");
    }

    return symbols_m;
}

vector<uint16_t> gray_coding(vector<int> &din, Param &p)
{
    vector<uint16_t> s;
    for(int i = 0; i < din.size(); i++)
    {
        s.push_back(uint16_t(floor(din[i] / 4)));
        s[i] = s[i] ^ (s[i] >> 1);
    }
    return s;

}

vector<uint16_t> diag_deinterleave(vector<uint16_t> symbols_g, int ppm)
{
    int N = symbols_g.size();
    // uint16_t order[] = {1,2,3,5,4,0,6,7};
    uint16_t order[] = {0,1,2,3,4,5,6,7}; //并没有交织，只是做了比特翻转

    vector<uint16_t> dout(ppm);
    int b = 0;
    for(int i = 0; i < ppm; i++) dout[i] = 0;
    for(int i = 0; i < ppm; i++)
    {
        for(int x = 0; x < N; x++)
        {
            b = (((-x+i) % ppm + ppm * (int(ppm / N) + 1)) % ppm);
            dout[i] += ((symbols_g[x] >> b) & 1) << order[x];
            // cout<<b;
            // cout<<((symbols_g[x] >> b) & 1);
        }
        // cout<<endl;
    }
    return dout;
}

uint16_t parity_fix(uint16_t p)
{
    switch(p)
    {
        case 3: return 4;
        case 5: return 8;
        case 6: return 1;
        case 7: return 2;
        default: return 0;
    }
}


vector<uint16_t> hamming_decode(vector<uint16_t> code, int rdd)
{
    vector<uint16_t> p1,p2,p3,p4,p5;
    vector<uint16_t> nibbles;
    for(int i = 0; i < code.size(); i++)
    {
        uint16_t b0 = (code[i] >> 0) & 1;
        uint16_t b1 = (code[i] >> 1) & 1;
        uint16_t b2 = (code[i] >> 2) & 1;
        uint16_t b3 = (code[i] >> 3) & 1;
        uint16_t b4 = (code[i] >> 4) & 1;
        uint16_t b5 = (code[i] >> 5) & 1;
        uint16_t b6 = (code[i] >> 6) & 1;
        uint16_t b7 = (code[i] >> 7) & 1;
        p1.push_back(b7 ^ b3 ^ b2 ^ b0);
        p2.push_back(b6 ^ b3 ^ b1 ^ b0);
        p3.push_back(b4 ^ b2 ^ b1 ^ b0);
        p4.push_back(b4 ^ b3 ^ b2 ^ b1 ^ b0);
        p5.push_back(b5 ^ b3 ^ b2 ^ b1);
    }
    switch(rdd)
    {
        case 5:
        case 6:
            for(int i = 0; i < code.size(); i++) nibbles.push_back(code[i] % 16);
            break;
        case 7:
        case 8:
            vector<uint16_t> parity;
            for(int i = 0; i < code.size(); i++) 
            {
                uint16_t parity = (p2[i] * 4 + p3[i] * 2 + p5[i]);
                code[i] = code[i] ^ parity;
                nibbles.push_back(code[i] % 16);
            }
            break;

    }
    return nibbles;
}


vector<uint16_t> dewhiten(vector<uint16_t> &data, Param &p)
{
    int len = data.size();
    for(int i = 0; i < len; i++)
    {
        data[i] = data[i] ^ p.whitening_seq[i];
    }
    return data;
}

// vector<uint16_t> calc_crc(vector<uint16_t> &data, Param &p)
// {
//     vector<uint16_t> checksum(2);
//     uint16_t crcode = 0x1021; // 16 12 5 0
//     switch(data.size())
//     {
//         case 0:
//             checksum = {0, 0};
//             break;
//         case 1:
//             checksum = {data[0], 0};
//             break;
//         case 2:
//             checksum = {data[0], data[1]};
//             break;
//         defalut:
//             uint16_t b1, b2;
//             vector<uint16_t> sub_data(data.begin(), data.end() - 2);
//             b1 = 
//     }


// }


vector<uint16_t> decode(vector<int> &symbols_m, Param &p)
{
    vector<uint16_t> data_m;
    vector<int> data;
    vector<uint16_t> symbols_g = gray_coding(symbols_m, p);
    // for(int i = 0 ;i < symbols_g.size(); i++) cout<<symbols_g[i]<<endl;
    vector<uint16_t> head(symbols_g.begin(), symbols_g.begin() + 8);
    vector<uint16_t> codewords = diag_deinterleave(head, p.sf - 2);
    vector<uint16_t> nibbles;
        
    if(p.has_header)
    {
        vector<int> header_checksum;
        vector<int> header_checksum_calc(5);
        nibbles = hamming_decode(codewords, 8);
        // for(int opop = 0; opop < nibbles.size(); opop++) cout<<nibbles[opop]<<endl;
        p.payload_len = int(nibbles[0] * 16 + nibbles[1]);
        // cout<<p.payload_len<<endl;
        p.crc = int(nibbles[2] & 1);
        // cout<<p.crc<<endl;

        p.cr = int(nibbles[2] >> 1);
        // cout<<p.cr<<endl;

        header_checksum.push_back(nibbles[3] & 1);
        for(int i = 3; i >= 0; i--) header_checksum.push_back((nibbles[4] >> i) & 1);
        cout<<(nibbles[1]<<4)<<endl;
        int header_checksum_seq = (nibbles[0] << 8) + (nibbles[1] << 4) + nibbles[2];
        for(int j = 0; j < 5; j++)
        {
            header_checksum_calc[j] = p.header_checksum_matrix[j] & header_checksum_seq;
            int cnt = 0;
            for(int i = 0; i < 12; i++) if((header_checksum_calc[j] >> i) & 1 == 1) cnt += 1;
            header_checksum_calc[j] = cnt % 2;
        }
        // for(int opop = 0; opop < header_checksum_calc.size(); opop++) cout<<header_checksum_calc[opop]<<endl;
        if (header_checksum != header_checksum_calc) throw("Invalid header checksum!");
        nibbles.erase(nibbles.begin(), nibbles.begin() + 5); 


    }

    int rdd = p.cr + 4;
    for(int ii = 8; ii < symbols_g.size() - rdd + 1; ii += rdd)
    {
        vector<uint16_t> sub_sym_g(symbols_g.begin() + ii, symbols_g.begin() + ii + rdd);
        codewords = diag_deinterleave(sub_sym_g, p.sf - 2 * p.ldr);
        vector<uint16_t> new_nibbles = hamming_decode(codewords, rdd);
        nibbles.insert(nibbles.end(), new_nibbles.begin(), new_nibbles.end());

    }
    int len_bytes = min(255, int(floor(nibbles.size() / 2)));
    vector<uint16_t> bytes(len_bytes);
    for(int ii = 0; ii < len_bytes; ii++)
    {
        bytes[ii] = uint16_t(nibbles[2 * ii]) | (16*uint16_t(nibbles[2*ii + 1]));
    }

    vector<uint16_t> sub_bytes(bytes.begin(), bytes.begin() + p.payload_len - 1);
    // for(int i = 0; i < sub_bytes.size(); i++) cout<< "byte "<<sub_bytes[i]<<endl;

    if(p.crc)
    {
        data_m = dewhiten(sub_bytes, p);
        // vector<uint16_t> checksum = calc_crc(data_m);
        data_m.push_back(bytes[p.payload_len]);
        data_m.push_back(bytes[p.payload_len + 1]);
        
    }
    else
    {
        data_m = dewhiten(sub_bytes, p);
    }
    data_m.insert(data_m.end(), data.begin(), data.end());

return data_m;
}




// vector<int32_t> read_file(string filename)
// {
//     vector<int32_t> result;
    
//     ifstream infile(filename, ios::binary);
//     if (!infile.is_open()) {
//         cerr << "Cannot open file: " << filename << endl;
//         return result;
//     }

//     int32_t buffer;
//     int count = 0;
//     while (infile.read(reinterpret_cast<char*>(&buffer), sizeof(buffer))) {
//         if(count >= 8) result.push_back((buffer));
//         count++;
//     }
//     infile.close();
//     return result;

// };


vector<complex<double>> read_file(const string& filename) {
    vector<complex<double>> result;

    ifstream infile(filename, ios::binary);
    if (!infile.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return result;
    }

    int32_t real_part, imag_part;
    int count = 0;

    while (infile.read(reinterpret_cast<char*>(&real_part), sizeof(int32_t))) {
        if (!infile.read(reinterpret_cast<char*>(&imag_part), sizeof(int32_t))) {
            cerr << "Incomplete complex number at the end of file." << endl;
            break;
        }

        // 跳过前 4 个复数（等价于前 8 个 int32_t）
        if (count >= 4) {
            result.emplace_back(static_cast<double>(real_part), static_cast<double>(imag_part));
        }

        count++;
    }

    infile.close();
    return result;
}


wstring decode_gbk_from_uint16_bytes(const std::vector<uint16_t>& input) {
    // 构造一个 std::string，只取每个 uint16_t 的低8位
    std::string gbk_str;
    for (uint16_t val : input) {
        gbk_str.push_back(static_cast<char>(val & 0xFF)); // 只保留低8位作为有效字节
    }

    // 转换 GBK 到 Unicode（UTF-16 / std::wstring）
    int len = MultiByteToWideChar(936, 0, gbk_str.c_str(), static_cast<int>(gbk_str.size()), NULL, 0);
    std::wstring result(len, L'\0');
    MultiByteToWideChar(936, 0, gbk_str.c_str(), static_cast<int>(gbk_str.size()), &result[0], len);

    return result;
}

void export_to_csv(const vector<complex<double>>& signal, const string& filename)
{
    ofstream file(filename);
    for (const auto& val : signal) {
        file << real(val) << "," << imag(val) << "\n";
    }
    file.close();
}

int main()
{
    double start_time = clock();
    vector<complex<double>> signal = read_file("M:\\lora_modulation\\LoRaPHY\\iqdata\\ebyte_1.2025-07-31T08_04_15_136.sdriq");
    double fc = 850.150e6;
    double bw = 500e3;
    int sf = 11;
    double fs = 2e6;
    int has_header = 1;        // explicit header mode
    int cr = 1;                // code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
    int crc = 1;               // enable payload CRC checksum
    int preamble_len = 6;       // preamble: 6 basic upchirps
    bool ldr = 1;                // 强制使用ldr，具体在loraphy.m的init()
    Param p1(fc,sf,bw,fs,ldr);
    vector<int> data;
    
    // signal = resample(signal, 1, 2);
    // cout<<signal.size()<<endl;
    data = demodulate(signal, p1);
    
    for(int i = 0; i < data.size(); i++) cout<<data[i]<<endl;
    vector<uint16_t> decode_res;
    vector<char> decode_res_gbk;
    decode_res = decode(data, p1);
    double end_time = clock();
    cout<<"time:"<<(end_time - start_time) /  CLOCKS_PER_SEC<<endl;
    cout<<"decode:"<<endl;
    for(int i = 0; i < decode_res.size(); i++) 
    {
        decode_res[i] ^= 0x12;

        cout<<decode_res[i]<<endl;
    }
    std::wstring decoded = decode_gbk_from_uint16_bytes(decode_res);

    std::wcout << L"解码结果: " << decoded << std::endl;
    // int outp = detect(0, signal, p1);
    // int synp = sync(outp, signal, p1);
    // cout<<outp<<endl;
    // cout<<synp<<endl;
    // vector<complex<double>> upchirp = chirp(true, p1.sf, p1.bw, 2 * p1.bw, 0, p1.cfo, 0, 1);
    // export_to_csv(upchirp, "M:\\lora_modulation\\LoRaPHY\\shit.csv");

    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\newshit.csv");
    // tuple<double, double> pk = dechirp(0, 1, p1, signal);
    // cout<<get<0>(pk)<<'\t'<<get<1>(pk);
    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\topn_of_chirp.csv");

    return 0;
};