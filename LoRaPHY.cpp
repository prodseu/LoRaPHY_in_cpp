#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <complex>
#include <tuple>
#include <cmath>
#define M_PI 3.14159265358979323846

using namespace std;


// 递归 Cooley-Tukey FFT
std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& input, int N) {
    // 填充或截断输入长度
    std::vector<std::complex<double>> data(N, 0);
    int input_size = input.size();
    for (int i = 0; i < min(N, input_size); ++i)
        data[i] = input[i];

    if (N <= 1)
        return data;

    // 奇偶分离
    std::vector<std::complex<double>> even, odd;
    for (int i = 0; i < N / 2; ++i) {
        even.push_back(data[i * 2]);
        odd.push_back(data[i * 2 + 1]);
    }

    // 递归调用
    auto fft_even = fft(even, N / 2);
    auto fft_odd = fft(odd, N / 2);

    std::vector<std::complex<double>> output(N);
    for (int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / N) * fft_odd[k];
        output[k] = fft_even[k] + t;
        output[k + N / 2] = fft_even[k] - t;
    }

    return output;
}

tuple<double, int> topn(vector<double> data)
{
    double mm = 0;
    tuple<double, int> result;
    for(int i = 0; i < data.size(); i++)
    {
        if(data[i] > mm) mm = data[i];
        result = {data[i], i};
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
    vector<double> t;
    for(int i = 0; i < samp_per_sym; i++) t.push_back(i*(N-h)/N/fs*tscale + tdelta);
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

std::vector<double> resample(const std::vector<double>& x, int p, int q) {
    //  上采样：插0
    std::vector<double> upsampled;
    upsampled.reserve(x.size() * p);
    for (auto val : x) {
        upsampled.push_back(val);
        for (int i = 1; i < p; ++i) upsampled.push_back(0.0);
    }

    //  下采样
    std::vector<double> y;
    for (size_t i = 0; i < upsampled.size(); i += q) {
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
            zero_padding_ratio = 10;
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

            // header_checksum_matrix = gf([
            //     1 1 1 1 0 0 0 0 0 0 0 0
            //     1 0 0 0 1 1 1 0 0 0 0 1
            //     0 1 0 0 1 0 0 1 1 0 1 0
            //     0 0 1 0 0 1 0 1 0 1 1 1
            //     0 0 0 1 0 0 1 0 1 1 1 1
            // ]);

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

tuple<double, int> dechirp(int x, bool is_up, Param p, vector<complex<double>> signal);


int detect(int start_idx, vector<complex<double>> signal, Param p)
{
    int ii = start_idx;
    vector<int> pk_bin_list;
    int x;
    while(ii <= signal.size() - p.sample_num * p.preamble_len)
    {
        if(pk_bin_list.size() == p.preamble_len - 1)
        {
            x = ii - round((pk_bin_list.back() - 1) / p.zero_padding_ratio * 2);
            return x;
        }
        tuple<double, int> pk0 = dechirp(ii, 1, p, signal);
        if(!(pk_bin_list.empty()))
        {
            int bin_diff = (pk_bin_list.back() - get<1>(pk0)) % p.bin_num;
            if(bin_diff > p.bin_num / 2) bin_diff = p.bin_num - bin_diff;
            
            
            if(bin_diff <= p.zero_padding_ratio) pk_bin_list.push_back(get<1>(pk0));
            else
            {
                pk_bin_list.clear();
                pk_bin_list.push_back(get<1>(pk0));
            }
        }
        else 
        {
            pk_bin_list.clear();
            pk_bin_list.push_back(get<1>(pk0));
        }

        ii += p.sample_num;
        cout<<pk_bin_list.size();
    }
    x = -1;
    for(int pp = 0; pp < pk_bin_list.size();pp++) cout<<pk_bin_list[pp]<<endl;
    return x;
}





// };

// int sync(int start)
// {
//     bool found = false;
//     while(x < signal.size() - p.sample_num)
//     {
//         up_peak = dechirp(x);
//         down_peak = dechirp(x, false);
//         if(abs(doun_peak))


// };




// bool parse_header(const vector<int>& symbols)
// {

// };
// int calc_sym_num(int payload_len)
// {

// };
// vector<int> dynamic_compensation(const vector<int>& symbols)
// {

// };
tuple<double, int> dechirp(int x, bool is_up, Param p, vector<complex<double>> signal)
{
    vector<complex<double>> c;
    if(is_up) c = p.upchirp;
    else c = p.downchirp;
    vector<complex<double>> subvec(signal.begin() + x, signal.begin() + x + p.sample_num - 1);
    vector<complex<double>> mult;
    for(int i = 0; i < c.size(); i++) mult.push_back(subvec[i] * c[i]);
    vector<complex<double>> ft = fft(mult, p.fft_len);
    vector<complex<double>> slice1(ft.begin(), ft.begin()+p.bin_num);
    vector<complex<double>> slice2(ft.begin()+p.fft_len - p.bin_num+1, ft.begin()+p.fft_len);
    vector<double> ft_1;
    for(int i = 0; i<p.bin_num; i++) ft_1.push_back(abs(slice1[i])+abs(slice2[i]));
    for(int i = 0; i<p.bin_num; i++) ft_1.push_back(i);
    tuple<double, int>pk = topn(ft_1);
    return pk;
};


vector<int32_t> read_file(string filename)
{
    vector<int32_t> result;
    
    ifstream infile(filename, ios::binary);
    if (!infile.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return result;
    }

    int32_t buffer;
    int count = 0;
    while (infile.read(reinterpret_cast<char*>(&buffer), sizeof(buffer))) {
        if(count >= 8) result.push_back((buffer));
        count++;
    }
    infile.close();
    return result;

};

// vector<int> demodulate(vector<complex<double>> signal, Param pr)
// {
//     vector<int> data;
//     pr.cfo = 0;
//     vector<int> symbols_m;
//     vector<double> cfo_m;
//     vector<double> netid_m;
//     int x = 1;
//     while(x < signal.size())
//     {
//         x = detect(x);
//         if(x < 0) break;

//         x = sync(x);

//         // double pk_netid1 = dechirp(round(x - 4.25*sample_num));
//         // double pk_netid2 = dechirp(round(x - 4.25*sample_num));
//         // netid_m.push_back([mod((pk_netid1(2)+ bin_num- preamble_bin)/ zero_padding_ratio, 2^ sf), mod((pk_netid2(2)+ bin_num- preamble_bin)/ zero_padding_ratio, 2^ sf)])
//         vector<int> symbols;
//         vector<tuple<double, int>> pk_list;
//         if (x > signal.size() - 8*pr.sample_num + 1) return vector<int>(-1);
//         for(int ii = 0; ii < 8; ii++)
//         {
//             tuple<double, int> pk = dechirp(x+ii*pr.sample_num);
//             pk_list.push_back(pk);
//             symbols.push_back(((get<1>(pk)+pr.bin_num-pr.preamble_bin)/pr.zero_padding_ratio) % 2^pr.sf);
//         };

//         if(pr.has_header)
//         {
//             bool is_valid = parse_header(symbols);
//             if(~is_valid)
//             {
//                 x = x + 7*pr.sample_num;
//                 continue;
//             };
//         };
//         int sym_num = calc_sym_num(pr.payload_len);
//         if (x > signal.size() - sym_num * pr.sample_num + 1) return vector<int>(-1);
//         for(int ii = 8; ii<sym_num; ii++)
//         {
//             tuple<double, int> pk = dechirp(x+ii*pr.sample_num);
//             pk_list.push_back(pk);
//             symbols.push_back(((get<1>(pk)+pr.bin_num-pr.preamble_bin)/pr.zero_padding_ratio) % 2^pr.sf);
//         };
//         x += sym_num * pr.sample_num;
//         symbols = dynamic_compensation(symbols);
//         for(int jj = 0; jj < symbols.size(); jj++)
//         {
//             symbols_m.push_back(round(symbols[jj] % (2^pr.sf)));
//             cfo_m.push_back(pr.cfo);
//         };
//         if(symbols_m.empty()) printf("No preamble detected!");
//     }

//     return data;
// };


void export_to_csv(const vector<complex<double>>& signal, const string& filename) {
    ofstream file(filename);
    for (const auto& val : signal) {
        file << real(val) << "," << imag(val) << "\n";
    }
    file.close();
}

int main()
{
    // vector<int32_t> signal_int32 = read_file("M:\\lora_modulation\\LoRaPHY\\ebyte_1.2025-07-31T08_04_15_136.sdriq");
    // vector<complex<double>> signal;
    // for(int i = 0; i < signal_int32.size(); i+=2)
    // {
    //     signal.push_back(complex<double>(static_cast<double>(signal_int32[i])/65536 / 65536, static_cast<double>(signal_int32[i+1])/65536 / 65536));
    // }
    // cout<<signal.size();
    vector<complex<double>> signal;
    double fc = 850.150e3;
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
    // data = demodulate(signal, p1);
    // vector<complex<double>> upchirp = chirp(true, p1.sf, p1.bw, 2 * p1.bw, 0, p1.cfo, 0, 1);
    // export_to_csv(upchirp, "M:\\lora_modulation\\LoRaPHY\\shit.csv");

    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\newshit.csv");
    signal = p1.upchirp;
    for(int ppp = 0; ppp < 19; ppp++) signal.insert(signal.end(), p1.upchirp.begin(), p1.upchirp.end());
    int detect_res = 98739;
    detect_res = detect(0, signal, p1);
    cout<<detect_res<<endl;
    cout << signal.size()<<endl;
    cout << p1.upchirp.size();
    // tuple<double, int> pk = dechirp(0, 1, p1, signal);
    // cout<<get<0>(pk)<<'\t'<<get<1>(pk);
    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\topn_of_chirp.csv");

    return 0;
};