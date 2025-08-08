#include "loracpp.h"

int main()
{
    double start_time = clock();
    vector<complex<double>> signal = read_file("M:\\lora_modulation\\LoRaPHY\\test_0731_hanzi.2025-07-31T01_35_49_110.sdriq");
    // vector<complex<double>> signal = read_file("M:\\lora_modulation\\LoRaPHY\\0808-1.2025-08-08T00_42_18_552.sdriq");
    // vector<complex<double>> signal = read_file("M:\\lora_modulation\\LoRaPHY\\0808-1.2025-08-08T01_27_31_935.sdriq");
    complex<double> sum = std::accumulate(signal.begin(), signal.end(), complex<double>(0.0, 0.0));
    complex<double> mean = sum / static_cast<double>(signal.size());
    for(int i = 0; i < signal.size(); i++) signal[i] -= mean;

    double fc = 800e6;
    double bw = 250e3;
    int sf = 10;
    double fs = 2.4e6;
    int has_header = 1;        // explicit header mode
    int cr = 1;                // code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
    int crc = 1;               // enable payload CRC checksum
    int preamble_len = 6;       // preamble: 6 basic upchirps
    bool ldr = 1;                // 强制使用ldr，具体在loraphy.m的init()
    
    vector<int> data;
    vector<uint16_t> decode_res;
    // int sp_gcd = gcd(int(fs), int(2*bw));

    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\shit_signal.csv");

    sf = estimate_sf(signal,128,64,128,5,fs,bw);
    bw = bw_round(bw);
    cout<<"Estimate sf:"<<sf<<endl;
    cout<<"new bw: "<<bw<<endl;
    // bw = 203.125e3;
    // sf = 10;


    signal = resample(signal, 2 * bw / fs);
    Param p1(fc,sf,bw,fs,ldr);
    data = demodulate(signal, p1);
    if(data.empty()) goto mark;
    decode_res = decode(data, p1);
    if(check_data_crc(decode_res))
    {
        cout<<"Data CRC Valid! Decode:"<<endl;
        for(int i = 0; i < decode_res.size(); i++) 
        {
            // decode_res[i] ^= 0x12;
            cout<<decode_res[i]<<endl;
        }
        string decoded = decode_gbk_from_uint16_bytes(decode_res);
        cout << decoded <<endl;
    }
    else
    {
        cout<<"Data CRC Invalid!"<<endl;
    }
        
    // int i = 7;
    // for(i = 7; i < 13; i++)
    // {
    //     Param p1(fc,i,bw,fs,ldr);
    //     data.clear();
    //     data = demodulate(signal, p1);
    //     cout<<"try "<<i<<endl;
    //     // for(int pp = 0; pp < data.size(); pp++) cout<<data[pp]<<endl;
    //     if(data.empty() ||  data[0] == -1) continue;
    //     decode_res = decode(data, p1);

    //     cout<<"decode:"<<endl;
    //     for(int i = 0; i < decode_res.size(); i++) 
    //     {
    //         decode_res[i] ^= 0x12;

    //         cout<<decode_res[i]<<endl;
    //     }
    //     string decoded = decode_gbk_from_uint16_bytes(decode_res);
    //     cout << decoded <<endl;
    //     break;
    // }
mark:
    double end_time = clock();
    cout<<"time:"<<(end_time - start_time) /  CLOCKS_PER_SEC<<endl;
    // cout<<"你好"<<endl;
    // wofstream fout("M:\\lora_modulation\\LoRaPHY\\decoded_output.txt");
    // // fout.imbue(locale(locale(), new codecvt_utf8<wchar_t>));
    // fout << decoded;
    // fout.close();
    // int outp = detect(0, signal, p1);
    // int synp = sync(outp, signal, p1);
    // cout<<outp<<endl;
    // cout<<synp<<endl;
    // vector<complex<double>> upchirp = chirp(true, p1.sf, p1.bw, 2 * p1.bw, 0, p1.cfo, 0, 1);
    // export_to_csv(upchirp, "M:\\lora_modulation\\LoRaPHY\\shit.csv");

    // tuple<double, double> pk = dechirp(0, 1, p1, signal);
    // cout<<get<0>(pk)<<'\t'<<get<1>(pk);
    // export_to_csv(signal, "M:\\lora_modulation\\LoRaPHY\\topn_of_chirp.csv");

    return 0;
};