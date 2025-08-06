function [decode_data] = demod_lora(signal, fc, fs, bw, sf, preamble_len, has_header)
    phy = LoRaPHY(fc, sf, bw, fs);
    phy.preamble_len = preamble_len;
    phy.has_header = has_header;
    [symbols_d, ~, ~] = phy.demodulate(signal);
    [data, checksum] = phy.decode(symbols_d);
    if ~isempty(checksum)
        if (checksum == data(end-1:end))
            fprintf("Valid!");
            key = uint8(18);  % 0x12
            decode_data = bitxor(uint8(data), key);
%             str = native2unicode(uint8(decode_data(1:end)'), 'GBK');
%             decode_ascii = char(decode_data(1:end))';
        else
            fprintf("Invalid!!")
            decode_data = -1;
        end
    end











end