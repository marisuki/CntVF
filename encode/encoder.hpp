#ifndef ENCODER_H
#define ENCODER_H

#include <immintrin.h>
#include <pmmintrin.h>
#include <cstdio>
#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <filesystem>
#include <set>

namespace encoding {
    const int CHAR_LEN = 7;
    const int LONG_LEN = 63;

    inline int max(int a, int b) {return (a > b? a : b);}

    class Delta32 {
    protected:
        int pre;
        int maximum;
        //bool write;
        int posv;
    public:
        //Delta32(bool needWrite) { write=needWrite, posv=0;};
        Delta32() {maximum = 0;}
        void reset() {maximum = 0;}
        void encode_epi32(void* bufferx, long length) {
            long* buffer = (long*) bufferx;
            int tmp;
            pre = buffer[0];
            maximum = pre;
            for(long i=1;i<length;i++) {
                tmp = *(buffer + i);
                *(buffer + i) -= pre;
                int check = *(buffer + i);
                maximum = encoding::max(maximum, *(buffer + i));
                pre = tmp;
            }
        }
        int exportMaxDelta() {return maximum; }
        void decode(int* encodings, int length) {
            // serial decoding of delta codes.
            for(int i=1;i<length;i++) {
                encodings[i] += encodings[i-1];
            }
        }
        //int decode(char* encodings, int pos);
    };

    class BitPack {
    protected:
        int bit_len;
        long mask[65];
    public:
        void setBitLen(int length) {
            this->bit_len = length;
            mask[0]=1L;
            for(int i=1;i<=63;i++){
                mask[i] = (mask[i-1] << 1) | 1;
            }
        }
        void encode(int value, char* array, int &bias, long &pos) {
            int rest = this->bit_len;
            int x = value;
            char* curr = array + pos;
            while(rest) {
                int tmp = (CHAR_LEN - bias);
                if(rest <= tmp) {
                    //printf("%d\n", x << (tmp - rest + 1));
                    (*curr) |= (x << (tmp - rest));
                    bias += rest;
                    rest = 0;
                }
                else {
                    (*curr) |= (x >> (rest - tmp));
                    curr += 1;
                    (*curr) = 0;
                    rest -= tmp;
                    x &= mask[rest-1];
                    bias = 0;
                }
                if(bias == CHAR_LEN) {bias = 0; curr ++; }
                //printf("%d %d %ld\n", bias, rest, curr - array);
            }
            pos = curr - array;
        }

        void encode_l4epi32(int value, long* array, int &bias, long &pos) {
            int rest = this->bit_len;
            long x = value;
            long* curr = array + pos;
            while(rest) {
                int tmp = (LONG_LEN - bias);
                if(rest <= tmp) {
                    //printf("%d\n", x << (tmp - rest + 1));
                    (*curr) |= (((long)x) << (tmp - rest));
                    bias += rest;
                    rest = 0;
                }
                else {
                    (*curr) |= (x >> (rest - tmp));
                    curr += 1;
                    (*curr) = 0;
                    rest -= tmp;
                    x &= mask[rest-1];
                    bias = 0;
                }
                if(bias == LONG_LEN) {bias = 0; curr ++; }
                //printf("%d %d %ld\n", bias, rest, curr - array);
            }
            pos = curr - array;
        }

        

        void decode_int(int* array, char* bits, long &pos, int length) {
            int rest = 0;
            int added = 0;
            int bb = CHAR_LEN;
            int* curr = array + pos;
            int encp = 0;
            while(encp<length) {
                int tmp = bb + added;
                //printf("dec %d %d\n", added, (int) bits[encp]);
                if(tmp >= this->bit_len) {
                    rest <<= (this->bit_len - added);
                    rest |= (bits[encp] >> (bb - this->bit_len + added));
                    *curr = rest;
                    curr ++;
                    rest = 0; //bits[i] & mask[CHAR_LEN - this->bit_len + added - 1];
                    //added = 0;//CHAR_LEN - this->bit_len + added;
                    bb = bb - this->bit_len + added;
                    added = 0;
                    if(bb == 0) {encp ++; bb = CHAR_LEN;}
                    else {
                        bits[encp] &= mask[bb - 1];
                    }
                } else {
                    rest = (rest << bb) | bits[encp++];
                    added += bb;
                    bb = CHAR_LEN;
                }
            }
            pos = curr - array;
        }

        void decode_64xepi32(int* buffer, long* bits, long &encp, long &pos, int &bias, long length) {

        }

        
        void decode_l4epi32(int* array, long* bits, long &pos, int length) {
            long rest = 0;
            int added = 0;
            int bb = LONG_LEN;
            int* curr = array + pos;
            int encp = 0;
            while(encp<length) {
                int tmp = bb + added;
                //printf("dec %d %d\n", added, (int) bits[encp]);
                if(tmp >= this->bit_len) {
                    rest <<= (this->bit_len - added);
                    rest |= (bits[encp] >> (bb - this->bit_len + added));
                    *curr = rest;
                    curr ++;
                    rest = 0; //bits[i] & mask[CHAR_LEN - this->bit_len + added - 1];
                    //added = 0;//CHAR_LEN - this->bit_len + added;
                    bb = bb - this->bit_len + added;
                    added = 0;
                    if(bb == 0) {encp ++; bb = LONG_LEN; }
                    else {
                        bits[encp] &= mask[bb - 1];
                    }
                } else {
                    rest = (rest << bb) | bits[encp++];
                    added += bb;
                    bb = LONG_LEN;
                }
            }
            pos = curr - array;
        }

        void decode_vec8i(__m256i* array, char* bits, int pos) {

        }
    };

    class RunLength32 {
    protected:
        int last_val;
        int runs;
        int* val_array;
        int* rle_array;
        int posv;
        int maximum;
    public:
        RunLength32() {}
        RunLength32(int* vals, int* rles) {
            val_array =vals; rle_array =rles;
            last_val = -1; runs = 0;
            posv = 0;
            maximum = 0;
        }
        void encode(int value, int &pos) {
            if(posv) {
                runs = 1; last_val = value;
                return;
            }
            if(value == last_val) runs ++;
            else {
                *(val_array+pos) = last_val;
                *(rle_array+pos) = runs;
                maximum = encoding::max(runs, maximum);
                runs = 1;
                last_val = value;
                pos ++;
            }
        }
        void encodeFinalize(int &pos) {
            *(val_array+pos) = last_val;
            *(rle_array+pos) = runs;
            pos++;
            runs = 1;
        }
        void decode(int value, int runs, int* buffer, int &pos) {
            for(auto i=pos;i<runs+pos;i++) {
                buffer[i] = value;
            }
            pos = pos + runs;
            return;
        }
        int exportRLEMaximum() {return this->maximum; }
    };

    class Writer {
    protected:
        const char* file_loc;
    public:
        Writer(){}
        Writer(const char* outputFile) {this->file_loc = outputFile; }
        void writeAsLongArray(long* data, long length) {
            
            FILE* fptr = fopen(this->file_loc, "wb");
            //long len[] = {length};
            fwrite(&length, sizeof(long), 1, fptr);
            fwrite(data, length*sizeof(long), 1, fptr);
            fclose(fptr);
        }
        void writeAsCharArray(char* data, long length) {
            FILE* fptr = fopen(this->file_loc, "wb");
            length = length*sizeof(char)/sizeof(long) + 1;
            fwrite(&length, sizeof(char), 1, fptr);
            fwrite(data, length*sizeof(char), 1, fptr);
            fclose(fptr);
        }
    };

    class Reader {
    protected:
        const char* file_loc;
    public:
        Reader(const char* inputFile) {file_loc = inputFile;}
        long readIntoLongArray(long* buffer) {
            FILE* fptr = fopen(this->file_loc, "rb");
            long length;
            fread(&length, sizeof(long), 1, fptr);
            fread(buffer, length*sizeof(long), 1, fptr);
            fclose(fptr);
            return length;
        }
        long readIntoCharArray(char* buffer) {
            FILE* fptr = fopen(this->file_loc, "rb");
            long length;
            fread(&length, sizeof(long), 1, fptr);
            length = length*sizeof(long)/sizeof(char);
            fread(buffer, length*sizeof(char), 1, fptr);
            fclose(fptr);
            return length;
        }
    };

    namespace Math {
        int log2(long x) {
            int ans = 1;
            while(x) {
                x >>= 1;
                ans += 1;
            }
            return ans;
        }
    }

    namespace Encoder {
        const int BUFFER_SZ = 1000000;
        long BUFFER[BUFFER_SZ];
        int BUFF_POS = 0;

        void clear_glb() {
            memset(encoding::Encoder::BUFFER, 0, sizeof encoding::Encoder::BUFFER);
            BUFF_POS = 0;
        }

        int BUFF_DELTA[1000000];
        int BUFF_RLE[1000000];
        const int FQUALIFY = 1000;
        const int BUFF_SZ = 1000000;
        int RLE_POS = 0;
        void reset() {
            RLE_POS = 0;
            memset(encoding::Encoder::BUFF_DELTA, 0, sizeof(encoding::Encoder::BUFF_DELTA));
            memset(encoding::Encoder::BUFF_RLE, 0, sizeof(encoding::Encoder::BUFF_RLE));
        }
    }

    namespace Blocked {
        class LSMBlockReader {
        protected:
            std::string filePrefix;
            //std::string* rleFs;
            //std::string* dtaFs;
            std::set<int> usedFs;
            std::string nextDtaF;
            std::string nextRLEF;
            // +all candidate, curr.
            int dtablen;
            int rleblen;
            long* blkbuffer_dta;
            long* blkbuffer_rle; 
            long posdta; long posrle; 
            long tot; 
            int biasdta; int biasrle;
            const int buffsz = 1000000;
            __m256i delta[8];
            __m256i rle[8];
            long simdLen; long encpdta; long encprle;
        private:
            void init() {
                posdta = 0; posrle = 0;
                tot = 0; 
                biasdta = 0; biasrle = 0;
                encpdta = 0; encprle = 0;
                this->blkbuffer_dta = (long*) malloc(buffsz*sizeof(long));
                this->blkbuffer_rle = (long*) malloc(buffsz*sizeof(long));
            }
            bool hasNextBlock() {
                for(auto const &entry: std::filesystem::directory_iterator(this->filePrefix)) {
                    std::string fileName = entry.path().string();
                    int hash = std::stoi(fileName.substr(0, fileName.length()-4));
                    if(!usedFs.count(hash)) {
                        if(fileName.at(fileName.length()-1) == 'e') {
                            this->nextRLEF = fileName;
                            //usedFs.insert(this->nextRLEF);
                            std::string pref = fileName.substr(0, fileName.length()-4);
                            this->nextDtaF = pref.append(".dta");
                            //usedFs.insert(this->nextDtaF);
                            usedFs.insert(hash);
                        }
                        return true;
                    }
                }
                return false;
            }
            bool loadNextBlock() {
                init();
                FILE* fdta; FILE* frle;
                fdta=fopen(this->nextDtaF.c_str(), "rb");
                frle=fopen(this->nextRLEF.c_str(), "rb");
                if(!fdta || !frle) {
                    printf("Locating data blocks failed.\n");
                    exit(1);
                }
                long param;
                fread(&param, sizeof(long), 1, fdta);
                readHeader(param, 'd');
                fread(&param, sizeof(long), 1, fdta);
                readHeader(param, 'r');
                if(this->tot <= this->buffsz) {
                    fread(this->blkbuffer_dta, this->tot*sizeof(long), 1, fdta);
                    fread(this->blkbuffer_rle, this->tot*sizeof(long), 1, frle);
                } else {
                    printf("Allocation excceed.\n");
                    exit(1); // +handle.
                }
            }
            void readHeader(long params, char mode) {
                if(mode == 'd') {
                    this->dtablen = (int) (mode >> 32);
                    this->tot = mode & ((((long) 1) << 32) - 1);
                    this->simdLen = this->tot - this->tot%64;
                } else if (mode == 'r') {
                    this->rleblen = (int) (mode >> 32);
                } else {
                    printf("Mode detection failed.\n");
                    exit(1);
                }
            }
        public:
            LSMBlockReader() { init(); }
            LSMBlockReader(std::string prefix) {
                this->filePrefix = prefix; init();
            }
            int* restUnparallelDta() {
                int tmpd[64];int tmpr[64];
                encoding::BitPack bp = encoding::BitPack();
                bp.setBitLen(this->dtablen);
                bp.decode_64xepi32(
                    tmpd, this->blkbuffer_dta, this->encpdta,
                    this->posdta,this->biasdta,this->simdLen
                );
                bp.setBitLen(this->rleblen);
                bp.decode_64xepi32(
                    tmpr, this->blkbuffer_rle, this->encprle,
                    this->posrle,this->biasrle,this->simdLen
                );
                
            }
            __m256i* next8x32VectorDta() {
                if(this->posdta == this->simdLen) {
                    if(hasNextBlock()) {
                        loadNextBlock();
                    }
                }
                int dta_buff[128];
                encoding::BitPack bp = encoding::BitPack();
                
                bp.decode_64xepi32(
                    dta_buff,
                    this->blkbuffer_dta, 
                    this->encpdta, 
                    this->posdta,
                    this->biasdta, 
                    this->tot);
                __m256i* target = (__m256i*) dta_buff;
                for(auto i=0;i<8;i++)
                    delta[i] = _mm256_loadu_si256(target + i);
                return delta;
            }
            __m256i* next8x32VectorRLE() {
                int rle_buff[128];
                encoding::BitPack bp = encoding::BitPack();
                bp.decode_64xepi32(
                    rle_buff,
                    this->blkbuffer_rle, 
                    this->encprle, 
                    this->posrle,
                    this->biasrle, 
                    this->tot);
                __m256i* target = (__m256i*) rle_buff;
                for(auto i=0;i<8;i++)
                    rle[i] = _mm256_loadu_si256(target + i);
                return rle;
            }
        };
    }
    

    namespace BitPacked {
        int bias = 0; // bit bias inside a long/char 
        long pos = 1; // position of current writer of a bitpacked class.
        void reset() {
            bias = 0; pos = 1;
        }
    }

    class CSVEncoder {
    protected:
        int col; int QUALIFY; int BLK_SZ;
        long* buffer; long buffer_sz;
        long buffer_pos;
        encoding::BitPack bp;
        encoding::Delta32 delta;
        encoding::RunLength32 rle;
        int file_count = 0;
        std::string prefix;
        std::string postfix_delta, postfix_rle;
        std::string deltafile;
        std::string rlefile;
        encoding::Writer writev;
        encoding::Writer writer;
        long accumulated_points = 0;
        long all_encoded_data = 0;
        long all_used_bits = 0;
        int bit_dta=0; int bit_rle = 0;
    public:
        CSVEncoder(int col, const int QUALIFY, const int BLK_SZ, std::string encpath) {
            this->col = col;
            this->QUALIFY = QUALIFY;
            this->BLK_SZ = BLK_SZ;
            this->buffer = encoding::Encoder::BUFFER;
            this->buffer_sz = encoding::Encoder::BUFFER_SZ;
            this->buffer_pos = 0;
            this->bp = encoding::BitPack();
            this->delta = encoding::Delta32();
            this->rle = encoding::RunLength32(Encoder::BUFF_DELTA, Encoder::BUFF_RLE);
            encoding::Encoder::clear_glb();
            encoding::Encoder::reset(); 
            encoding::BitPacked::reset();
            this->prefix = encpath;
            if(mkdir(encpath.c_str(), S_IRWXU)) {
                //System();
                printf("Try to remove existing path.\n");
                if(remove(encpath.c_str())) {
                    printf("Create target encoding files failed.\n");
                    exit(1);
                }
            }
            this->file_count = 0;
            this->postfix_delta = ".dta"; this->postfix_rle = ".rle";
            this->deltafile = prefix;
            this->deltafile.append(std::to_string(file_count)).append(postfix_delta);
            this->rlefile = prefix;
            this->rlefile.append(std::to_string(file_count)).append(postfix_rle);
            this->writev = encoding::Writer(deltafile.c_str());
            this->writer = encoding::Writer(rlefile.c_str());
            encoding::Encoder::clear_glb();
        }
        void reset_flush() {
            this->buffer_pos = 0;
            this->bp = encoding::BitPack();
            this->delta = encoding::Delta32();
            this->rle = encoding::RunLength32(Encoder::BUFF_DELTA, Encoder::BUFF_RLE);
            encoding::Encoder::reset(); 
            encoding::BitPacked::reset();
            encoding::Encoder::clear_glb();
            this->file_count ++;
            this->deltafile = prefix;
            this->deltafile.append(std::to_string(file_count)).append(postfix_delta);
            this->rlefile = prefix;
            this->rlefile.append(std::to_string(file_count)).append(postfix_rle);
            this->writev = encoding::Writer(deltafile.c_str());
            this->writer = encoding::Writer(rlefile.c_str());
            this->accumulated_points = 0;
            this->buffer_pos = 0;
        }
        void flush() {
            this->delta.encode_epi32(this->buffer, this->buffer_pos);
            for(int i=0;i<this->buffer_pos;i++) {
                this->rle.encode(this->buffer[i], encoding::Encoder::RLE_POS);
            }
            this->rle.encodeFinalize(encoding::Encoder::RLE_POS);
            int delta_blen = encoding::Math::log2(this->delta.exportMaxDelta());
            int rle_blen = encoding::Math::log2(this->rle.exportRLEMaximum());
            
            // encode delta.
            bp.setBitLen(delta_blen);
            this->buffer[0] = (((long) delta_blen) << 32) | (this->buffer_pos);
            for(int i=0; i< encoding::Encoder::RLE_POS;i++) {
                bp.encode_l4epi32(encoding::Encoder::BUFF_DELTA[i], this->buffer, encoding::BitPacked::bias, encoding::BitPacked::pos);
            }
            long length = encoding::BitPacked::pos;
            this->all_used_bits += delta_blen*length;
            writev.writeAsLongArray(this->buffer, length);
            encoding::BitPacked::reset();
            
            //encode rle
            bp.setBitLen(rle_blen);
            this->buffer[0] = (((long) rle_blen) << 32) | (this->buffer_pos);
            for(int i=0;i<encoding::Encoder::RLE_POS;i++) {
                bp.encode_l4epi32(encoding::Encoder::BUFF_RLE[i], this->buffer, encoding::BitPacked::bias, encoding::BitPacked::pos);
            }
            length = encoding::BitPacked::pos;
            this->all_used_bits += rle_blen*length;
            writer.writeAsLongArray(this->buffer, length);
            encoding::BitPacked::reset();
            
            // reset flush
            this->all_encoded_data += this->buffer_pos;
            this->bit_dta = encoding::max(this->bit_dta, delta_blen);
            this->bit_rle = encoding::max(this->bit_rle, rle_blen);
            this->reset_flush();
        }
        int exportDtaBit() { return this->bit_dta; }
        int exportRleBit() { return this->bit_rle; }
        double finalize() {
            return 64.0*((double)this->all_encoded_data)*1.0/((double) this->all_used_bits);
        }
        void encode(const char* file) {
            bool line0 = true;
            std::string line;
            std::ifstream infile(file);
            this->buffer_pos = 0;
            while(std::getline(infile, line)) {
                if(line0) {line0 =false; continue;}
                std::vector<long> tmp;
                std::stringstream ss(line);
                int p = 0;
                for(double i; ss>>i;p++) {
                    if(p == col) tmp.push_back((long)(i * QUALIFY));
                    else tmp.push_back((long) i);
                    if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n' || ss.peek() == '-' || ss.peek() == ':') {
                        ss.ignore();
                    }
                }
                if(this->buffer_pos == BLK_SZ - 1) {
                    flush();
                    this->buffer_pos = 0;
                }
                this->buffer[this->buffer_pos++] = tmp.at(col);
            }
            if(this->buffer_pos!=0) {
                flush();
            }
        }
    };
}

#endif
//#pragma once
