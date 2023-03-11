#include "optcntVFq.hpp"

//int blk_sz = 100;
//int buf_pos = 0;
int bufferd[100010];
int bufferr[100010];
long tot_T = 0;
long tot_F = 0;
int st = -1;
//int posv = 1;
//double selectivity = 0.10;
long cntTup = 0;

const __m256i INV = _mm256_setr_epi32(3, 2, 1, 0, 7, 6, 5, 4);

class sboost {
private:
    int lastest = 0;
    long count = 0;
    char op1, op2;
    int c1, c2;
public:
    void decode_and_filter(int* data);
    void init(int c1, int c2, char op1, char op2);
    long extract();
    void init_lastest(int start);
};



__m256i cumsum(__m256i b) {
    __m256i bp = permute1(b);
    __m256i s1 = _mm256_hadd_epi32(b, bp);
    __m256i s2 = permute2(s1);
    __m256i s3 = _mm256_hadd_epi32(s1, s2);
    __m256i s4 = permute4(s3);
    __m256i result = _mm256_add_epi32(s3, s4);
    return _mm256_permutevar8x32_epi32(result, INV);
}

void sboost::decode_and_filter(int* data) {
    __m256i b = _mm256_load_si256((__m256i*) data);
    __m256i sum = cumsum(b);
    __m256i decoded = _mm256_add_epi32(_mm256_set1_epi32(this->lastest), sum);
    this->lastest = _mm256_extract_epi32(decoded, 7);
    __m256i res1 = predicates(_mm256_sub_epi32(decoded, _mm256_set1_epi32(c1)), ZERO, this->op1);
    __m256i res2 = predicates(_mm256_sub_epi32(decoded, _mm256_set1_epi32(c2)), ZERO, this->op2);
    __m256i merge = _mm256_and_epi32(_mm256_and_epi32(res1, res2), POS1);
    this->count = this->count + sumAll(merge);
}


void sboost::init(int c1, int c2, char op1, char op2) {
    this->c1 = c1;
    this->c2 = c2;
    this->op1 = op1;
    this->op2 = op2;
}

void sboost::init_lastest(int start) {
    this->lastest = start;
}

long sboost::extract() {
    return this->count;
}

int process_sboost(int* delta, int* rle, bool isSingleBitRLE, int length,
    char op1, int c1, char op2, int c2, bool enablePrune) {
    sboost bst = sboost();
    bst.init(c1, c2, op1, op2);
    int rle_rest = 0;
    int currd = -1;
    int data[8], dpos = 0;
    for (int i = 0; i < length; i += 1) {
        currd = bufferd[i];
        rle_rest = bufferr[i];
        while (rle_rest) {
            data[dpos++] = currd;
            rle_rest--;
            if (dpos == 8) {
                bst.decode_and_filter(data);
                dpos = 0;
            }
        }
    }
    if (dpos != 0) {
        for (int i = dpos; i < 8; i++) {
            data[i] = 0;
        }
    }
    bst.decode_and_filter(data);
    return bst.extract();
}

__m256i sumDelta(__m256i* deltas, __m256i* rles) {
    //__m256i vc1 = _mm256_set1_epi32(c1);
    //__m256i vc2 = _mm256_set1_epi32(c2);
    //int* delta_inc = deltas + 1;
    //int* rle_inc = rles + 1;
    __m256i d = _mm256_load_si256(deltas);
    //__m256i d0 = _mm256_load_si256((__m256i*) deltas);
    __m256i rle = _mm256_load_si256(rles);
    //__m256i rle0 = _mm256_load_si256((__m256i*) rles);

    __m256i bias;
    __m256i bias0;
    bias = _mm256_mullo_epi32(d, rle);
    bias0 = permute1(bias);
    __m256i inc = _mm256_add_epi32(bias, bias0);
    __m256i inc00 = permute2(inc);
    __m256i psum = _mm256_add_epi32(inc, inc00);
    __m256i psumH = permute4(psum);
    __m256i sum = _mm256_add_epi32(psum, psumH);
    return sum;
}

__m256i sumDelta(__m256i* deltas) {
    //__m256i vc1 = _mm256_set1_epi32(c1);
    //__m256i vc2 = _mm256_set1_epi32(c2);
    //int* delta_inc = deltas + 1;
    //int* rle_inc = rles + 1;
    __m256i d = _mm256_load_si256(deltas);
    //__m256i d0 = _mm256_load_si256((__m256i*) deltas);
    //__m256i rle = _mm256_load_si256(rles);
    //__m256i rle0 = _mm256_load_si256((__m256i*) rles);

    __m256i bias = d;
    __m256i bias0;
    //bias = _mm256_mullo_epi32(d, rle);
    bias0 = permute1(bias);
    __m256i inc = _mm256_add_epi32(bias, bias0);
    __m256i inc00 = permute2(inc);
    __m256i psum = _mm256_add_epi32(inc, inc00);
    __m256i psumH = permute4(psum);
    __m256i sum = _mm256_add_epi32(psum, psumH);
    return sum;
}


__m256i buffer[1000100];
int base[250010];
int rleBit = 1;
int maxd = -0x3f3f3f3f, mind = 0x3f3f3f3f;
int maxv = -1, minv = -1;
int maxr = -1, minr = -1;

int process_scalar(int* delta, int* rle, bool isSingleBitRLE, int length,
    char op1, int c1, char op2, int c2, bool enablePrune) {
    int acc = 0;
    int count = 0;
    for (int i = 0; i < length; i++) {
        int rest = rle[i];
        while (rest) {
            rest--;
            acc += delta[i];
            if (predicate_scalar(acc, op1, c1) && predicate_scalar(acc, op2, c2)) {
                count += 1;//random_i32();
            }
        }
    }
    return count;
}

// now only singleBit
int process_cntVF(int* delta, int* rle, bool isSingleBitRLE, int length,
    char op1, int c1, char op2, int c2, bool enablePrune) {
    int len = length;
    length = length / 8;
    if (length > 1000000) {
        printf("Mem not enough");
        return -1;
    }
    __m256i* delta_src = (__m256i*) delta;
    __m256i* rle_src = (__m256i*) rle;
    int maxCompGroup = length / 4;

#pragma prefetch
#pragma vector
    if (isSingleBitRLE) {
#pragma prefetch
#pragma omp parallel for
        for (int i = 0; i < length; i++) {
            buffer[i] = sumDelta(delta_src + i);
        }

    }
    else {
#pragma prefetch
#pragma omp parallel for
        for (int i = 0; i < length; i++) {
            buffer[i] = sumDelta(delta_src + i, rle_src + i);
        }
    }
#pragma prefetch
#pragma omp parallel for
    for (int i = 0; i < maxCompGroup; i++) {
        buffer[i * 4 + 1] = _mm256_add_epi32(buffer[i * 4 + 1], _mm256_permutevar8x32_epi32(buffer[i * 4], POS7));
        buffer[i * 4 + 2] = _mm256_add_epi32(buffer[i * 4 + 2], _mm256_permutevar8x32_epi32(buffer[i * 4 + 1], POS7));
        buffer[i * 4 + 3] = _mm256_add_epi32(buffer[i * 4 + 3], _mm256_permutevar8x32_epi32(buffer[i * 4 + 2], POS7));
        base[i + 1] = _mm256_extract_epi32(buffer[i * 4 + 3], 7);
    }

#pragma prefetch
    for (int i = 2; i < maxCompGroup; i++) {
        base[i] += base[i - 1];
        if (enablePrune) {
            int lb = base[i] + mind * (maxr) * (length - 4 * i) * 8;//*(selectivity);
            int ub = base[i] + maxd * (maxr) * (length - 4 * i) * 8;
            if (lb >= c1 && ub <= c2) {
                int xcnt = 0;
                //printf("Prune:%d %d", i, maxCompGroup);
#pragma prefetch
                for (int i = 0; i < len; i++) {
                    xcnt += rle[i];
                }
                return xcnt;
            }
            if (lb > c2 || ub < c1) {
                //printf("Prune:%d %d", i, maxCompGroup);
                maxCompGroup = i;
                break;
            }
        }
    }
#pragma omp parallel for
    for (int i = 0; i < maxCompGroup; i++) {
        __m256i vc1 = _mm256_set1_epi32(c1 - base[i]);
        __m256i vc2 = _mm256_set1_epi32(c2 - base[i]);
        __m256i tmp[4];
#pragma omp parallel for
        for (int j = 0; j < 4; j++) {
            __m256i left = _mm256_and_si256(predicates(buffer[i * 4 + j], vc1, op1),
                predicates(buffer[i * 4 + j], vc2, op2));
            tmp[j] = _mm256_and_si256(POS1, left);
        }
        base[i] = sumAll(
            _mm256_add_epi32(_mm256_add_epi32(tmp[0], tmp[1]),
                _mm256_add_epi32(tmp[2], tmp[3])));
    }
#pragma prefetch
    for (int i = 2; i < maxCompGroup; i++) {
        base[i] += base[i - 1];
    }
    return base[maxCompGroup - 1];
}



int process_cntVFnp(int* delta, int* rle, bool isSingleBitRLE, int length,
    char op1, int c1, char op2, int c2, bool enablePrune) {
    int len = length;
    length = length / 8;
    if (length > 1000000) {
        printf("Mem not enough");
        return -1;
    }
    __m256i* delta_src = (__m256i*) delta;
    __m256i* rle_src = (__m256i*) rle;
    int maxCompGroup = length / 4;

#pragma prefetch
#pragma vector
    if (isSingleBitRLE) {
#pragma prefetch
#pragma omp parallel for
        for (int i = 0; i < length; i++) {
            buffer[i] = sumDelta(delta_src + i);
        }

    }
    else {
#pragma prefetch
#pragma omp parallel for
        for (int i = 0; i < length; i++) {
            buffer[i] = sumDelta(delta_src + i, rle_src + i);
        }
    }
#pragma prefetch
#pragma omp parallel for
    for (int i = 0; i < maxCompGroup; i++) {
        buffer[i * 4 + 1] = _mm256_add_epi32(buffer[i * 4 + 1], _mm256_permutevar8x32_epi32(buffer[i * 4], POS7));
        buffer[i * 4 + 2] = _mm256_add_epi32(buffer[i * 4 + 2], _mm256_permutevar8x32_epi32(buffer[i * 4 + 1], POS7));
        buffer[i * 4 + 3] = _mm256_add_epi32(buffer[i * 4 + 3], _mm256_permutevar8x32_epi32(buffer[i * 4 + 2], POS7));
        base[i + 1] = _mm256_extract_epi32(buffer[i * 4 + 3], 7);
    }

#pragma prefetch
    for (int i = 2; i < maxCompGroup; i++) {
        base[i] += base[i - 1];
    }
#pragma omp parallel for
    for (int i = 0; i < maxCompGroup; i++) {
        __m256i vc1 = _mm256_set1_epi32(c1 - base[i]);
        __m256i vc2 = _mm256_set1_epi32(c2 - base[i]);
        __m256i tmp[4];
#pragma omp parallel for
        for (int j = 0; j < 4; j++) {
            __m256i left;
            if (op2 == 'x') {
                left = predicates(buffer[i * 4 + j], vc1, op1);
            }
            else {
                left = _mm256_and_si256(
                    predicates(buffer[i * 4 + j], vc1, op1),
                    predicates(buffer[i * 4 + j], vc2, op2));
            }
            tmp[j] = _mm256_and_si256(POS1, left);
        }
        base[i] = sumAll(
            _mm256_add_epi32(_mm256_add_epi32(tmp[0], tmp[1]),
                _mm256_add_epi32(tmp[2], tmp[3])));
    }
#pragma prefetch
    for (int i = 1; i < maxCompGroup; i++) {
        base[i] += base[i - 1];
    }
    return base[maxCompGroup - 1];
}

long tot_F2 = 0, tot_F3 = 0;
int bit_delta = 0, bit_rle = 0;

long test_sys() {
    //time count 
    // value 
    // time --data--> value
    return -1;
}

int buffermst[100000100];
int mst_plc = 0;

long decompress_MST(int* delta, int* rle, int st, int length) {
    int acc = st;
    int count = 0;
    buffermst[mst_plc++] = acc;
    for (int i = 0; i < length; i++) {
        int rest = rle[i];
        while (rest) {
            rest--;
            acc += delta[i];
            buffermst[mst_plc++] = acc;
        }
    }
    return 0;
}

int query_quantile(double quant) {
    std::sort(buffermst, buffermst + mst_plc);
    return buffermst[(int)( quant*mst_plc)];
}

long test(std::string input, int dup, int blk_sz, double selectivity,
    int posv) { //int rle_bit_restrict
    std::ifstream infile(input);
    //printf("<");
    int blkID = 0;
    //std::ofstream outFile(output_path.append("" + blk));
    std::string line;
    int pos = 0;
    int pre = -1;
    int delta, rle = 1;
    bool line0 = true;
    bool proceeded = false;
    int bufferPos = 0;
    while (std::getline(infile, line)) {
        //printf("<");
        if (line0) {
            line0 = false;
            continue;
        }
        std::stringstream ss(line);
        std::vector<int> tmp;
        for (double i; ss >> i; ) {
            tmp.push_back((int)(i * QUALIFY));
            if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n' || ss.peek() == '-' || ss.peek() == ':') {
                ss.ignore();
            }
        }
        if (tmp.size() < posv + 1) continue;
        //printf("%d", tmp.at(1));
        if (pos == 0) {
            //std::string tmp = line.substr(line.find(","), line.length());

            pre = tmp.at(posv);
            st = pre;
            maxv = minv = pre;
        }
        else if (pos == 1) {
            int tmpv = tmp.at(posv);
            delta = tmpv - pre;
            //delta = delta >> delta_div;
            pre = tmpv;
            rle = 1;// std::max(1, std::abs(random_i32()));
            //rle = std::max(1, std::abs(random_i32()));;
            maxd = mind = delta;
            maxr = minr = rle;
        }
        else {
            int tmpv = tmp.at(posv);
            int tmpd = tmpv - pre;
            //tmpd = tmpd >> delta_div;
            if (delta == tmpd) {
                //&& rle < (1<<rle_bit_restrict)
                rle += 1;// std::max(1, std::abs(random_i32()));
                //rle += std::max(1, std::abs(random_i32()));;
            }
            else {
                cntTup += 1;//rle;
                bufferd[bufferPos] = delta;
                bufferr[bufferPos++] = rle;
                maxr = std::max(maxr, rle);
                minr = std::min(minr, rle);

                rle = 1; // std::max(1, std::abs(random_i32()));
                //rle = std::max(1, std::abs(random_i32()));;
                delta = tmpd;
                maxd = std::max(maxd, delta);
                mind = std::min(mind, delta);
            }
            pre = tmpv;
        }
        pos++;
        maxv = std::max(maxv, pre);
        minv = std::min(minv, pre);
        if (bufferPos == blk_sz - 1) {
            proceeded = true;
            bit_delta = std::max(bit_delta, log2(std::max(maxd, std::abs(mind))));
            bit_rle = std::max(bit_rle, log2(std::max(maxr, std::abs(minr))));
            int c2 = (int)(minv + (maxv - minv) * (selectivity)) - st;
            int c1 = (int)(minv)-st;
            //int c1 = -5*QUALIFY - st;
            //int c2 = -5*QUALIFY + (int)(10*QUALIFY*selectivity) - st;
            //printf("Start testing.");
            int tmp1, tmp2, tmp3;
            long start_ts;

            std::ofstream ofile("./data/tmp.csv");
            start_ts = timeSinceEpochMillisec();
            for (int i = 0; i < dup; i++) {
                for (int j = 0; j < blk_sz; j++) bufferd[j] += 100;
                tmp3 = process_scalar(bufferd, bufferr, true, bufferPos + 1,
                    '[', c1, '<', c2, false);
                ofile << tmp3;
                //printf("Scalar: %d, ", tmp3);
            }
            //printf("Scalar: %d, ", tmp3);
            tot_F3 += (timeSinceEpochMillisec() - start_ts) * dup;


            for (int i = 0; i < dup; i++) {
                start_ts = timeSinceEpochMillisec();
                tmp1 = process_sboost(bufferd, bufferr, true, bufferPos + 1,
                    '[', c1, '<', c2, true);
                tot_T += (timeSinceEpochMillisec() - start_ts);
                ofile << tmp1;
            }
            //printf("SBoost: %d, ", tmp1);



            for (int i = 0; i < dup; i++) {
                start_ts = timeSinceEpochMillisec();
                tmp2 = process_cntVF(bufferd, bufferr, true, bufferPos + 1,
                    '[', c1, '<', c2, true);
                tot_F += (timeSinceEpochMillisec() - start_ts);
                memset(base, 0, sizeof(base));
                ofile << tmp2;
            }
            //printf("cntvf: %d, ", tmp2);



            for (int i = 0; i < dup; i++) {
                start_ts = timeSinceEpochMillisec();
                tmp3 = process_cntVFnp(bufferd, bufferr, true, bufferPos + 1,
                    '[', c1, '<', c2, false);
                tot_F2 += (timeSinceEpochMillisec() - start_ts);
                memset(base, 0, sizeof(base));
                memset(buffer, 0, sizeof(buffer));
                ofile << tmp3;
            }
            //printf("cntvf-p: %d, ", tmp3);
            ofile.close();

            //printf("%d,%d: count %d-%d \n", c1, c2, tmp1, tmp2);
            pos = 0;
            st = pre;
            bufferPos = 0;
            //break;
        }
    }
    infile.close();
    if (!proceeded) {
        proceeded = true;
        bit_delta = std::max(bit_delta, log2(std::max(maxd, std::abs(mind))));
        bit_rle = std::max(bit_rle, log2(std::max(maxr, std::abs(minr))));
        int c2 = (int)(minv + (maxv - minv) * (selectivity)) - st;
        int c1 = (int)(minv)-st;
        //int c1 = -5*QUALIFY - st;
        //int c2 = -5*QUALIFY + (int)(10*QUALIFY*selectivity) - st;
        //printf("Start testing.");
        int tmp1, tmp2, tmp3;
        long start_ts;

        std::ofstream ofile("./data/tmp.csv");
        start_ts = timeSinceEpochMillisec();
        for (int i = 0; i < dup; i++) {
            for (int j = 0; j < blk_sz; j++) bufferd[j] += 100;
            tmp3 = process_scalar(bufferd, bufferr, true, bufferPos + 1,
                '[', c1, '<', c2, false);
            ofile << tmp3;
            //printf("Scalar: %d, ", tmp3);
        }
        //printf("Scalar: %d, ", tmp3);
        tot_F3 += (timeSinceEpochMillisec() - start_ts) * dup;


        for (int i = 0; i < dup; i++) {
            start_ts = timeSinceEpochMillisec();
            tmp1 = process_sboost(bufferd, bufferr, true, bufferPos + 1,
                '[', c1, '<', c2, true);
            tot_T += (timeSinceEpochMillisec() - start_ts);
            ofile << tmp1;
        }
        //printf("SBoost: %d, ", tmp1);



        for (int i = 0; i < dup; i++) {
            start_ts = timeSinceEpochMillisec();
            tmp2 = process_cntVF(bufferd, bufferr, true, bufferPos + 1,
                '[', c1, '<', c2, true);
            tot_F += (timeSinceEpochMillisec() - start_ts);
            memset(base, 0, sizeof(base));
            ofile << tmp2;
        }
        //printf("cntvf: %d, ", tmp2);



        for (int i = 0; i < dup; i++) {
            start_ts = timeSinceEpochMillisec();
            tmp3 = process_cntVFnp(bufferd, bufferr, true, bufferPos + 1,
                '[', c1, '<', c2, false);
            tot_F2 += (timeSinceEpochMillisec() - start_ts);
            memset(base, 0, sizeof(base));
            memset(buffer, 0, sizeof(buffer));
            ofile << tmp3;
        }
        //printf("cntvf-p: %d, ", tmp3);
        ofile.close();

        //printf("%d,%d: count %d-%d \n", c1, c2, tmp1, tmp2);
        pos = 0;
        st = pre;
        bufferPos = 0;
    }
    //bufferd[buf_pos] = delta;
    //bufferr[buf_pos++] = rle;
    //bufferd2[buf_pos] = (short)delta;
    //bufferr2[buf_pos++] = (short)rle;
    //write(label);
    return tot_T;
}



long throughput(long time, int round) {
    return (long)cntTup * round / (time / (1.0 * std::pow(10, 9)));
}

int query_filter_cntvf(std::string input, int dup, int blk_sz, double selectivity,
    int posv1, int posv2) {
    test(input, dup, blk_sz, selectivity, posv1);
    test(input, dup, blk_sz, selectivity, posv2);
    return 0;
}

void run3(){
    int blk_sz = 1600;
    double b = 0.1;
    int round = 1; // 100 for sine
    int posv;
    for (int i =9; i <= 10; i++) {
        // various selectivity
        //double selectivity = b*(i*1.0);
        double selectivity = 0.1*(i*1.0);

        // keep not changed
        tot_T = tot_F = tot_F2 = tot_F3 = 0;

        // change here for different dataset
        //posv = 1;// 1:iot.climate, sine, tpch
        posv = 7;//, 8,9,10,... // gas, wind
        //posv = 0; QUALIFY =1;// timestamp
        //QUALIFY = (1<<i);
        //int rle_bit_rest = 15;
        //int delta_div = i;

        //printf("%lf\n", selectivity);
        //test("./data/sin-1.csv", round, blk_sz, selectivity, delta_div);

        // // use following line for latebncy
        // printf("%d %ld %ld %ld %ld\n", i*10, tot_T, tot_F, tot_F2, tot_F3);
        // //selectivity = b*i;
        // round*=10;

        query_filter_cntvf("./data/iot.climate.csv", 1, blk_sz, selectivity, 0, 1);

        // use following line for throughput
        //tot_T = throughput(tot_T, round);
        //tot_F = throughput(tot_F, round);
        //tot_F2 = throughput(tot_F2, round);
        //tot_F3 = throughput(tot_F3, round);
        printf("%d %ld %ld %ld %ld\n", i*10, tot_T, tot_F, tot_F2, tot_F3);
        //printf("%d %ld %ld %ld %ld\n", blk_sz, tot_T, tot_F, tot_F2, tot_F3);
        //blk_sz *= 5;
        //the output printed result > copy to .csv in overleaf project: ./exp/xxx/data/???.csv
    }
}

int rank_cntvf(int x, bool singleBit, int length, bool less) {
    int rk = 0;
    if (less) {
        rk = process_cntVFnp(bufferd, bufferr, singleBit, length,
            '<', x, 'x', -1, false);
    }
    else {
        rk = process_cntVFnp(bufferd, bufferr, singleBit, length,
            '[', x, 'x', -1, false);
    }
    return rk;
}



long quantile_cntvf(double quant, std::string input, int blk_sz, int posv) {
    
    bool less = true;
    if (quant > 0.5) {
        less = false; quant = 1 - quant;
    }
    std::ifstream infile(input);
    //printf("<");
    int blkID = 0;
    int glb_max, glb_min = -1;
    //std::ofstream outFile(output_path.append("" + blk));
    std::string line;
    int pos = 0;
    int pre = -1;
    int delta, rle = 1;
    bool line0 = true;
    bool proceeded = false;
    int bufferPos = 0;
    while (std::getline(infile, line)) {
        if (line0) {
            line0 = false;
            continue;
        }
        std::stringstream ss(line);
        std::vector<int> tmp;
        for (double i; ss >> i; ) {
            tmp.push_back((int)(i * QUALIFY));
            if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n' || ss.peek() == '-' || ss.peek() == ':') {
                ss.ignore();
            }
        }
        if (tmp.size() < posv + 1) continue;
        int curr = (int) QUALIFY * tmp.at(posv);
        if (pos == 0) {
            glb_min = glb_max = curr;
        }
        else {
            glb_min = std::min(glb_min, curr);
            glb_max = std::max(glb_max, curr);
        }
        pos++;
    }
    int glb_sz = pos;
    int find_plc = (int) (quant * 1.0 * glb_sz);
    infile.close();
    line0 = true; pos = 0;
    int prev_quant = 0, prev_x= -1;
    int curr_quant = 0;
    int lb = glb_min, ub = glb_max;
    int curr_x = lb; // (int)(glb_min + glb_max) / 2;
    // counting time start here:
    long time_cost=0;
    while (true) {
        if (prev_x == curr_x || curr_quant == find_plc) {
            return curr_x;
        }
        else {
            prev_x = curr_x; 
            if (curr_quant < find_plc) {
                lb = curr_x;
                curr_x = (int)(curr_x + ub) / 2;
            }
            else {
                ub = curr_x;
                curr_x = (int)(curr_x + lb) / 2;
            }
            prev_quant = curr_quant;
        }
        std::ifstream infile(input);
        while (std::getline(infile, line)) {
            //printf("<");
            if (line0) {
                line0 = false;
                continue;
            }
            long exclude_ts = timeSinceEpochMillisec();
            std::stringstream ss(line);
            std::vector<int> tmp;
            for (double i; ss >> i; ) {
                tmp.push_back((int)(i * QUALIFY));
                if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n' || ss.peek() == '-' || ss.peek() == ':') {
                    ss.ignore();
                }
            }
            if (tmp.size() < posv + 1) continue;
            //printf("%d", tmp.at(1));
            if (pos == 0) {
                //std::string tmp = line.substr(line.find(","), line.length());
                pre = tmp.at(posv);
                st = pre;
                maxv = minv = pre;
            }
            else if (pos == 1) {
                int tmpv = tmp.at(posv);
                delta = tmpv - pre;
                pre = tmpv;
                rle = 1;// std::max(1, std::abs(random_i32()));
                //rle = std::max(1, std::abs(random_i32()));;
                maxd = mind = delta;
                maxr = minr = rle;
            }
            else {
                int tmpv = tmp.at(posv);
                int tmpd = tmpv - pre;
                if (delta == tmpd) {
                    //&& rle < (1<<rle_bit_restrict)
                    rle += 1;// std::max(1, std::abs(random_i32()));
                    //rle += std::max(1, std::abs(random_i32()));;
                }
                else {
                    cntTup += 1;//rle;
                    bufferd[bufferPos] = delta;
                    bufferr[bufferPos++] = rle;
                    maxr = std::max(maxr, rle);
                    minr = std::min(minr, rle);

                    rle = 1; // std::max(1, std::abs(random_i32()));
                    //rle = std::max(1, std::abs(random_i32()));;
                    delta = tmpd;
                    maxd = std::max(maxd, delta);
                    mind = std::min(mind, delta);
                }
                pre = tmpv;
            }
            pos++;
            maxv = std::max(maxv, pre);
            minv = std::min(minv, pre);
            long start_ts = timeSinceEpochMillisec();
            if (bufferPos == blk_sz - 1) {
                pos = 0; 
                
                curr_quant = rank_cntvf(curr_x, maxr < 2, bufferPos, less);
                bufferPos = 0;
            }
            time_cost += timeSinceEpochMillisec() - start_ts;
        }
        infile.close();
    }
    return time_cost;
}

long quantile_mst(double quant, std::string input, int blk_sz, int posv) {
    
    int pos = 0;
    int pre = -1;
    int delta, rle = 1;
    bool line0 = true;
    bool proceeded = false;
    int bufferPos = 0;
    //int prev_quant = 0, prev_x= -1;
    //int curr_quant = 0;
    //int lb = glb_min, ub = glb_max;
    //int curr_x = lb; // (int)(glb_min + glb_max) / 2;
    // counting time start here:
    long time_cost = 0;
    int st;
    std::string line;
    std::ifstream infile(input);
        while (std::getline(infile, line)) {
            //printf("<");
            if (line0) {
                line0 = false;
                continue;
            }
            long exclude_ts = timeSinceEpochMillisec();
            std::stringstream ss(line);
            std::vector<int> tmp;
            for (double i; ss >> i; ) {
                tmp.push_back((int)(i * QUALIFY));
                if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n' || ss.peek() == '-' || ss.peek() == ':') {
                    ss.ignore();
                }
            }
            if (tmp.size() < posv + 1) continue;
            //printf("%d", tmp.at(1));
            if (pos == 0) {
                //std::string tmp = line.substr(line.find(","), line.length());
                pre = tmp.at(posv);
                st = pre;
                maxv = minv = pre;
            }
            else if (pos == 1) {
                int tmpv = tmp.at(posv);
                delta = tmpv - pre;
                pre = tmpv;
                rle = 1;// std::max(1, std::abs(random_i32()));
                //rle = std::max(1, std::abs(random_i32()));;
                maxd = mind = delta;
                maxr = minr = rle;
            }
            else {
                int tmpv = tmp.at(posv);
                int tmpd = tmpv - pre;
                if (delta == tmpd) {
                    //&& rle < (1<<rle_bit_restrict)
                    rle += 1;// std::max(1, std::abs(random_i32()));
                    //rle += std::max(1, std::abs(random_i32()));;
                }
                else {
                    cntTup += 1;//rle;
                    bufferd[bufferPos] = delta;
                    bufferr[bufferPos++] = rle;
                    maxr = std::max(maxr, rle);
                    minr = std::min(minr, rle);

                    rle = 1; // std::max(1, std::abs(random_i32()));
                    //rle = std::max(1, std::abs(random_i32()));;
                    delta = tmpd;
                    maxd = std::max(maxd, delta);
                    mind = std::min(mind, delta);
                }
                pre = tmpv;
            }
            pos++;
            maxv = std::max(maxv, pre);
            minv = std::min(minv, pre);
            if (bufferPos == blk_sz - 1) {
                pos = 0; 
                printf("decoded\n");
                //start_ts += timeSinceEpochMillisec() - exclude_ts;
                long start_ts = timeSinceEpochMillisec();
                decompress_MST(bufferd, bufferr, st, bufferPos);
                time_cost += timeSinceEpochMillisec() - start_ts;
                bufferPos = 0;
                proceeded = true;
            }
        }
        infile.close();
    if (!proceeded) {
        pos = 0; 
        printf("decoded\n");
        //start_ts += timeSinceEpochMillisec() - exclude_ts;
        long start_ts = timeSinceEpochMillisec();
        decompress_MST(bufferd, bufferr, st, bufferPos);
        time_cost += timeSinceEpochMillisec() - start_ts;
        bufferPos = 0;
    }
    std::ofstream ofile("./data/tmp.csv");
    long start_ts = timeSinceEpochMillisec();
    int val = query_quantile(quant);
    time_cost += timeSinceEpochMillisec() - start_ts;
    ofile << val;
    mst_plc = 0;
    return time_cost;
}

void run4(){
    int blk_sz = 1600;
    double b = 0.1;
    int round = 1; // 100 for sine
    int posv;
    for (int i =1; i <= 10; i++) {
        // various selectivity
        //double selectivity = b*(i*1.0);
        double selectivity = 0.1*(i*1.0);

        // keep not changed
        //tot_T = tot_F = tot_F2 = tot_F3 = 0;

        // change here for different dataset
        //posv = 1;// 1:iot.climate, sine, tpch
        //posv = 7;//, 8,9,10,... // gas, wind
        posv = 0; QUALIFY =1;// timestamp
        //QUALIFY = (1<<i);
        //int rle_bit_rest = 15;
        //int delta_div = i;

        //printf("%lf\n", selectivity);
        long time = quantile_mst(selectivity, "./data/iot.climate.csv", blk_sz, posv);
        //long time2 = quantile_cntvf(selectivity, "./data/sin-1.csv", blk_sz, posv);

        //test("./data/sin-1.csv", round, blk_sz, selectivity, delta_div);
        printf("%d %ld\n", 10*i, time);

        //query_filter_cntvf("./data/iot.climate.csv", 1, blk_sz, selectivity, 0, 1);
        //blk_sz *= 5;
        //the output printed result > copy to .csv in overleaf project: ./exp/xxx/data/???.csv
    }
}

int run2() {
    int blk_sz = 1600;
    std::string input = "./data/sin-1.csv";
    QUALIFY = 1000;
    for (int i = 1; i <= 3; i++) {
        double quantile = i * 0.10;
        int posv;

        posv = 1;// 1:iot.climate, sine, tpch
        //posv = 7;//, 8,9,10,... // gas, wind
        //posv = 0; QUALIFY =1;// timestamp
        
        long tmp = quantile_cntvf(quantile, input, blk_sz, 1);

        printf("%lf %ld \n", quantile, tmp);
    }
    return 0;
}

int run1() {
    int blk_sz = 1600;
    double b = 0.1;
    int round = 1; // 100 for sine
    int posv;
    for (int i = 12; i >= 9; i--) {
        // various selectivity
        //double selectivity = b*(i*1.0);
        double selectivity = 0.1;//*(i*1.0);

        // keep not changed
        tot_T = tot_F = tot_F2 = tot_F3 = 0;

        // change here for different dataset
        posv = 1;// 1:iot.climate, sine, tpch
        //posv = 7;//, 8,9,10,... // gas, wind
        //posv = 0; QUALIFY =1;// timestamp
        //QUALIFY = (1<<i);
        //int rle_bit_rest = 15;
        int delta_div = i;

        //printf("%lf\n", selectivity);
        test("./data/sin-1.csv", round, blk_sz, selectivity, delta_div);

        // use following line for throughput
        //tot_T = throughput(tot_T, round);
        //tot_F = throughput(tot_F, round);
        //tot_F2 = throughput(tot_F2, round);
        //tot_F3 = throughput(tot_F3, round);
        printf("%d %ld %ld %ld %ld\n", 14 - delta_div, tot_T, tot_F, tot_F2, tot_F3);
        //printf("%d %ld %ld %ld %ld\n", blk_sz, tot_T, tot_F, tot_F2, tot_F3);
        //blk_sz *= 5;
        //the output printed result > copy to .csv in overleaf project: ./exp/xxx/data/???.csv
    }
    printf("%ld\n", cntTup * 10);
    printf("%d %d\n", bit_delta, bit_rle);
    tot_T = throughput(tot_T, round);
    tot_F = throughput(tot_F, round);
    tot_F2 = throughput(tot_F2, round);
    tot_F3 = throughput(tot_F3, round);
    printf("%ld %ld %ld %ld\n", tot_T, tot_F, tot_F2, tot_F3);
    tot_T = tot_F = tot_F2 = tot_F3 = 0;
    //printf("Time: %ld ns - %ld ns\t Tup: %ld\n", tot_T, tot_F, cntTup);


    return 0;
}

int main() {
    //run1();
    //run2();
    //run3();
    run4();
}