#include "encoder.hpp"



//

int test3(std::string prefix) {
    //assert();
    //std::string file = "./data/tpc-h.csv";
    //std::string encLocation = "./data/tpc/";
    encoding::Blocked::LSMBlockReader reader = encoding::Blocked::LSMBlockReader(prefix);
    reader.next8x32VectorDta();
    return 0;
}

int test2(std::string file, std::string encLocation, int col, const int QUALIFY, const int BLK_SZ) {
    //assert();
    //std::string file = "./data/tpc-h.csv";
    //std::string encLocation = "./data/tpc/";
    encoding::CSVEncoder encode = encoding::CSVEncoder(col, QUALIFY, BLK_SZ, encLocation);
    encode.encode(file.c_str());
    printf("Compress ratio: %lf\n", encode.finalize());
    return 0;
}

int test() {
        encoding::BitPack bp = encoding::BitPack();
        bp.setBitLen(6);
        int a[] = {1, 31, 23, 21, 12, 12, 12, 12, 13, 14, 14, 15, 16, 17, 18, 19, 22};
        long buff[100];
        memset(buff, 0, sizeof buff);
        int bias = 0;
        long pos = 0;
        //for(int j=0;j<4;j++) printf("%d ", buff[j]);
        printf("\n");
        for(int i=0;i<17;i++) {
            bp.encode_l4epi32(a[i], buff, bias, pos);
            for(int j=0;j<10;j++) printf("%ld ", buff[j]);
            printf("%ld %d \n", pos, bias);
        }
        char file[] = "./data.ts";
        encoding::Writer write = encoding::Writer(file);
        write.writeAsLongArray(buff, pos+1);
        memset(buff, 0, sizeof buff);

        encoding::Reader read = encoding::Reader(file);
        pos = read.readIntoLongArray(buff);
        printf("%ld\n", pos);

        int buffer[100];
        long posv = 0;

        bp.decode_l4epi32(buffer, buff, posv, pos);
        for(int i=0;i< posv;i++)printf("%d ", buffer[i]);
        printf("\n");
        return 0;
}

int main() {
    test3("./data/sin1/");
    return 0;
}