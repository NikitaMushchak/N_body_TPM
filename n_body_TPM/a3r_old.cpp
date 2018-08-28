#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>


namespace ai{
    template<typename T>
    std::string string(const T value){
        std::ostringstream stream;

        stream << value;

        return stream.str();
    }

    template<typename T>
    inline bool saveA3R(
        const std::string filename,
        std::vector< std::vector<T> > positions,
        const double radius = 1
    ){
        const std::size_t numberOfParticles = positions.size();

        const int integerNumberOfParticles = (int) numberOfParticles;

        std::ofstream a3r(filename, std::ios::binary | std::ios::out);

        const int startByte = 36;
        const int futureValue = 0;

        const int intSize = 4;
        const int floatSize = 4;
        const int doubleSize = 8;

        if(sizeof(int) != intSize){
            throw std::runtime_error(
                ai::string("int size is ") + ai::string(sizeof(int))
                + ai::string(" instead of ") + ai::string(intSize)
            );
        }
        if(sizeof(float) != floatSize){
            throw std::runtime_error(
                ai::string("float size is ") + ai::string(sizeof(float))
                + ai::string(" instead of ") + ai::string(floatSize)
            );
        }
        if(sizeof(double) != doubleSize){
            throw std::runtime_error(
                ai::string("double size is ") + ai::string(sizeof(double))
                + ai::string(" instead of ") + ai::string(doubleSize)
            );
        }

        a3r.write("a3r\0", 4);
        a3r.write((char*) &integerNumberOfParticles, intSize);
        a3r.write((char*) &startByte, intSize);
        a3r.write("AiLib 1.1.0\0", 12);
        a3r.write((char*) &radius, doubleSize);
        a3r.write((char*) &futureValue, intSize);

        for(std::size_t i = 0; i < numberOfParticles; ++i){
            float x = (float) positions[i][0];
            float y = (float) positions[i][1];
            float z = (float) positions[i][2];

            //float x = (float) positions[i].x;
            //float y = (float) positions[i].y;
            //float z = (float) positions[i].z;

            a3r.write((char*) &x, floatSize);
            a3r.write((char*) &y, floatSize);
            a3r.write((char*) &z, floatSize);
        }

        a3r.close();

        return true;
    }
}

//int main(){
//    std::vector< std::vector<double> > particles;
//
//    for(size_t i = 0; i < 5; ++i){
//        particles.push_back(std::vector<double>{(double) i, 5. - i, 1.7 * i});
//    }
//
//    ai::saveA3R("./test.a3r", particles);
//
//    return 0;
//}
