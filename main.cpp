/************************************************************
Author: Wangyi     Version: C++11  Date: 2020/11/26
Theme: convolutional code & viterbi decode
Channel: fading channel
Moduulaiton: BPSK
Note: Viterbi Decoder by Hamming Distance or Probability
Rate = 1/3 -> PunchPoints = 0;
Rate = 1/2 -> PunchPoints = 1;
***********************************************************/
#include "header.h"

int main(int argc, char *argv[]){
    
    if(argc < 3){
        cout<<"1st argument: beamforming algorithms"<<endl;
        cout<<"1: hard decoder based on hamming distance"<<endl;
        cout<<"2: soft decoder based on probability"<<endl;
        
        cout<<"2rd argument: PunchPoints"<<endl;
        cout<<"0: Rate = 1/3, PunchPoints = 0"<<endl;
        cout<<"1: Rate = 1/2, PunchPoints = 1"<<endl;
        
        return 0;
    }

    srand((unsigned)time(NULL));
    Initialize(argv);
    time_t start = time(NULL);

    for(int EbN0dB = MinEbN0dB; EbN0dB <= MaxEbN0dB; EbN0dB = EbN0dB + Step){

        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; ++loop){
            BitSource(Source);
            ConvEncoder(Source, Code); 
            Modulation(Code, Modu);
            ReceiverAWGN(Modu, RecAWGN);
            switch (*argv[1])
            {
                case '1':
                    ViterbiHardDecoder(Source, RecAWGN, Decode);
                    break;
                case '2':
                    ViterbiSoftDecoder(Source, RecAWGN, Decode);
                    break;
                
                default:
                    break;
            }

            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }

        #ifndef DebugMode
            BER = static_cast<double>(BER_TOTAL/(NLoop * LenBit));
            cout<<EbN0dB<<setw(20)<<BER<<endl;
            outfile<<EbN0dB<<setw(20)<<BER<<endl;
        #endif 

    }
    #ifndef DebugMode
        cout<<"time(s): "<<time(NULL) - start<<endl;
        outfile<<"time(s): "<<time(NULL) - start<<endl;
        outfile.close();
    #endif

    return 0;
}
