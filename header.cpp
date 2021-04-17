/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: header.cpp
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                          /* total number of error bits */
double BER = 0.0;                          /* total number of error symbols */
fstream outfile;

SourceMatrix Source;                                /* source codewords */
CodeMatrix Code;                                    /* codewords after coding */
ModuMatrix Modu;                                    /* symbols after modulation */
SourceMatrix Decode;
RecAWGNMatrix RecAWGN;


void Initialize(char* argv[]){

    cout<<"Convolution_Code_with_Viterbi_Decoder"<<endl;
    outfile.open("viterbi.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"Convolution_Code_with_Viterbi_Decoder"<<endl;
    
    cout<<"AWGN"<<endl;
    outfile<<"AWGN"<<endl;
    
    switch (*argv[1]){
        case '1':
            cout<<"hard decoder based on hamming distance"<<endl;
            outfile<<"hard decoder based on hamming distance"<<endl;
            break;
        
        case '2':
            cout<<"soft decoder based on probability"<<endl;
            outfile<<"soft decoder based on probability"<<endl;
            break;
        
        default:
            break;
    }
    cout<<"Length of coded bits: "<<NJ<<setw(15)<<"Rate: "<<ReciRate<<setw(15)<<"Punched Points: "<<PunchPoints<<setw(15)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Length of coded bits: "<<NJ<<setw(15)<<"Rate: "<<ReciRate<<setw(15)<<"Punched Points: "<<PunchPoints<<setw(15)<<"NLoop: "<<NLoop<<endl;   
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;
}

/* key point */
void ChannelInitialize(int ebN0dB){

    double EbN0 = pow(10, (double)ebN0dB/10); 
    double EsN0 = EbN0/static_cast<double>(ReciRate - PunchPoints);
    double N0 = 1.0 / EsN0;
    N_Var = N0 / 2;
    BER_TOTAL = 0;               
    BER = 0;  
#ifdef DebugMode
    cout<<"N_var: "<<N_Var<<endl;
#endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int i = 0; i < LenBit; i++){
        source(i)= rand()%2;
    }
#ifdef DebugMode
    cout<<"source: "<<endl<<source<<endl;
#endif
}


void ConvEncoder(SourceMatrix& source, CodeMatrix& code){

    /* 5 7 7 */
    int reg1=0,reg2=0;
    for(int i = 0; i < LenBit; i++){
        code(ReciRate*i + 0) = source(i) xor reg2; 
        code(ReciRate*i + 1) = source(i) xor reg1 xor reg2;
        code(ReciRate*i + 2) = source(i) xor reg1 xor reg2;

        reg2 = reg1;
        reg1 = source(i);
    }
#ifdef DebugMode
    cout<<"code: "<<endl<<code<<endl;
#endif
}


void Modulation(CodeMatrix& code, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int nj = 0; nj < Nj; nj++){
        for(int j = 0; j < J; j++){
            modu(nj, j) = 2 * code(nj * J + j) - 1;
        }
    }
#ifdef DebugMode
    cout<<"bpsk signal: "<<endl<<modu<<endl;
#endif
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void ReceiverAWGN(ModuMatrix& modu, RecAWGNMatrix& recAWGN){

    for(int nj = 0; nj < Nj; nj++){
        for(int j = 0; j < J; j++) {
            ComplexD tmp = AWGN(N_Var);
            recAWGN(nj * J + j) = modu(nj, j) + tmp;
        }
    }

#ifdef DebugMode
    cout<<"recAWGN: "<<endl;
    for(size_t index = 0; index < ReciRate; index++){
        for(size_t index2 = 0; index2 < LenBit; index2++){
            cout<<recAWGN(index * LenBit + index2)<<" ";
        }
        cout<<endl;
    }
#endif
}

ComplexD operator-(ComplexD comp, int b){
    return ComplexD(comp.real()-b, comp.imag());
}