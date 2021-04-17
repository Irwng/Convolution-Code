#ifndef CONV_VITERBI_HEADER_H
#define CONV_VITERBI_HEADER_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <complex>
#include <Eigen/SVD>
#include <Eigen/Dense>

#define randN() (rand()/(double)RAND_MAX)     /* Random value in [0,30000] */

using namespace std;
using namespace Eigen;

/* basic model parameters */ 
constexpr int U = 6;                              /* Number of users */
constexpr int K = 50;                             /* Length of per user's bit stream,length of vector(b)=K*U */
constexpr int LenBit = U * K;                     /* number of bits of all users */
constexpr int ReciRate = 3;                       /* reciprocal of the rate */
constexpr int PunchPoints = 0;                    /* points bunched */
constexpr int NODE = 4;
constexpr int Route = NODE;                       /* number of the routes to keep in viterbi decoder */
constexpr int NJ = LenBit * ReciRate;             /* number of coded bits, J*Nj = K*U*ReciRate,one block per transmissi  on  */
constexpr int J = 6;                              /* bits per block, length of vector(c) = J*Nj */
constexpr int Nj = NJ/J;                          /* number of blocks,J*Nj = K*U, one block per transmission  */

constexpr int MinEbN0dB = 0;
#ifdef DebugMode
    constexpr long NLoop = pow(10, 0);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = MinEbN0dB;
#else
    constexpr long NLoop = pow(10, 5);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = 10;           
#endif//DebugMode
constexpr int Step = 1;    

typedef complex<double> ComplexD;

/* source codewords */
typedef Matrix<int, 1, LenBit> SourceMatrix;
extern SourceMatrix Source;

/* codewords after coding */
typedef Matrix<int, 1, NJ> CodeMatrix;
extern CodeMatrix Code;

/* symbols after modulation */
typedef Matrix<ComplexD, Nj, J> ModuMatrix;
extern ModuMatrix Modu;

/* codewords with AWGN noise */
typedef Matrix<ComplexD, 1, NJ> RecAWGNMatrix;
extern RecAWGNMatrix RecAWGN;

/* source codewords */
extern SourceMatrix Decode;

extern double N_Var;                        /* variance of Noise*/
extern double BER_TOTAL;                    /* total number of error bits*/
extern double BER;                          /* error bits rate */
constexpr double PI = 3.141592653589793;

extern fstream outfile;

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void Initialize(char *argv[]);

/* functions declartions */
ComplexD AWGN(double nvar);
ComplexD operator-(ComplexD comp, int b);

void ChannelInitialize(int ebN0dB);
void BitSource(SourceMatrix& source);
void ConvEncoder(SourceMatrix& source, CodeMatrix& code);
void Modulation(CodeMatrix& code, ModuMatrix& modu);
void ReceiverAWGN(ModuMatrix& modu, RecAWGNMatrix& recAWGN);

void ViterbiHardDecoder(SourceMatrix& source, RecAWGNMatrix& recAWGN, SourceMatrix& decode);
void ViterbiSoftDecoder(SourceMatrix& source, RecAWGNMatrix& recAWGN, SourceMatrix& decode);

#endif //CONV_VITERBI_HEADER_H
