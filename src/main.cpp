#include "readData.h"
#include "CustoIn.h"
#include "infoSeq.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <random>
#include <sys/timeb.h>
#include <sys/resource.h>

#include <iomanip>
#include <string>

using namespace std;

double ** mJobs; // matriz de adjacencia
double ** mSetupTimes; // matriz reorganizada;
int n; // quantidade total de vertices
vector<vector<int>> arrangedMatrix;
unsigned seed;

const double EulerConstant = std::exp(1.0);

void printData(int n, double ** mJobs, double ** mSetupTimes);
vector<int> construction(int n, double ** mJobs, double ** mSetupTimes, double alfa, double &cost);
double sequenceTime(vector<int> s, double ** mJobs, double **mSetupTimes);
vector<CustoIn> calculaCusto(vector<int> listaCandidatos, vector<int> &s, double custoAtual);
bool comp(const CustoIn& a, const CustoIn& b);
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes);
void arrangeMatrix(int dimension, double **adjMatrix, vector<vector<int>> &arrangedMatrix);
void removeSelected(vector<int> &s, vector<int> &absentJobs, int l_t, int j_t, int index_j_t ,double alfa, double beta);
int genRandomInteger(int min, int max);
void ruin(vector<int> &s, vector<int> &absentJobs, int l_s_max, double alfa, double beta);
void recreate(vector<int> &s, double &mkspan, vector<int> &absentJobs);
double costInsertionJob(vector<int> &s, int job, int pos);
vector<int> localSearch(vector<int> &s, double &mkspan, int nIter, int t_init, int l_s_max, double alfa, double beta);
vector<int> sisr_ruin_recreate(vector<int> s, double &mkspan ,int l_s_max, int alfa, int beta);
double realcostInsertionJob(vector<int> &s, int job, int pos);

int main(int argc, char** argv) {

    vector<int> s;
    vector<int> absentJobs;
    double cost;
    readData(argc, argv, &n, &mJobs, &mSetupTimes);
    arrangeMatrix(n, mSetupTimes, arrangedMatrix);

    sequencesMatrix = new infoSequence*[n+1];
    for(int i = 0; i <= n; ++i){
        sequencesMatrix[i] = new infoSequence[n+1];
    }

    seed = time(0);
    //unsigned seed = 1611170583;
    //seed = 1611495465; 
    //seed = 1611497795;
    //seed = 1611500157;
    //seed =  1611503946; // seg fault
    //seed = 1611504278; // seg fault again
    //seed = 1611778693;
    cout << "\nseed: " << seed << endl;
    srand(seed);

    double makespan;
    s = construction(n, mJobs, mSetupTimes, 0.1, makespan);
    printSolution(s, mJobs, mSetupTimes);
    //ruin(s, absentJobs, 10, 0.5, 0.1);
    //printSolution(s, mJobs, mSetupTimes);
    //recreate(s, absentJobs);
    //printSolution(s, mJobs, mSetupTimes);
    
    s = localSearch(s, makespan ,1000000, 100, 5, 1, 0.01);
    printSolution(s, mJobs, mSetupTimes);
    cout << "custo = " << makespan << endl;

    /*s = {1,2,3,4,5,6};
    s.insert(s.begin()+6, 7);
    printSolution(s, mJobs, mSetupTimes);*/


    return 0;
}

 void printData(int n, double ** mJobs, double ** mSetupTimes){
    cout <<"Jobs QTD: "<< n << "\n";

    for (int i = 1; i <= n; i++){
        for(int j = 0; j < 4; j++){
            cout << mJobs[i][j] << " ";
        }
    }
    
   for (int i = 0; i <= n; i++){
        for(int j = 1; j <= n; j++){
            cout << mSetupTimes[i][j] << " ";
        }
    }
 } 

vector<int> construction(int n, double ** mJobs, double ** mSetupTimes, double alfa, double &cost){
    vector<int> s;
    vector<int> listaCandidatos;
    double custoAtual = 0;

    for(int i = 1; i <= n; i++){
        listaCandidatos.push_back(i); // insere todos os nós na lista de candidatos
    }

    int k; // indice do job a ser adicionado

    for(int j = 0; j < 1; j++){ // tamanho subsequencia inicial
        if(alfa == 0)
            k = 0; // escolhe o primeiro job caso alfa = 0
        else
            int k = rand() % ((int)std::floor(alfa *listaCandidatos.size())); // seleciona um job dentre os alfa menores release dates
        s.push_back(listaCandidatos[k]); // adiciona o job a solução
        listaCandidatos.erase(listaCandidatos.begin() + k); // apaga da lista de candidatos oq ja foi pra solução
        if(j == 0){
            custoAtual += mJobs[s[0]][1] + mSetupTimes[0][s[j]] + mJobs[s[0]][2]; // computa o tempo de conclusão do primeiro job
        }
        else{ // tempo de conclusão do resto da sequencia
            if(custoAtual >= mJobs[s[j]][1]){
                custoAtual += mSetupTimes[s[j-1]][s[j]] + mJobs[s[j]][2];
            }
            else{
                custoAtual = mSetupTimes[s[j-1]][s[j]] + mJobs[s[j]][2] + mJobs[s[j]][1];
            }
        }
    }
   // cout << custoAtual << endl;

    vector<CustoIn> custoInsertion = calculaCusto(listaCandidatos, s, custoAtual); // calcula os custo de inserção dos candidatos
    std::sort(custoInsertion.begin(), custoInsertion.end(), comp); // ordena de forma crescente de acordo com os custos
    int sel;
    while(!listaCandidatos.empty()){
        if(alfa == 0){
            sel = 0;
        }
        else{
            sel = rand() % ((int)std::floor(alfa * (custoInsertion.size() - 1)) + 1); // escolhe um nó dentro de uma faixa definida por alfa
        }

        s.push_back(custoInsertion[sel].noIn); // insere nó na sequencia
        custoAtual += max(custoInsertion[sel].custo,0.0) + mSetupTimes[s[s.size()-2]][s[s.size()-1]] + mJobs[s[s.size()-1]][2]; // atualiza tempo de conclusão da sequencia
        
        for(int i = 0; i < listaCandidatos.size(); i++){
            if(listaCandidatos[i]==custoInsertion[sel].noIn)
                listaCandidatos.erase(listaCandidatos.begin() + i); // exclui o nó da lista de candidatos
            }
        
        custoInsertion.erase(custoInsertion.begin(), custoInsertion.end()); // exclui o nó da lista de custos
        custoInsertion = calculaCusto(listaCandidatos, s, custoAtual); // calcula os novos custos de inserção
        std::sort(custoInsertion.begin(), custoInsertion.end(), comp); // ordena os custos
        
    }
    
    setSequencesMatrix(sequencesMatrix,s,n,mJobs,mSetupTimes); // construção da matriz de subsequencias
    
    infoSequence initialSolution;
    initialSolution = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][n-1]); // concatena a sequencia com tarefa dummy para computar initial setup
    //idleTime = initialSolution.waitingTime;
    cost = initialSolution.initialTime + initialSolution.duration; // custo da solução inicial
    //printSolution(s, mJobs, mSetupTimes);
    //cost = sequenceTime(s, mJobs, mSetupTimes);

    return s; // retorna a sequencia
}

vector<CustoIn> calculaCusto(vector<int> listaCandidatos, vector<int> &s, double custoAtual){
  vector<CustoIn> custoInsertion (listaCandidatos.size());
    //for(int i = 0, l = 0; i < s.size()-1; i++){
        int l = 0;
        for(auto k : listaCandidatos){
            custoInsertion[l].custo = max((mJobs[k][1] - custoAtual), 0.0) + mSetupTimes[s[s.size() - 1]][k]; // custo de inserção definido pelo Idle Time gerado
            custoInsertion[l].noIn = k; // nó inserido
            //custoInsertion[l].arestaOut = i; // posicao de inserção;
            l++;
        }

  return custoInsertion;

}

inline double sequenceTime(vector<int> s, double ** mJobs, double **mSetupTimes){ // calcula o custo percorrendo a sequencia
    
    double cTime = mSetupTimes[0][s[0]] + mJobs[s[0]][2] + mJobs[s[0]][1]; // tempo de conclusão da primeira tarefa
    double totalWT = mJobs[s[0]][1]; // idle time inicial
    /*for(int i = 0; i < s.size(); i++){
        cout << s[i] << " ";
    }*/
    for(int i = 0, j = 1; j < s.size(); i++ ,j++){
        //cout<< i << " " << j << endl;
        if(cTime >= mJobs[s[j]][1]){
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2];
        }
        else{
            totalWT += mJobs[s[j]][1] - cTime;
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2] + (mJobs[s[j]][1] - cTime);
        }
    }
    //cout << "\ntotalWT: " << totalWT - mJobs[s[0]][1]<< endl;
    return cTime;
} 


 bool comp(const CustoIn& a, const CustoIn& b) // comparação dos custos utilizada para ordenar os objetos
  {
    return a.custo < b.custo;
  }

void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes){
    for(int i = 0; i < solucao.size(); i++){
        cout << solucao[i] << " ";
    }
    cout << "\n";
    /*for(int i = 0; i < n; i++){
        cout << mJobs[solucao[i]][1] << " ";
    }*/

    double cTime = mSetupTimes[0][solucao[0]] + mJobs[solucao[0]][2] + mJobs[solucao[0]][1];
    double totalWT = mJobs[solucao[0]][1];
    /*for(int i = 0; i < s.size(); i++){
        cout << s[i] << " ";
    }*/
    //cout << totalWT;
    for(int i = 0, j = 1; j < solucao.size(); i++ ,j++){
        //cout<< i << " " << j << endl;
        if(cTime >= mJobs[solucao[j]][1]){
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2];
            //cout << " " << 0;
        }
        else{
            totalWT = mJobs[solucao[j]][1] - cTime;
            //cout << " " << totalWT;
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2] + (mJobs[solucao[j]][1] - cTime);
        }
    }
}

void arrangeMatrix(int dimension, double **adjMatrix, vector<vector<int>> &arrangedMatrix)
{
  struct arrange
  {
    int counter;
    double **adjMatrix;
    bool operator()(const int &a, const int &b) const
    {
      return adjMatrix[counter][a] < adjMatrix[counter][b];
    }
  } arranged;

 
  arranged.counter = 0;
  arranged.adjMatrix = adjMatrix;

  /*for (int i = 1; i <= dimension; i++)
  {
    optimal.push_back(i);
  }*/

    for (int i = 0; i <= dimension; i++){
        vector<int> optimal;
        for (int j = 1; j <= dimension; j++){
            if(i != j)
                optimal.push_back(j);
        }      
        arrangedMatrix.push_back(optimal);
        sort(arrangedMatrix[i].begin(), arrangedMatrix[i].end(), arranged);
        arranged.counter++;
    }

}

void ruin(vector<int> &s, vector<int> &absentJobs, int l_s_max, double alfa, double beta){
   // int l_s_max = L_max;
    int seedJob = genRandomInteger(1,n);

    //for(int i = 0; i < arrangedMatrix[seedJob].size(); i++){
     //   int j = arrangedMatrix[seedJob][i];
     //   if(absentJobs.empty() || (j != seedJob && absentJobs.end() != find(absentJobs.begin(), absentJobs.end(), j))){
           int j_t = seedJob;

           auto it = find(s.begin(), s.end(), j_t);
           int index_j_t = it - s.begin();

           int l_t = genRandomInteger(1, min((int)s.size()-index_j_t, l_s_max));
           removeSelected(s, absentJobs, l_t, j_t, index_j_t, alfa, beta);
           setSequencesMatrix(sequencesMatrix, s, s.size(), mJobs, mSetupTimes);
       // }

    //}
}

void removeSelected(vector<int> &s, vector<int> &absentJobs, int l_t, int j_t, int index_j_t, double alfa, double beta){

    int stringEnd;
    int stringBegin;
    //int m;
    //double beta = 2;
    //cout << j_t << endl;

    //random_device rd;
    mt19937 gen(seed);
    bernoulli_distribution d(alfa);


    if(d(gen) == true){ // split string removal
        int m = 1; // qtd jobs preservados

        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0,1.0);

       // if(distribution(generator) > beta){

        while(distribution(generator) < beta && m < s.size() - l_t){
            m++;
        }
        

        do{
            stringBegin = genRandomInteger(max(0,index_j_t-l_t), index_j_t); //indice do inicio da string na soluçao
            //cout << "begin: " << stringBegin << endl;
            stringEnd = stringBegin + l_t; // indice do fim da string
            //cout << "end: " << stringEnd << endl;
        }while(stringEnd >= s.size());

        int preservedStringBegin = genRandomInteger(stringBegin, stringEnd - m);

       // cout <<"m =  " << m << " begin = " << stringBegin << " end = " << stringEnd << " mBegin = " << preservedStringBegin << endl;

        if(preservedStringBegin == stringBegin){
            absentJobs.insert(absentJobs.end(), s.begin() + preservedStringBegin + m, s.begin() + stringEnd + 1);
            s.erase(s.begin() + preservedStringBegin + m, s.begin() + stringEnd + 1);
        }
        if(preservedStringBegin - 1 == stringEnd - m){
            absentJobs.insert(absentJobs.end(), s.begin() + stringBegin, s.begin() + stringEnd - m + 1);
            s.erase(s.begin() + stringBegin, s.begin() + stringEnd - m + 1);
        }
        if(preservedStringBegin > stringBegin && preservedStringBegin <= stringEnd - m){
            absentJobs.insert(absentJobs.end(), s.begin() + stringBegin, s.begin() + preservedStringBegin);
            absentJobs.insert(absentJobs.end(), s.begin() + preservedStringBegin + m, s.begin() + stringEnd + 1);
            s.erase(s.begin() + stringBegin, s.begin() + preservedStringBegin);
            s.erase(s.begin() + stringBegin + m, s.begin() + stringEnd - preservedStringBegin + stringBegin + 1);

        }

        /*cout << "absent = ";
        for (int i = 0; i < absentJobs.size(); i++){
            cout << absentJobs[i] << " ";
        }
        cout << endl;*/


    }
    else{ //string removal

        do{
            stringBegin = genRandomInteger(max(0,index_j_t-l_t), index_j_t);
            //cout << "begin: " << stringBegin << endl;
            stringEnd = stringBegin + l_t;
            //cout << "end: " << stringEnd << endl;
        }while(stringEnd >= s.size());

        //cout <<" begin = " << stringBegin << " end = " << stringEnd << endl;


        absentJobs.insert(absentJobs.end(), s.begin() + stringBegin, s.begin() + stringEnd + 1);
        s.erase(s.begin() + stringBegin, s.begin() + stringEnd + 1);

    
    }


    

}

int genRandomInteger(int min, int max){
   return min + (rand() % (max - min + 1));
}


void recreate(vector<int> &s, double &mkspan, vector<int> &absentJobs){

    int blinkRate = 0;
    
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0,1.0);

    //cout <<  distribution(generator) << endl;
    
    sort(absentJobs.begin(), absentJobs.end());
    
    for (auto j : absentJobs){
        int best_pos = -1;

        if(s.size() == 0)
            best_pos = 0;

        double costInsertionInPos;
        double costInsertionBestPos;
        
        //printSolution(s, mJobs, mSetupTimes);
        for(int k = 0; k <= s.size(); k++){
            //cout << "k = " << k << endl;
            //cout << "j = " << j << endl;
            if(distribution(generator) < 1 - blinkRate){
               /* double calc = costInsertionJob(s, j, k);
                double real = realcostInsertionJob(s, j, k);
                    if(calc != real)
                        cout << "calc: " << calc << " real: " << real <<" pos =  " << k << " j = " << j << " size = " << s.size() << endl;*/
                costInsertionInPos = costInsertionJob(s, j, k);
                if(best_pos == -1 || costInsertionInPos < costInsertionBestPos){   
                    best_pos = k;
                    costInsertionBestPos = costInsertionInPos;
                }
            } 
        }
        //cout << s.size() << endl;
        s.insert(s.begin() + best_pos, j);
        //cout << "best_pos: "<< best_pos << endl;
        absentJobs.pop_back();
        mkspan = costInsertionBestPos;

        //atualização da matriz de subsequencias

        /*infoSequence jobSeq; // subsequencia unitária com o job a ser inserido;

        jobSeq.initialTime = mJobs[j][1];
        jobSeq.duration = mJobs[j][2];
        jobSeq.firstJob = j;
        jobSeq.lastJob = j;*/

        /* for(int i = 0; i < n; i++){
            sequencesMatrix[i][i].firstJob = s[i];
            sequencesMatrix[i][i].lastJob = s[i];
            sequencesMatrix[i][i].duration = mJobs[s[i]][2];
            sequencesMatrix[i][i].initialTime = mJobs[s[i]][1];
        }

        for(int l = 0; l < best_pos; l++){
            for(int m = best_pos; m < s.size(); m++){
                sequencesMatrix[l][m] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[l][m-1], sequencesMatrix[m][m]);
                sequencesMatrix[m][l] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[m][m], sequencesMatrix[l][m-1]);
            }
            
        }*/

        setSequencesMatrix(sequencesMatrix, s, s.size(), mJobs, mSetupTimes);
    }

}

double realcostInsertionJob(vector<int> &s, int job, int pos){
    vector<int> s_1 = s;
    s_1.insert(s_1.begin() + pos, job);

    //printSolution(s_1, mJobs, mSetupTimes);

    return sequenceTime(s_1, mJobs, mSetupTimes);
}

double costInsertionJob(vector<int> &s, int job, int pos){
    //vector<int> s_1 = s;
    //s_1.insert(s_1.begin() + pos, job);

    //printSolution(s_1, mJobs, mSetupTimes);

    //return sequenceTime(s_1, mJobs, mSetupTimes);
    
    infoSequence dummyJob;
    infoSequence seq;
    infoSequence jobSeq; // subsequencia unitária com o job a ser inserido;

    jobSeq.initialTime = mJobs[job][1];
    jobSeq.duration = mJobs[job][2];
    jobSeq.firstJob = job;
    jobSeq.lastJob = job;

    dummyJob.initialTime = 0;
    dummyJob.duration = 0;
    dummyJob.firstJob = 0;
    dummyJob.lastJob = 0;


   /* if(pos != 0)
        seq = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[0][pos-1]);
    else
        seq = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, jobSeq);

    if(pos != 0 && pos != s.size() - 1)
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, jobSeq);

    if(pos != s.size() - 1)
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, sequencesMatrix[pos][s.size()]);
    else
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, jobSeq);
    */


    if(pos != 0  && pos != s.size()){
        seq = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[0][pos-1]);
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, jobSeq);
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, sequencesMatrix[pos][s.size()-1]);
   }

    if(pos == 0){
        seq = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, jobSeq);
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, sequencesMatrix[pos][s.size()-1]);
    }

    if(pos == s.size()){
        seq = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[0][pos-1]);
        seq = concatenateSequencev2(mSetupTimes, mJobs, seq, jobSeq);
    }


    return seq.initialTime + seq.duration;

}

vector<int> localSearch(vector<int> &s, double &mkspan, int nIter, int t_init, int l_s_max, double alfa, double beta){
    
    vector<int> bestSolution = s;
    double bestCost = mkspan;
    vector<int> iterationSolution;
    double iterationCost;
    //vector<int> absentJobs;
    int t = t_init;

    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0,1.0);
    
    for(int i = 0; i < nIter; i++){
        vector<int> absentJobs;
        //cout << "iter "<< i << endl;
        iterationSolution = sisr_ruin_recreate(s, iterationCost, l_s_max, alfa, beta);
        if(iterationCost < mkspan - t*log(distribution(generator))){
            s = iterationSolution;
            mkspan = iterationCost;
        }
        if(iterationCost < bestCost){
            bestSolution = iterationSolution;
            bestCost = iterationCost;
        }
        //iterationSolution = s;
        
        t = t*0.01;
    }

    return bestSolution;

}

vector<int> sisr_ruin_recreate(vector<int> s, double &mkspan ,int l_s_max, int alfa, int beta){
    vector<int> absentJobs;

    ruin(s, absentJobs, l_s_max, alfa, beta);
    recreate(s, mkspan, absentJobs);


    return s;
}