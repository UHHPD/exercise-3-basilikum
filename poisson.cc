#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

double poisson(double mu, int k) {
    return pow(mu, k) * exp(-mu) / tgamma(k+1); // gamma function for factorial
}

int main() {
    using namespace std;

    vector<int> zaehler(11); // vector with 11 entries
    ifstream fin("datensumme.txt");
    ofstream fout("hist.txt");

    // calculate the histogram: 
    // index of the vector for the number; stored value as the frequency
    int N = 234;
    int n_i;
    for(int i = 0 ; i < N ; ++i) {
        fin >> n_i;
        zaehler[n_i] += 1;
    }    
    
    for(unsigned int k=0; k<zaehler.size(); k++){
        cout << k << ":" << zaehler[k] << endl;
        fout << k << " " << zaehler[k] << endl;
    }

    fin.close();
    fout.close();

    ofstream fout_comparison("histpoi.txt");
    for(unsigned int k=0; k<zaehler.size(); k++){
        fout_comparison << k << " " << zaehler[k] << " " << N*poisson(3.11538, k) << endl;
    }

    fout_comparison.close();

    return 0;
}