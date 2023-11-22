#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


double poisson(double mu, int k) {
    return pow(mu, k) * exp(-mu) / tgamma(k+1);
}

double chi2(double x, double dof){
    return pow(x, dof/2 - 1) * exp(-x/2) / (pow(2, dof/2) * tgamma(dof/2 + 1));
}

// evaluate Poisson that takes vector
double prob(std::vector<int> daten, double mu){
    double res = 1;
    for(int k : daten){
        res *= poisson(mu, k);
    }

    return res;
}

int main() {
    using namespace std;

    int N = 234;
    // double mu = 
    vector<int> daten;

    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < N ; ++i) {
        fin >> n_i;
        // cout << n_i;
        daten.push_back(n_i);
    }

    double mu_exact = 3.11538;
    cout << "joint probability (single Poisson) = " << prob(daten, mu_exact) << endl;
    double nll_min = -2 * log(prob(daten, mu_exact));

    // Scan μ with a step size of 0.1 from 0~6 --> storing into "likelihood.txt"
    ofstream fout_mu_scan("likelihood.txt");
    ofstream fout_mu_scan_nll("nll.txt");
    ofstream fout_mu_scan_deltanll("deltanll.txt");

    double mu_scan_start = 0;
    double mu_scan_end = 6;
    double mu_scan_step = 0.1;
    int mu_scan_N = int((mu_scan_end - mu_scan_start)/mu_scan_step);

    for(int i=0; i<mu_scan_N; i++){
        double mu_scan = mu_scan_start + i * mu_scan_step;
        fout_mu_scan << mu_scan << " " << prob(daten, mu_scan) << endl;
        fout_mu_scan_nll << mu_scan << " " << -2 * log(prob(daten, mu_scan)) << endl;
        fout_mu_scan_deltanll << mu_scan << " " << -2 * log(prob(daten, mu_scan)) - nll_min << endl;
    }
    
    fout_mu_scan.close();
    fout_mu_scan_nll.close();
    fout_mu_scan_deltanll.close();
    fin.close();

    // Wilks' theorem --> chi2 distribution
    double probDenominator = 1;
    for(int k : daten){
        probDenominator *= poisson(k, k);
    }
    cout << "joint probability (234 Poissons) = " << probDenominator << endl;
    // note that the mean value is not the most probable in Poisson, this is probably just for convenience

    double neg_log_likelihood_ratio = -2 * (log(prob(daten, mu_exact)) - log(probDenominator));
    
    cout << "negative log-likelihood ratio = " << neg_log_likelihood_ratio << endl; // this should be close to the degrees of freedom 233!!!
    int n_dof = N - 1;
    // cout << "chi2(DOF=234-1, neg_log_likelihood_ratio) at = " << chi2(neg_log_likelihood_ratio, 233) << endl;
    
    // transform to a variable that obeys standard normal distribution
    // value far away from 0 would indicate the Null hypothesis is likely to be wrong
    double z = (neg_log_likelihood_ratio - n_dof) / sqrt(2 * n_dof);
    
    // cout << "chi2(DOF=234-1, neg_log_likelihood_ratio) at = " << chi2(neg_log_likelihood_ratio, 233) << endl;

    // fron Wilks' theorem, the neg_log_likelihood_ratio should follow the chi2 distribution
    // probability of actually getting this value from chi2 distribution?
    cout << "z = " << z << endl;

    
    return 0;
}
