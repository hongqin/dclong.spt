#include<algorithm>

template<typename InputIt> void shuffle_left(InputIt first, InputIt last, int k) {
    int n = last - first;
    int n_minus_one = n - 1;
    int upperBound = k;
	//no need to shuffle if only 1 element is left
	if (upperBound>=n) {
		upperBound = n_minus_one;
	}
	for (int i=0; i<upperBound; ++i) {
		//generate an index between i and n - 1 (inclusive)
		int index = static_cast<int>(Rf_runif(i, n));
        index = std::min(index, n);
        index = std::max(index, i-1);
        // size_t index = k;
		//swap x[i] and x[index]
        std::swap(*(first + i), *(first + index));
	}
}

double abs_mean_diff0(arma::rowvec & data, std::size_t n1){
    std::size_t n = data.size();
    return std::abs(std::accumulate(data.begin(), data.begin()+n1, 0.0)/n1 - std::accumulate(data.begin()+n1, data.end(), 0.0)/(n-n1));
}

double abs_mean_diff(arma::rowvec & data, std::size_t n1){
    // permute data
    shuffle_left(data.begin(), data.end(), n1);
    // recalculating test statistic
    return abs_mean_diff0(data, n1);
}


double max_abs_corr0(arma::rowvec & y, const arma::mat & X){
    arma::rowvec corrs = y * X;
    for(arma::rowvec::iterator it=corrs.begin(); it!=corrs.end(); ++it){
        *it = std::abs(*it);
    }
    return *std::max_element(corrs.begin(), corrs.end());
}


double max_abs_corr(arma::rowvec & y, const arma::mat & X){
    // permute data
    shuffle_left(y.begin(), y.end(), y.size());
    // recalculating test statistic
    return max_abs_corr0(y, X);
}

