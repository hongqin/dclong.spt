
template<typename T0, typename T, typename... Args> double spt_one(arma::rowvec data, T0 t0, T t, size_t h, size_t n, Args... args){
    double t_null = t0(data, args...);
    size_t seq_count = 0;
    for(size_t i=1; i<n; ++i){
        seq_count += t(data, args...) >= t_null;
        if(seq_count >= h){
            return seq_count / static_cast<double>(i);
        }
    }
    ++ seq_count;
    return seq_count / static_cast<double>(n);
}

template<typename T0, typename T, typename... Args> std::vector<double> spt(arma::mat & data, T0 t0, T t, size_t h, size_t n, Args... args){
    size_t data_size = data.n_rows;
    std::vector<double> pvalues(data_size);
    for(size_t i=0; i<data_size; ++i){
        arma::rowvec data_i = data.row(i);
        pvalues[i] = spt_one(data_i, t0, t, h, n, args...);
    }
    return pvalues;
}


