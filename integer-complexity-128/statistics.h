namespace statistics{
	template<class T>
	double mean(const vector<T> &a){
		double s=0;
		for (auto &x:a)s+=x;
		return s/a.size();
	}
	template<class T>
	double variance(const vector<T> &a){
		double ave=mean(a),s=0;
		for (auto &x:a)s+=(x-ave)*(x-ave);
		return s/a.size();
	}
	template<class T>
	double stddev(const vector<T> &a){
		return sqrt(variance(a));
	}
} using namespace statistics;

