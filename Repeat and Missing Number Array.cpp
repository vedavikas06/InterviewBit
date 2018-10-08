vector<int> Solution::repeatedNumber(const vector<int> &A) {
    
    long long int sum =0,xr =0,t_xr=0,t_sum=0,n=A.size();
    
    for(long long int i=1;i<=n;i++){
        sum+=(long long int)i;
        t_sum+=(long long int)A[i-1];
        xr += (long long int)i*i;
        t_xr += (long long int)A[i-1]*A[i-1];
        
    }
    long long int ac_pr = sum-t_sum;
    
    long long int xd = xr - t_xr;
    
    long long int xa = xd/ac_pr;
    
    long long int a = (xa+ac_pr)/2,b = (xa-ac_pr)/2;
    
    vector<int> v;
    
    v.push_back((int)b);
    v.push_back((int)a);
    
    return v;
}
