int Solution::maxSpecialProduct(vector<int> &A) {
    
    int n = A.size();
    
    
    vector<long long> l(n,0),r(n,0);
    stack<long long> s;
    for(long long i= 0;i<n;i++){
        
        if(s.empty()){
            s.push(i);
        }else{
            
            while(!s.empty() && A[s.top()] <= A[i]){
                s.pop();
            }
            
            if(!s.empty()){
                l[i] = s.top();
                
            }
            
            s.push(i);
            

        }
    }
    while(!s.empty()){
    s.pop();
    }
    
     for(long long i= n-1;i>=0;i--){
        
        if(s.empty()){
            s.push(i);
        }else{
            
            while(!s.empty() && A[s.top()] <= A[i]){
                s.pop();
            }
            
            if(!s.empty()){
                r[i] = s.top();
                
            }
            s.push(i);

        }
        
    }
    
    long long maxx = 0;
    long long mod = 1000000007;
    for(int i =0;i<n;i++){
        //cout << l[i] << " " << r[i] << endl;
        long long val = (l[i]*r[i]*1LL); 
        maxx= max(val,maxx);
        
        
        
    }
    
    
    return maxx%mod;
    
    
    
    
    
}
