vector<int> Solution::solve(vector<int> &A, vector<int> &B) {
    
    
    int n = A.size(),m=B.size();
    vector<pair<int,int> > v;
    v.push_back({INT_MIN,2});
    
    for(int i=0;i<n;i++){
        v.push_back({A[i],0});
    }
    for(int i=0;i<m;i++){
        v.push_back({B[i],1});
    }
    
    sort(v.begin(),v.end());
    
    int diff = INT_MIN,fir,sec;
    
    int cnt1=0,cnt2=0;
    
    for(int i=v.size()-1;i>=0;i--){
        pair<int,int> p =v[i];
        
        int tot = 3*cnt1 + (2*(n-cnt1)) - (3*cnt2 + (2*(m-cnt2)));
        
        if(tot>diff){
            fir = 3*cnt1 + (2*(n-cnt1));
            sec = (3*cnt2 + (2*(m-cnt2)));
            diff =tot;
        }
         if(p.second ==0){
            cnt1++;
        }else{
            cnt2++;
        }
        
    }
    
    
   
    vector<int> v1;
    
    v1.push_back(fir);v1.push_back(sec);
    
    return v1;
    
}
