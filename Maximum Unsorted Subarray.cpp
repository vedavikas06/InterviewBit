vector<int> Solution::subUnsort(vector<int> &A) {
    
    int n = A.size();
    
    
    int st = 0,en = n-1;
    
    for(int i=1;i<n;i++){
        if(A[i]<A[i-1]){
            st=i-1;
            break;
        }
    }
    
    for(int i=n-2;i>=0;i--){
        if(A[i]>A[i+1]){
            en = i+1;
            break;
        }
    }
    
    if(st==0 && en == n-1){
        return vector<int>(1,-1);
    }
    
    int min_el = *min_element(A.begin()+st,A.begin()+en+1);
    
    int max_el = *max_element(A.begin()+st,A.begin()+en+1);
    
    int bg =-1,sm=n;
    for(int i = 0;i<st;i++){
        if(A[i]>min_el){
            bg =i;
            break;
        }
    }
    
    for(int i = n-1;i>en;i--){
        if(A[i]<max_el){
            sm =i;
            break;
        }
    }
    
    if(bg!=-1){st=bg;}
    if(sm!=n){en = sm;}
    
    vector<int> op;op.push_back(st);op.push_back(en);

    
    return op;
}
