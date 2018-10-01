int Solution::solve(vector<vector<int> > &A) {
    
    int siz1 = A.size(),siz2 = A[0].size();
    
    for(int i= 0;i<siz1;i++){
        for(int j= 1;j<siz2;j++){
            if(A[i][j]){
                A[i][j]+=A[i][j-1];
            }
        }
    }
    
    for(int i= 0;i<siz1;i++){
     sort(A[i].begin(),A[i].end(),greater<int>());  
    }
    
    int maxx = INT_MIN;
    for(int i= 0;i<siz1;i++){
        for(int j= 0;j<siz2;j++){
           
          maxx = max(maxx,A[i][j]*(j+1));
            
        }
    }
    
    return maxx;
    
}
