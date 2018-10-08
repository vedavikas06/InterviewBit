vector<int> Solution::plusOne(vector<int> &A) {
    
    // int no = 0,in=1;
    
    // for(int i =A.size()-1;i>=0;i--){
    //     no = no+(in*A[i]);
    //     in*=10;
    // }
    
    // no+=1;
    
    // vector<int> v;
    
    // while(no){
    //     v.push_back(no%10);
    //     no/=10;
    // }
    
    
    // reverse(v.begin(),v.end());
    
    // return v;
    
    
    
    while(!A.empty() && A[0]==0){
        A.erase(A.begin());
       }
       
       if(A.empty()){
           A.push_back(0);
       }
    int n = A.size();
    if(A[n-1]<9){
        A[n-1]++;
        
        
    }else{
        int cr = (A[n-1]+1)/10;
        A[n-1] = (A[n-1]+1)%10;
        
        int x = n-2;
        while(cr && x>=0){
            
            int y = A[x];
            A[x] = (A[x]+cr)%10;
            cr = (y+cr)/10;
            x--;
            
        }
        
        if(cr==1){
            A.insert(A.begin(),1);
        }
        
        
    }
    
    
    
    return A;
    
}
