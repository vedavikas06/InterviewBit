//have to find the first missing positive number which is n't there in the array. so I seperated {positive},{0,negative} numbers seperately by swapping

int Solution::firstMissingPositive(vector<int> &A) {
    
    int n = A.size();
    
    int j = 0;
    
    for(int i =0;i<n;i++){
        if(A[i]<=0){
            swap(A[j],A[i]);
            j++;
        }
    }
    //start of +ve numbers
    int st = j;
    //maximum +ve numbers to be present
    int cnt = n-st;
    //Do the standard procedure (make the corresponding elements in corresponding array indexes -ve if they are present in this array.
    
    for(int i = st;i<n;i++){
        if(abs(A[i]) <= cnt){
            A[st+abs(A[i])-1] = -abs(A[st+abs(A[i])-1]);
        }
    }
    //check the missing element in the array by the above mentioned procedure. 
    for(int i = st;i<n;i++){
        if(A[i]>=0){
            return i-st+1;
        }
    }
    
    return cnt+1;
}
