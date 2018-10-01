vector<int> Solution::solve(vector<int> &A, vector<int> &B) {
    
    priority_queue<pair<int,pair<int,int> >>p;
    
    sort(A.begin(),A.end());
    
    sort(B.begin(),B.end());
    
    int end1 = A.size()-1,end2 = B.size()-1;
    
    map<pair<int,int>,int> m;
    
    p.push({A[end1]+B[end2],{end1,end2}});
    
    m[{end1,end2}] = 1;
    
    int st = 1;
    
    vector<int> v;
    while(st<=A.size()){
        
        pair<int,pair<int,int>> tp = p.top();
        
        v.push_back(tp.first);
        
        p.pop();
        
        pair<int,int> p1 = {tp.second.first-1,tp.second.second};
        
        if(!m[p1]){
            p.push({A[p1.first]+B[p1.second],p1});
            m[p1] =1;
            
        };
        
        p1 = {tp.second.first,tp.second.second-1};
        
        if(!m[p1]){
            p.push({A[p1.first]+B[p1.second],p1});
            m[p1] =1;
        };
        
        
        st++;
        
        
        
        
        
    }
    
    
    return v;
    
    
}
